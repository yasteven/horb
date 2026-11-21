#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>

static constexpr double G_KPC = 4.3009e-6;

// State IDs
// 0  = 1s
// 10 = 2p_z
// 11 = 2p_x
// 12 = 2p_y
// 20 = 3d_z2
// 21 = 3d_x2_y2
// 22 = 3d_xy
// 23 = 3d_xz
// 24 = 3d_yz

struct MassCell {
    double x;
    double y;
    double z;
    double m;
};

struct OrbitalTerm {
    int state_id;
    double coeff;
};

struct EulerRotation {
    double r00, r01, r02;
    double r10, r11, r12;
    double r20, r21, r22;
};

__device__ __forceinline__
void rotate_point(
    const EulerRotation rot,
    double x,
    double y,
    double z,
    double* xr,
    double* yr,
    double* zr
) {
    *xr = rot.r00 * x + rot.r01 * y + rot.r02 * z;
    *yr = rot.r10 * x + rot.r11 * y + rot.r12 * z;
    *zr = rot.r20 * x + rot.r21 * y + rot.r22 * z;
}

__device__ __forceinline__
double orbital_amplitude(
    int state_id,
    double x,
    double y,
    double z,
    double a0
) {
    double r2 = x*x + y*y + z*z;
    double r = sqrt(r2);

    if (r <= 1.0e-14) {
        if (state_id == 0) {
            return 1.0;
        }
        return 0.0;
    }

    double inv_r = 1.0 / r;
    double inv_r2 = 1.0 / r2;

    if (state_id == 0) {
        return exp(-r / a0);
    }

    if (state_id == 10) {
        return r * exp(-r / (2.0 * a0)) * (z * inv_r);
    }

    if (state_id == 11) {
        return r * exp(-r / (2.0 * a0)) * (x * inv_r);
    }

    if (state_id == 12) {
        return r * exp(-r / (2.0 * a0)) * (y * inv_r);
    }

    double radial_3d = r2 * exp(-r / (3.0 * a0));

    if (state_id == 20) {
        double c = z * inv_r;
        return radial_3d * (3.0 * c * c - 1.0);
    }

    if (state_id == 21) {
        return radial_3d * ((x*x - y*y) * inv_r2);
    }

    if (state_id == 22) {
        return radial_3d * ((x*y) * inv_r2);
    }

    if (state_id == 23) {
        return radial_3d * ((x*z) * inv_r2);
    }

    if (state_id == 24) {
        return radial_3d * ((y*z) * inv_r2);
    }

    return 0.0;
}

__device__ __forceinline__
double orbital_density_single(
    int state_id,
    double x,
    double y,
    double z,
    double a0,
    EulerRotation rot
) {
    double xr, yr, zr;
    rotate_point(rot, x, y, z, &xr, &yr, &zr);

    double psi = orbital_amplitude(state_id, xr, yr, zr, a0);
    return psi * psi;
}

__device__ __forceinline__
double orbital_density_superposition(
    const OrbitalTerm* terms,
    int term_count,
    double x,
    double y,
    double z,
    double a0,
    EulerRotation rot
) {
    double xr, yr, zr;
    rotate_point(rot, x, y, z, &xr, &yr, &zr);

    double psi = 0.0;

    for (int i = 0; i < term_count; i++) {
        psi += terms[i].coeff * orbital_amplitude(terms[i].state_id, xr, yr, zr, a0);
    }

    return psi * psi;
}

__global__
void fill_single_orbital_density_kernel(
    double* rho,
    int n_side,
    double extent_kpc,
    int state_id,
    double a0,
    EulerRotation rot
) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    size_t n_total = (size_t)n_side * (size_t)n_side * (size_t)n_side;

    if (idx >= n_total) return;

    int ix = idx % n_side;
    int iy = (idx / n_side) % n_side;
    int iz = idx / ((size_t)n_side * (size_t)n_side);

    double dx = 2.0 * extent_kpc / (double)n_side;

    double x = -extent_kpc + ((double)ix + 0.5) * dx;
    double y = -extent_kpc + ((double)iy + 0.5) * dx;
    double z = -extent_kpc + ((double)iz + 0.5) * dx;

    rho[idx] = orbital_density_single(state_id, x, y, z, a0, rot);
}

__global__
void fill_superposition_density_kernel(
    double* rho,
    int n_side,
    double extent_kpc,
    const OrbitalTerm* terms,
    int term_count,
    double a0,
    EulerRotation rot
) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    size_t n_total = (size_t)n_side * (size_t)n_side * (size_t)n_side;

    if (idx >= n_total) return;

    int ix = idx % n_side;
    int iy = (idx / n_side) % n_side;
    int iz = idx / ((size_t)n_side * (size_t)n_side);

    double dx = 2.0 * extent_kpc / (double)n_side;

    double x = -extent_kpc + ((double)ix + 0.5) * dx;
    double y = -extent_kpc + ((double)iy + 0.5) * dx;
    double z = -extent_kpc + ((double)iz + 0.5) * dx;

    rho[idx] = orbital_density_superposition(terms, term_count, x, y, z, a0, rot);
}

__global__
void reduce_sum_kernel(
    const double* values,
    double* partial,
    size_t n
) {
    extern __shared__ double shared[];

    size_t tid = threadIdx.x;
    size_t global = blockIdx.x * blockDim.x + threadIdx.x;
    size_t stride = blockDim.x * gridDim.x;

    double sum = 0.0;

    for (size_t i = global; i < n; i += stride) {
        sum += values[i];
    }

    shared[tid] = sum;
    __syncthreads();

    for (size_t s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            shared[tid] += shared[tid + s];
        }
        __syncthreads();
    }

    if (tid == 0) {
        partial[blockIdx.x] = shared[0];
    }
}

__global__
void materialize_mass_cells_kernel(
    const double* rho,
    MassCell* cells,
    int n_side,
    double extent_kpc,
    double mass_scale
) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    size_t n_total = (size_t)n_side * (size_t)n_side * (size_t)n_side;

    if (idx >= n_total) return;

    int ix = idx % n_side;
    int iy = (idx / n_side) % n_side;
    int iz = idx / ((size_t)n_side * (size_t)n_side);

    double dx = 2.0 * extent_kpc / (double)n_side;

    double x = -extent_kpc + ((double)ix + 0.5) * dx;
    double y = -extent_kpc + ((double)iy + 0.5) * dx;
    double z = -extent_kpc + ((double)iz + 0.5) * dx;

    cells[idx].x = x;
    cells[idx].y = y;
    cells[idx].z = z;
    cells[idx].m = rho[idx] * mass_scale;
}

__global__
void disk_plane_curve_tiled_kernel(
    const MassCell* cells,
    size_t cell_count,
    const double* radii,
    double* v_out,
    size_t radius_count,
    double softening_kpc
) {
    int radius_idx = blockIdx.x;
    if ((size_t)radius_idx >= radius_count) return;

    int tid = threadIdx.x;
    int threads = blockDim.x;

    extern __shared__ double shared_sum[];

    double R = radii[radius_idx];
    double eps2 = softening_kpc * softening_kpc;

    double local_a_in = 0.0;

    // Correct cooperative reduction:
    // each thread handles a disjoint subset of cells.
    for (size_t j = (size_t)tid; j < cell_count; j += (size_t)threads) {
        double dx = cells[j].x - R;
        double dy = cells[j].y;
        double dz = cells[j].z;

        double d2 = dx*dx + dy*dy + dz*dz + eps2;
        double d = sqrt(d2);
        double d3 = d2 * d;

        double ax = G_KPC * cells[j].m * dx / d3;

        // inward radial acceleration at test point (R,0,0)
        local_a_in += -ax;
    }

    shared_sum[tid] = local_a_in;
    __syncthreads();

    for (int stride = threads / 2; stride > 0; stride >>= 1) {
        if (tid < stride) {
            shared_sum[tid] += shared_sum[tid + stride];
        }
        __syncthreads();
    }

    if (tid == 0) {
        double a_in = shared_sum[0];
        double v2 = R * a_in;
        v_out[radius_idx] = (v2 > 0.0) ? sqrt(v2) : 0.0;
    }
}


static EulerRotation identity_rotation() {
    EulerRotation r;
    r.r00 = 1.0; r.r01 = 0.0; r.r02 = 0.0;
    r.r10 = 0.0; r.r11 = 1.0; r.r12 = 0.0;
    r.r20 = 0.0; r.r21 = 0.0; r.r22 = 1.0;
    return r;
}

static bool check_cuda(cudaError_t err, const char* label) {
    if (err != cudaSuccess) {
        fprintf(stderr, "%s failed: %s\n", label, cudaGetErrorString(err));
        return false;
    }
    return true;
}

static double device_sum_double(const double* d_values, size_t n) {
    int threads = 256;
    int blocks = 1024;

    double* d_partial = nullptr;
    cudaMalloc(&d_partial, blocks * sizeof(double));

    reduce_sum_kernel<<<blocks, threads, threads * sizeof(double)>>>(d_values, d_partial, n);
    cudaDeviceSynchronize();

    double* h_partial = (double*)malloc(blocks * sizeof(double));
    cudaMemcpy(h_partial, d_partial, blocks * sizeof(double), cudaMemcpyDeviceToHost);

    double sum = 0.0;
    for (int i = 0; i < blocks; i++) {
        sum += h_partial[i];
    }

    free(h_partial);
    cudaFree(d_partial);

    return sum;
}

extern "C" void compute_single_orbital_disk_curve_cuda(
    int state_id,
    double a0,
    double total_mass_msun,
    int n_side,
    double extent_kpc,
    double softening_kpc,
    const double* host_radii,
    double* host_v_out,
    size_t radius_count,
    EulerRotation rot
) {
    size_t cell_count = (size_t)n_side * (size_t)n_side * (size_t)n_side;

    double* d_rho = nullptr;
    MassCell* d_cells = nullptr;
    double* d_radii = nullptr;
    double* d_v_out = nullptr;

    double rho_sum = 0.0;
    double mass_scale = 0.0;

    cudaMalloc(&d_rho, cell_count * sizeof(double));
    cudaMalloc(&d_cells, cell_count * sizeof(MassCell));
    cudaMalloc(&d_radii, radius_count * sizeof(double));
    cudaMalloc(&d_v_out, radius_count * sizeof(double));

    int threads = 256;
    int blocks = (int)((cell_count + threads - 1) / threads);

    fill_single_orbital_density_kernel<<<blocks, threads>>>(
        d_rho,
        n_side,
        extent_kpc,
        state_id,
        a0,
        rot
    );

    if (!check_cuda(cudaDeviceSynchronize(), "fill_single_orbital_density_kernel")) goto cleanup;

    rho_sum = device_sum_double(d_rho, cell_count);

    if (rho_sum <= 0.0) {
        fprintf(stderr, "rho_sum non-positive\n");
        goto cleanup;
    }

    mass_scale = total_mass_msun / rho_sum;

    materialize_mass_cells_kernel<<<blocks, threads>>>(
        d_rho,
        d_cells,
        n_side,
        extent_kpc,
        mass_scale
    );

    if (!check_cuda(cudaDeviceSynchronize(), "materialize_mass_cells_kernel")) goto cleanup;

    cudaMemcpy(d_radii, host_radii, radius_count * sizeof(double), cudaMemcpyHostToDevice);

    {
        int curve_blocks = (int)radius_count;
        int curve_threads = 256;
        size_t shared_bytes = curve_threads * sizeof(double);

        disk_plane_curve_tiled_kernel<<<curve_blocks, curve_threads, shared_bytes>>>(
            d_cells,
            cell_count,
            d_radii,
            d_v_out,
            radius_count,
            softening_kpc
        );
    }

    if (!check_cuda(cudaDeviceSynchronize(), "disk_plane_curve_tiled_kernel")) goto cleanup;

    cudaMemcpy(host_v_out, d_v_out, radius_count * sizeof(double), cudaMemcpyDeviceToHost);

cleanup:
    cudaFree(d_rho);
    cudaFree(d_cells);
    cudaFree(d_radii);
    cudaFree(d_v_out);
}

extern "C" void compute_superposition_disk_curve_cuda(
    const OrbitalTerm* host_terms,
    int term_count,
    double a0,
    double total_mass_msun,
    int n_side,
    double extent_kpc,
    double softening_kpc,
    const double* host_radii,
    double* host_v_out,
    size_t radius_count,
    EulerRotation rot
) {
    size_t cell_count = (size_t)n_side * (size_t)n_side * (size_t)n_side;

    double* d_rho = nullptr;
    MassCell* d_cells = nullptr;
    OrbitalTerm* d_terms = nullptr;
    double* d_radii = nullptr;
    double* d_v_out = nullptr;

    double rho_sum = 0.0;
    double mass_scale = 0.0;

    cudaMalloc(&d_rho, cell_count * sizeof(double));
    cudaMalloc(&d_cells, cell_count * sizeof(MassCell));
    cudaMalloc(&d_terms, term_count * sizeof(OrbitalTerm));
    cudaMalloc(&d_radii, radius_count * sizeof(double));
    cudaMalloc(&d_v_out, radius_count * sizeof(double));

    cudaMemcpy(d_terms, host_terms, term_count * sizeof(OrbitalTerm), cudaMemcpyHostToDevice);

    int threads = 256;
    int blocks = (int)((cell_count + threads - 1) / threads);

    fill_superposition_density_kernel<<<blocks, threads>>>(
        d_rho,
        n_side,
        extent_kpc,
        d_terms,
        term_count,
        a0,
        rot
    );

    if (!check_cuda(cudaDeviceSynchronize(), "fill_superposition_density_kernel")) goto cleanup;

    rho_sum = device_sum_double(d_rho, cell_count);

    if (rho_sum <= 0.0) {
        fprintf(stderr, "rho_sum non-positive\n");
        goto cleanup;
    }

    mass_scale = total_mass_msun / rho_sum;

    materialize_mass_cells_kernel<<<blocks, threads>>>(
        d_rho,
        d_cells,
        n_side,
        extent_kpc,
        mass_scale
    );

    if (!check_cuda(cudaDeviceSynchronize(), "materialize_mass_cells_kernel")) goto cleanup;

    cudaMemcpy(d_radii, host_radii, radius_count * sizeof(double), cudaMemcpyHostToDevice);

    {
        int curve_blocks = (int)radius_count;
        int curve_threads = 256;
        size_t shared_bytes = curve_threads * sizeof(double);

        disk_plane_curve_tiled_kernel<<<curve_blocks, curve_threads, shared_bytes>>>(
            d_cells,
            cell_count,
            d_radii,
            d_v_out,
            radius_count,
            softening_kpc
        );
    }

    if (!check_cuda(cudaDeviceSynchronize(), "disk_plane_curve_tiled_kernel")) goto cleanup;

    cudaMemcpy(host_v_out, d_v_out, radius_count * sizeof(double), cudaMemcpyDeviceToHost);

cleanup:
    cudaFree(d_rho);
    cudaFree(d_cells);
    cudaFree(d_terms);
    cudaFree(d_radii);
    cudaFree(d_v_out);
}

extern "C" EulerRotation make_identity_rotation_cuda() {
    return identity_rotation();
}

extern "C" void compute_disk_plane_rotation_curve(
    const MassCell* host_cells,
    size_t cell_count,
    const double* host_radii,
    double* host_v_out,
    size_t radius_count,
    double softening_kpc
) {
    MassCell* d_cells = nullptr;
    double* d_radii = nullptr;
    double* d_v_out = nullptr;

    cudaMalloc(&d_cells, cell_count * sizeof(MassCell));
    cudaMalloc(&d_radii, radius_count * sizeof(double));
    cudaMalloc(&d_v_out, radius_count * sizeof(double));

    cudaMemcpy(d_cells, host_cells, cell_count * sizeof(MassCell), cudaMemcpyHostToDevice);
    cudaMemcpy(d_radii, host_radii, radius_count * sizeof(double), cudaMemcpyHostToDevice);

    int curve_blocks = (int)radius_count;
    int curve_threads = 256;
    size_t shared_bytes = curve_threads * sizeof(double);

    disk_plane_curve_tiled_kernel<<<curve_blocks, curve_threads, shared_bytes>>>(
        d_cells,
        cell_count,
        d_radii,
        d_v_out,
        radius_count,
        softening_kpc
    );

    if (!check_cuda(cudaDeviceSynchronize(), "compute_disk_plane_rotation_curve wrapper")) {
        goto cleanup;
    }

    cudaMemcpy(host_v_out, d_v_out, radius_count * sizeof(double), cudaMemcpyDeviceToHost);

cleanup:
    cudaFree(d_cells);
    cudaFree(d_radii);
    cudaFree(d_v_out);
}

// Signed density fits must preserve R*a_in, including negative values. The
// historical speed API above deliberately keeps its old behavior.
__global__ void disk_plane_signed_v2_kernel(
    const MassCell* cells, size_t count, const double* radii, double* out,
    size_t radius_count, double softening
) {
    size_t i = blockIdx.x;
    if (i >= radius_count) return;
    extern __shared__ double sums[];
    double r = radii[i];
    double total = 0.0;
    for (size_t j = threadIdx.x; j < count; j += blockDim.x) {
        double dx = r - cells[j].x;
        double d2 = dx*dx + cells[j].y*cells[j].y + cells[j].z*cells[j].z + softening*softening;
        total += G_KPC * r * cells[j].m * dx / (d2*sqrt(d2));
    }
    sums[threadIdx.x] = total;
    __syncthreads();
    for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
        if (threadIdx.x < stride) sums[threadIdx.x] += sums[threadIdx.x + stride];
        __syncthreads();
    }
    if (threadIdx.x == 0) out[i] = sums[0];
}

extern "C" int compute_disk_plane_signed_v2_cuda(
    const MassCell* cells, size_t count, const double* radii, double* out,
    size_t radius_count, double softening
) {
    if (!cells || !radii || !out || count == 0 || radius_count == 0 ||
        radius_count > 100001 || !std::isfinite(softening) || softening < 0.0) return 0;
    MassCell* dc = nullptr;
    double* dr = nullptr;
    double* dk = nullptr;
    bool ok = check_cuda(cudaMalloc(&dc, count*sizeof(MassCell)), "signed cells allocation");
    if (ok) ok = check_cuda(cudaMalloc(&dr, radius_count*sizeof(double)), "signed radii allocation");
    if (ok) ok = check_cuda(cudaMalloc(&dk, radius_count*sizeof(double)), "signed force allocation");
    if (ok) ok = check_cuda(cudaMemcpy(dc,cells,count*sizeof(MassCell),cudaMemcpyHostToDevice), "signed cells copy");
    if (ok) ok = check_cuda(cudaMemcpy(dr,radii,radius_count*sizeof(double),cudaMemcpyHostToDevice), "signed radii copy");
    if (ok) {
        disk_plane_signed_v2_kernel<<<(unsigned int)radius_count,256,256*sizeof(double)>>>(dc,count,dr,dk,radius_count,softening);
        ok = check_cuda(cudaGetLastError(), "signed force launch");
    }
    if (ok) ok = check_cuda(cudaDeviceSynchronize(), "signed force execution");
    if (ok) ok = check_cuda(cudaMemcpy(out,dk,radius_count*sizeof(double),cudaMemcpyDeviceToHost), "signed force result");
    cudaFree(dc); cudaFree(dr); cudaFree(dk);
    if (ok) for (size_t i=0; i<radius_count; ++i) if (!std::isfinite(out[i])) return 0;
    return ok ? 1 : 0;
}
