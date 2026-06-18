#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>

struct MassCell {
    double x;
    double y;
    double z;
    double m;
};

static constexpr double G_KPC = 4.3009e-6; // kpc (km/s)^2 Msun^-1

__global__
void disk_plane_curve_kernel(
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
    int block_size = blockDim.x;

    extern __shared__ double shared_sum[];

    double R = radii[radius_idx];
    double eps2 = softening_kpc * softening_kpc;

    // Test particle at (R, 0, 0), disk plane.
    // Acceleration from cell:
    //   a_x = G m (x_cell - R) / d^3
    // inward radial acceleration at +x is:
    //   a_in = -a_x
    double local_a_in = 0.0;

    for (size_t i = tid; i < cell_count; i += block_size) {
        double dx = cells[i].x - R;
        double dy = cells[i].y;
        double dz = cells[i].z;

        double d2 = dx * dx + dy * dy + dz * dz + eps2;
        double d = sqrt(d2);
        double d3 = d2 * d;

        double ax = G_KPC * cells[i].m * dx / d3;
        local_a_in += -ax;
    }

    shared_sum[tid] = local_a_in;
    __syncthreads();

    for (int stride = block_size / 2; stride > 0; stride >>= 1) {
        if (tid < stride) {
            shared_sum[tid] += shared_sum[tid + stride];
        }
        __syncthreads();
    }

    if (tid == 0) {
        double a_in = shared_sum[0];
        double v2 = R * a_in;

        if (v2 > 0.0) {
            v_out[radius_idx] = sqrt(v2);
        } else {
            v_out[radius_idx] = 0.0;
        }
    }
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

    cudaError_t err;

    err = cudaMalloc(&d_cells, cell_count * sizeof(MassCell));
    if (err != cudaSuccess) {
        fprintf(stderr, "cudaMalloc d_cells failed: %s\n", cudaGetErrorString(err));
        return;
    }

    err = cudaMalloc(&d_radii, radius_count * sizeof(double));
    if (err != cudaSuccess) {
        fprintf(stderr, "cudaMalloc d_radii failed: %s\n", cudaGetErrorString(err));
        cudaFree(d_cells);
        return;
    }

    err = cudaMalloc(&d_v_out, radius_count * sizeof(double));
    if (err != cudaSuccess) {
        fprintf(stderr, "cudaMalloc d_v_out failed: %s\n", cudaGetErrorString(err));
        cudaFree(d_cells);
        cudaFree(d_radii);
        return;
    }

    err = cudaMemcpy(d_cells, host_cells, cell_count * sizeof(MassCell), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy cells failed: %s\n", cudaGetErrorString(err));
        cudaFree(d_cells);
        cudaFree(d_radii);
        cudaFree(d_v_out);
        return;
    }

    err = cudaMemcpy(d_radii, host_radii, radius_count * sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy radii failed: %s\n", cudaGetErrorString(err));
        cudaFree(d_cells);
        cudaFree(d_radii);
        cudaFree(d_v_out);
        return;
    }

    int threads = 256;
    int blocks = (int)radius_count;
    size_t shared_bytes = threads * sizeof(double);

    disk_plane_curve_kernel<<<blocks, threads, shared_bytes>>>(
        d_cells,
        cell_count,
        d_radii,
        d_v_out,
        radius_count,
        softening_kpc
    );

    err = cudaDeviceSynchronize();
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA kernel failed: %s\n", cudaGetErrorString(err));
        cudaFree(d_cells);
        cudaFree(d_radii);
        cudaFree(d_v_out);
        return;
    }

    err = cudaMemcpy(host_v_out, d_v_out, radius_count * sizeof(double), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy v_out failed: %s\n", cudaGetErrorString(err));
    }

    cudaFree(d_cells);
    cudaFree(d_radii);
    cudaFree(d_v_out);
}
