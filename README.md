# horb - Hydrogenic Orbital Radial Basis

Rust/CUDA Simulation Suite For Testing Hydrogen Electron Orbitals As Spiral Galaxy Dark Matter Halos

## Current Best Result

The current best Milky Way candidate is a morphology-constrained hydrogenic `3d_z2` halo tested against the Sofue 2020 unified Milky Way rotation curve using a literature bulge+disk baryonic model and a CUDA direct-summation disk-plane force calculation.

| Parameter | Value |
|---|---:|
| Orbital state | `3d_z2` |
| Orbital scale `a0_star` | `1.0 kpc` |
| HORB dark halo mass | `2.5e11 M_sun` |
| HORB force model | CUDA disk-plane direct summation |
| Grid size | `96^3` mass cells |
| Grid extent | `±80 kpc` |
| Softening length | `0.25 kpc` |
| Baryon model | Sofue literature bulge+disk |
| Stellar disk mass | `3.41e10 M_sun` |
| Stellar disk scale | `3.19 kpc` |
| Bulge mass | `1.652e10 M_sun` |
| Bulge scale | `0.522 kpc` |

Using the same literature baryonic model for every halo, the CUDA HORB candidate is compared directly against the previous spherical HORB approximation and against NFW, pseudo-isothermal, and Burkert halo baselines, with all curves plotted against the Sofue 2020 unified Milky Way rotation curve.

![CUDA HORB disk-plane force candidate versus classical halo models and the Milky Way rotation curve](images/mw_cuda_horb_diskplane_total_vs_sofue2020.png)

This benchmark replaces the earlier spherical-enclosed-mass approximation. The CUDA calculation evaluates the actual disk-plane gravitational response of the non-spherical `3d_z2` orbital density field. The result lowers the mid-radius HORB peak relative to the spherical approximation, which is exactly the region where the spherical model overshot the observed Milky Way rotation curve.

The plot is generated reproducibly by:

```bash
LD_LIBRARY_PATH="$PWD/cuda_kernels/cuda:${LD_LIBRARY_PATH:-}" \
cargo run -p cuda_kernels --bin test_horb_disk_plane_curve -- \
  1.0 2.5e11 96 80 0.25

cargo run -q -p curve_fitter -- compare_dm_fixed_baselines \
  3d_z2 1.0 2.5e11 2.5e11 \
  5.0 15.0 8.0 80.0 \
  > curves/spherical_horb_compare_a1.0_m2.5e11.csv

./scripts/plot_cuda_total_vs_mw.py \
  --cuda-csv curves/cuda_horb_diskplane_a1.000_m2.500e11_n96.csv \
  --dm-compare-csv curves/spherical_horb_compare_a1.0_m2.5e11.csv \
  --rc-csv data/milky_way/sofue_2020_standard_rc.csv \
  --baryons-csv data/milky_way/sofue_literature_baryons.csv \
  -o plots/mw_cuda_horb_diskplane_total_vs_sofue2020.png
```


## Motivation

Standard computational astrophysics models spiral-galaxy dark-matter halos using smooth, parametric density curves (e.g., NFW, Burkert, Einasto). While these phenomenological profiles function as flexible data-fitting forms, they are rarely constrained by independent, non-spherical geometric observables. 

The **HORB** (`Hydrogenic Orbitals Radial Basis`) engine is built to test a radically different, morphology-first alternative: evaluating whether a scaled hydrogenic orbital density field centered on the Galactic core can simultaneously resolve macroscopic geometry (Fermi Bubbles) and dynamics (Flat Rotation Curves). Preliminary morphology indicates that the real angular probability density of the hydrogen $3d_{z^2}$ orbital:

$$\rho_{\Omega}(\theta) \propto \left(3\cos^2\theta - 1\right)^2$$

reproduces the distinctive bipolar lobe geometry, equatorial waist, and $54.74^\circ$ nodal cone opening structure of the Milky Way's Fermi bubbles. 

### The Fractal Connection

- **Step 0:** Modern philosiphy takes Physics models and derives ontology, a system built by complete morons. We need the correct ontology first.
  
- **Step 1:** The final fundamental laws of physics must be self-contained and cannot rely on external “flying spaghetti monster” concepts. This recursion forces dynamic fractals to be the only viable geometric ontology.

- **Step 2:** Classical physics is built from averages and should fail at fractal scale boundaries.

- **Step 3:** The clearest boundary objects where classical physics fails are atoms with quantized energy and galaxies with flat rotation curves.

- **Step 4:** Assume Galaxy → Atom to define the scale factor.

- **Step 5:** Apply that independently defined scale factor to the Sun’s Schwarzschild–de Sitter balance radius:
  - Solar balance radius: $r_{\mathrm{SdS},\odot}=\left(\frac{3GM_\odot}{\Lambda c^2}\right)^{1/3}\approx3.42\times10^{18}\,\mathrm{m}$
  - Scaled prediction: $\frac{r_{\mathrm{SdS},\odot}}{S_{G\to A}}\approx3.83\times10^{-13}\,\mathrm{m}$
  - Electron reduced Compton wavelength: $\bar{\lambda}_e=\frac{\hbar}{m_ec}\approx3.86\times10^{-13}\,\mathrm{m}$
  - Therefore: **Star ↔ Electron**

- **Step 6:** Apply extremal three-charge microstate counting to Sagittarius A*:
  - Black-hole entropy: $S_{\mathrm{BH}}=\frac{4\pi GM_{\mathrm{SgrA^*}}^2}{\hbar c}\approx1.94\times10^{90}$
  - Three-charge entropy: $S_{\mathrm{BH}}=2\pi\sqrt{Q_1Q_5N}$
  - Equal-channel depth: $\ln Q_i=\frac{1}{3}\ln(Q_1Q_5N)\approx137.37$
  - Quantized register depth: $\left\lfloor\ln Q_i\right\rfloor=137$
  - This gives a finite, quantized crossing condition and the photonic mechanism.

- **Step 7:** The fractal continues beyond atoms and galaxies; the next lower scales appears as dark matter.

- **Step 8:** HORB tests the claim by fitting galactic rotation curves with hydrogenic orbital structure. And they beat every PhD-developed model.



## CUDA Build Note

The CUDA shared library `libhorb_cuda.so` is intentionally ignored by git because it is a local build artifact. On a fresh clone, build it before running CUDA tests:

```bash
cd cuda_kernels/cuda
./build.sh
cd ../..
LD_LIBRARY_PATH="$PWD/cuda_kernels/cuda:${LD_LIBRARY_PATH:-}" cargo test
```

## License

This project is dual-licensed under either:

* Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
* MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.
