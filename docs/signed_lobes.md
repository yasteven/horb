# Two HORB density definitions

The comparison suite retains the original basis and fits a second, independent
family. Both sum signed gravitational forces before taking a circular speed.

| Model | Density basis | Reference mass |
|---|---|---|
| `squared` | \(b_i=\lvert\psi_i\rvert^2\) | \(\int b_i\,d^3x=1\) |
| `psi` | \(b_i=\psi_i/\int \lvert\psi_i\rvert\,d^3x\) | \(\int \lvert b_i\rvert\,d^3x=1\) |

In either family, \(\rho=\sum_i w_i M_{\rm ref,i} b_i\). The second definition
preserves **both angular-lobe signs and radial-node signs**. It uses the real
harmonics and positive radial normalization already defined in the Rust engine.
There is no absolute value or squaring of the gravitational source. The absolute
value appears only in the mass normalization and concentration statistic.

This defines a proposed signed gravitational density, not a conventional quantum
probability. The coherent quantum expression \(|\sum c_i\psi_i|^2\) is another
model: it has interference terms and remains nonnegative. Neither HORB family
implements time evolution, transitions, relaxation or a galaxy-scale energy law.
A global phase/sign reversal of a real basis is absorbed by reversing its weight;
it is not a separate physical orbital in the dictionary.

## Force and identifiable modes

The disk is in the x–y plane, with the orbital axis along z. For real `psi`, the
axisymmetric source has one multipole, L=ell, whereas squared orbitals need
L=0,2,...,2ell. The force includes **interior and exterior** radial integrals:

\[
v_i^2(R)=\frac{4\pi G}{2\ell+1}P_\ell(0)
\left[\frac{\ell+1}{R^{\ell+1}}\int_0^R\rho_\ell(r)r^{\ell+2}dr
-\ell R^\ell\int_R^\infty\rho_\ell(r)r^{1-\ell}dr\right].
\]

The exported response is signed; no per-component square root or clipping occurs.
Only the **total** baryon-plus-HORB force must support the fitted circular motion.
Positivity is checked on a dense recorded radial grid, not proved globally.

For signed psi, odd ell has zero equatorial radial response and m!=0 has zero
azimuthal mean. Such terms cannot be determined from these circular-speed data.
They are rejected in the rotation basis export instead of assigning arbitrary
weights. Odd modes may exert vertical forces; this fit does not establish their
vertical equilibrium. The default identifiable psi dictionary contains
`1s,2s,3s,3d_z2` at five scales (20 columns). The original squared dictionary
retains its six states (30 columns). Both searches allow up to four weights.
Explicit `n:l:0` states through n=10, ell<=2 are supported; larger n is not by
itself a calibrated jet energy. The full 3D real harmonics remain available in
Rust density evaluation, including orientations that rotation alone cannot fit.

## Mass and concentration

Every ell>0 psi basis has exactly zero net mass. Its positive and negative lobes
still source a potential. The normalization therefore never divides by net mass.

- `amplitude_mass_msun` = weight × reference absolute mass; retains coefficient sign.
- `absolute_component_mass_msun` = absolute value of that amplitude mass.
- `signed_mass_msun` = integrated **net** component mass, which can be zero.
- `net_mass_fraction` = net component mass / signed amplitude mass.
- Positive/negative lobe masses are recorded before overlap between components.

Concentration is \(\max_i|w_iM_i|/\sum_i|w_iM_i|\), using absolute component
strength. It is not a probability, confidence level, or fraction of radial force.
The AGN export uses this definition for both models, so a zero-monopole d orbital
does not disappear from the statistics. Negative mass fractions count negative
lobes, including those within a positive-coefficient psi orbital.

The empty subset is a valid candidate: all coefficients are zero, dark density
is zero, and the curve is exactly the fixed baryonic prediction. Concentration
is undefined for this solution. A Keplerian \(v\propto R^{-1/2}\) follows only
outside a finite enclosed baryonic mass, not throughout an extended disk.

## Commands and records

```sh
./run.sh --release do signed-fit                 # compare both definitions, MW
./run.sh --release do m31-signed-fit             # same for M31
./run.sh --release do signed-fit --density-model psi
./run.sh --release do sparc-fit --name psi_models --density-model both --jobs 2
./run.sh --release do sparc-fit --name higher_n --states 1s,2s,3d_z2,5:2:0,6:2:0
./run.sh do swells-fit                          # both rotation fits, then lensing
./run.sh --release do agc114905-fit
./run.sh --release do readme-highlights
./run.sh data agn --fits-dir reports/sparc/psi_models --density-model psi --out-dir reports/agn_psi
```

`--weights nonnegative` constrains coefficients. It ensures positive density
for the squared basis, but cannot remove the intrinsic negative lobes of psi.
SWELLS `--nonnegative` controls its squared-density comparison; the independent
signed-psi fit always allows signed coefficients.

Squared CSV schema 1 remains readable. Psi force libraries use schema 2 with
`density_model=psi`, `normalization=absolute_mass_full_space`, and net mass
fractions. Mixed or mislabeled libraries are rejected. New weight CSVs identify
their model explicitly; old weight CSVs retain their squared interpretation.

The CUDA exporter accepts the new signed mass cells; the same signed direct-force
kernel is used. CPU multipoles, Rust/SciPy agreement, node signs, normalization,
direct-force checks and inclined lens projection are tested. CUDA execution must
be checked on a CUDA-equipped host; it is not inferred from the CPU tests.

SWELLS projection rotates the density into the observed galaxy inclination and
retains signed surface density. Its lensing calculation assumes the same GR
weak-field thin-lens source law as the earlier model; this is an additional
assumption, not a derived relativistic theory for signed mass. Rotation-selected
and closest-lensing candidates remain labeled separately.

The signed-ψ activity export also works directly from the bundled 175-galaxy
fit records: `./run.sh data agn --density-model psi --out-dir reports/agn_psi`.
Use `--fits-dir` to select a different completed run.
