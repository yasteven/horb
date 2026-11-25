# Quantised inertia in the README comparisons

The curves use McCulloch (2017), *Galaxy rotations from quantised inertia and
visible matter only*, [DOI 10.1007/s10509-017-3128-6](https://doi.org/10.1007/s10509-017-3128-6).
Equation 8 gives the total acceleration from the fixed baryonic acceleration:

\[
g_{\rm QI}=\sqrt{g_b(g_b+a_{\rm QI})},\qquad
v_{\rm QI}=\sqrt{r g_{\rm QI}},\qquad
a_{\rm QI}=2c^2/\Theta.
\]

Use c=299792458 m/s and the paper's present cosmic diameter
Theta(0)=8.8e26 m. MW and M31 use the present epoch. J2141 uses its cosmological
redshift z=0.138 and the paper's Equation 9, Theta(z)=Theta(0)/(1+z).
This redshift scaling is the paper's stated prescription.

The input baryons, distances and geometry stay fixed. No QI parameter is fitted.
The result is the total speed; adding the baryons again would double count them.
Negative total inward baryonic force is outside this prescription and is rejected.
These models have nonnegative baryonic force; J2141's zero-radius grid node uses
the continuous central zero-speed limit of its cored baryons.

For J2141, compute intrinsic circular speed first, then apply exactly the same
inclination, seeing, slit and spatial-bin integration as the other models.
The plotted RMS is computed on the six outer LOS measurements after the refined
observing integration. The shaded inner measurements are displayed but excluded
from every score. QI does not receive a lensing prediction from this rotation law.

`./run.sh --release do readme-highlights` rebuilds all three rotation comparisons
and their companion figures. Per-point predictions, scores and constants are
tracked in `data/highlights/`. The single orbital remains implemented for legacy
calculations but is omitted from the displayed README comparisons.
