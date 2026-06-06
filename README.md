# horb - Hydrogenic Orbital Radial Basis

Rust/CUDA Simulation Suite For Testing Hydrogen Electron Orbitals As Spiral Galaxy Dark Matter Halos

## Motivation

Standard computational astrophysics models spiral-galaxy dark-matter halos using smooth, parametric density curves (e.g., NFW, Burkert, Einasto). While these phenomenological profiles function as flexible data-fitting forms, they are rarely constrained by independent, non-spherical geometric observables. 

The **HORB** (`Hydrogenic Orbitals Radial Basis`) engine is built to test a radically different, morphology-first alternative: evaluating whether a scaled hydrogenic orbital density field centered on the Galactic core can simultaneously resolve macroscopic geometry (Fermi Bubbles) and dynamics (Flat Rotation Curves). Preliminary morphology indicates that the real angular probability density of the hydrogen $3d_{z^2}$ orbital:

$$\rho_{\Omega}(\theta) \propto \left(3\cos^2\theta - 1\right)^2$$

reproduces the distinctive bipolar lobe geometry, equatorial waist, and $54.74^\circ$ nodal cone opening structure of the Milky Way's Fermi bubbles. 

### The Fractal Connection
The broader theoretical framework driving this investigation is detailed in the accompanying independent ontology paper:

> **[The Fact of Fractal Tennis: The Universe As A Fractal Computer Defined by the Star=Electron+Atom=Galaxy Quine (arXiv/viXra:2605.0003)](https://ai.vixra.org/abs/2605.0003)**

In the TFOFT framework, the universe operates as a recursive, self-similar scale hierarchy. Under this identity, what macro-astrophysics labels as "exotic dark matter" is physically composed of sub-scale hydrogen gas clouds, which is sitting at the threshold of electron-generating micro-fusion at the Fermi Bubbles. 

This hypothesis yields three concrete, independently falsifiable predictions that this software is engineered to test:
1. **Morphology:** Does the $3d_{z^2}$ boundary template quantitatively match the observed Fermi bubble outlines? Initial analysis says YES!
2. **Dynamics:** Does the morphology-constrained density field yield a competitive galactic circular-velocity curve ($v_c$) contribution without over-parameterization?
3. **Substrate Density:** Does the internal X-ray surface brightness distribution within the Fermi lobes trace the predicted $3d_{z^2}$ internal density matrix?

**HORB** does not assume the validity of the TFOFT ontology *a priori*. Instead, it provides the high-performance Rust/CUDA pipeline required to subject its core macro-scale predictions to rigorous algorithmic stress-testing, benchmarking the results directly against legacy astrophysics profiles.

## License

This project is dual-licensed under either:

* Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
* MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.
