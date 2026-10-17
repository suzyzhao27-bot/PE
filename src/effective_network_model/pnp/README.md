# Finite PNP with prescribed shared ion-entry energy

This local extension copies the v3 finite PNP core without changing originals.
`solver_v3_baseline.py` preserves the exact input core; `legacy_tests.py` is the
unaltered v3 test file, importing the extended solver with zero energy.

Use existing `C:\ProgramData\Miniconda3\python.exe`; NumPy 2.2.6 and SciPy 1.15.3
were verified locally. No software was installed.

```
python tests_energy.py
python legacy_tests.py
python run_energy_examples.py
python run_fitted_field_check.py
```

The extension uses dimensionless chemical potentials
`mu_plus/kBT=ln(c_plus)+u+w` and `mu_minus/kBT=ln(c_minus)-u+w`.
Every SG face, including a membrane/bath energy jump, uses `Delta(z*u+w)`.
External reservoirs and finite outer baths have w=0. A prescribed discontinuity
has exact discrete Boltzmann equilibrium. For non-equilibrium sharp interfaces,
the SG interpolation represents the jump over a center-to-center edge; the
continuum jump-resistance error decreases with refinement and is separately
checked against an independent 1D analytic resistance. Do not infer second-order
convergence across a discontinuity from SG's equilibrium consistency.

`CommonIonEnergy('rectangle', amplitude_kbt=w0)` gives w0 in the membrane water.
`CommonIonEnergy('two_sided_exponential', amplitude_kbt=ws, decay_nm=ell)` gives
`ws*(exp(-(x-xL)/ell)+exp(-(xR-x)/ell))` in that water. These are imposed effective
structure/interface entry energies U/kBT, not energies computed from microscopic
PE dielectric physics. Negative amplitudes are allowed as exploratory attraction.
The two amplitudes are not interchangeable: rectangle w0 is a plateau; exponential
ws is the one-face coefficient, so finite-width overlap contributes at each face.

The axial coordinate x is the original array first coordinate, called z in v3.
Defaults: xL=external_water_each_nm, xR=xL+network_span_nm. The field excludes
solid fiber cells. `fiber_side_nm=26.6` means the side of a representative square
solid cross-section. `nominal_open_gap_nm=53.2` is a periodic lateral gap.
`network_span_nm=(2*fiber_rows-1)*26.6`; five rows span 239.4nm.
None of these is a measured wet pore radius or the reduced straight-channel
model's radius a. Those models have distinct geometries.

Example API:
```
from solver import solve, Geometry, CommonIonEnergy
result, state, extra, mesh, history = solve(
    30., 1., sigma=0., geometry=Geometry(fiber_rows=5),
    ion_energy_model=CommonIonEnergy('two_sided_exponential',3.,26.6))
```

Fixed given sigma can be studied with the common field. Combining nonzero w
with the old unvalidated Na/Cl adsorption-chemistry candidate is explicitly
rejected, rather than silently treating their parameters as the same mechanism.
No flow, wetting dynamics, dry gaps, real electrodes, 3D pore reconstruction,
microscopic energy derivation, or fit to real PE samples is claimed.

Representative examples are 30/1mM at 298.15K, five idealized fiber rows, fixed
sigma=0 or -0.005C/m2, w0=3 or ws=3 with ell=26.6nm, grids8/16/24. The chosen
energy/charge parameters are demonstrations, not fitted or measured PE values.

`run_fitted_field_check.py` transfers root-supplied ws=7.03979054 and
ell=893.953885nm into this different five-row geometry for a spatial check.
It also checks ws=0/1/2 on the coarse grid. All six attempted cases converged
locally. This runner tightens the equation tolerance to 1e-13 because a very
small salt flux requires tighter absolute residuals to retain useful relative
current/ion-conservation errors. The default solver tolerance is inherited from
v3; always inspect relative conservation for strong barriers instead of treating
Newton convergence alone as physical or numerical acceptance.
