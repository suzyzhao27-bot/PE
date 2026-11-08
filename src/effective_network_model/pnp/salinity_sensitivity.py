"""Conditional salinity scan for the frozen five-row v5 PNP example.

The three bath pairs keep c_left/c_right=30. Results are model outputs,
not concentrations reconstructed from the old PE conductivity records.
"""

import argparse
import csv
import math
from pathlib import Path

from solver import CommonIonEnergy, Geometry, solve


BATH_PAIRS_MM = ((15.0, 0.5), (30.0, 1.0), (60.0, 2.0))
SIGMA_C_M2 = (0.0, -0.005)
DEFAULT_CELLS_PER_FIBER = 24
FIBER_ROWS = 5
WS_KBT = 3.0
ELL_NM = 26.6


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cells-per-fiber", type=int, default=DEFAULT_CELLS_PER_FIBER)
    parser.add_argument("--sigma", type=float, action="append", dest="sigmas",
                        help="Repeat to select fixed-charge values; default is 0 and -0.005 C/m2")
    args = parser.parse_args()

    geometry = Geometry(cells_per_fiber=args.cells_per_fiber, fiber_rows=FIBER_ROWS)
    energy = CommonIonEnergy("two_sided_exponential", WS_KBT, ELL_NM)
    rows = []
    for left, right in BATH_PAIRS_MM:
        for sigma in (SIGMA_C_M2 if args.sigmas is None else args.sigmas):
            base = solve(left, right, sigma, geometry=geometry)[0]
            field = solve(left, right, sigma, geometry=geometry, ion_energy_model=energy)[0]
            j0 = base["common_salt_flux_mol_m2_s"]
            j = field["common_salt_flux_mol_m2_s"]
            ratio = j / j0
            checks = (j0, j, ratio, base["max_scaled_equation_residual"],
                      field["max_scaled_equation_residual"],
                      base["plus_conservation_relative"], base["minus_conservation_relative"],
                      field["plus_conservation_relative"], field["minus_conservation_relative"])
            if not all(math.isfinite(value) for value in checks) or j0 <= 0 or j <= 0:
                raise ValueError(f"Invalid flux or diagnostics at {left}/{right} mM, sigma={sigma}")
            rows.append(dict(
                c_left_mM=left, c_right_mM=right, concentration_ratio=left / right,
                sigma_c_m2=sigma, cells_per_fiber=args.cells_per_fiber, fiber_rows=FIBER_ROWS,
                energy_kind="two_sided_exponential", ws_kbt=WS_KBT, ell_nm=ELL_NM,
                zero_energy_J_mol_m2_s=j0, energy_J_mol_m2_s=j,
                J_over_same_sigma_zero_energy=ratio,
                zero_newton_iterations=base["newton_iterations"],
                energy_newton_iterations=field["newton_iterations"],
                zero_max_scaled_equation_residual=base["max_scaled_equation_residual"],
                energy_max_scaled_equation_residual=field["max_scaled_equation_residual"],
                zero_max_ion_conservation_relative=max(base["plus_conservation_relative"], base["minus_conservation_relative"]),
                energy_max_ion_conservation_relative=max(field["plus_conservation_relative"], field["minus_conservation_relative"]),
            ))

    if not args.output.parent.is_dir():
        raise FileNotFoundError(f"Output directory must already exist: {args.output.parent}")
    with args.output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print(f"{row['c_left_mM']:g}/{row['c_right_mM']:g} mM "
              f"sigma={row['sigma_c_m2']:g} J/J0={row['J_over_same_sigma_zero_energy']:.9g}")


if __name__ == "__main__":
    main()
