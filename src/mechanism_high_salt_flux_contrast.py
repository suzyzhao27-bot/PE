"""Propagate fixed 25 nm PB coefficients through the initial flux root."""

from __future__ import annotations

import csv
from pathlib import Path

from mechanism_transport import PoreTable, salt_flux


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/processed/salt_transport/22_高盐初始通量荷电中性对照_20261029.csv"
FIELDS = ("h_nm", "sigma_mC_m2", "external", "b", "donor_mM", "receiver_mM",
          "charged_J_mol_m2_s", "neutral_J_mol_m2_s", "charged_to_neutral_J_ratio",
          "numeric_status", "interpretation")


def main() -> None:
    rows = []
    for h in (65, 700, 2600):
        for sigma_mc in (-2, -8, -20):
            table = PoreTable(a_nm=25, sigma=sigma_mc / 1000)
            for external in ("zero", "reference_2600_neutral"):
                neutral = salt_flux(600, 0, h, 0.01, 0, external)
                charged = salt_flux(600, 0, h, 0.01, sigma_mc / 1000, external, table)
                if not 0 < charged or not 0 < neutral:
                    raise ArithmeticError("Invalid initial flux")
                rows.append({"h_nm": h, "sigma_mC_m2": sigma_mc, "external": external,
                             "b": 0.01, "donor_mM": 600, "receiver_mM": 0,
                             "charged_J_mol_m2_s": charged,
                             "neutral_J_mol_m2_s": neutral,
                             "charged_to_neutral_J_ratio": charged / neutral,
                             "numeric_status": "GRID_CHECK_PASS_600mM",
                             "interpretation": "FIXED_INITIAL_CONDITION_NOT_MEASURED"})
    if len(rows) != 18:
        raise ArithmeticError("Expected 3 thickness x 3 charge x 2 external")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    for external in ("zero", "reference_2600_neutral"):
        subset = [r["charged_to_neutral_J_ratio"] for r in rows if r["external"] == external]
        print(f"{external}: ratio={min(subset):.8g}..{max(subset):.8g}")


if __name__ == "__main__":
    main()
