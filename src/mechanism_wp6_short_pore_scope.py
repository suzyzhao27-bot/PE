"""Map the fixed WP6 parameter candidates to the stated short-pore boundary."""

from __future__ import annotations

import csv
from itertools import product
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPLIT = ROOT / "data/processed/salt_transport/22_WP6校准与留出实验行_20261029.csv"
OUT = ROOT / "data/processed/salt_transport/22_WP6候选短孔适用范围_20261029.csv"
H_NM = {"thin": 65, "middle": 700, "thick": 2600}
SIGMAS = (0, -2, -8, -20)  # mC/m2, fixed in the 22-series WP6 plan
FIELDS = ("candidate_id", "candidate_scope", "a_nm", "tau", "sigma_B_mC_m2",
          "sigma_E_mC_m2", "thickness_level", "nominal_h_nm", "tau_h_over_a",
          "geometry_status", "measured_h_status")


def candidate_id(a: int, tau: int, sb: int, se: int) -> str:
    return f"a{a}_t{tau}_B{sb}_E{se}"


def main() -> None:
    with SPLIT.open(encoding="utf-8-sig", newline="") as stream:
        split = list(csv.DictReader(stream))
    if len(split) != 10 or {r["thickness_level"] for r in split} != set(H_NM):
        raise ArithmeticError("WP6 fixed condition split changed")
    if sum(r["analysis_role"].startswith("CHARGE_CALIBRATION") for r in split) != 2:
        raise ArithmeticError("WP6 calibration split changed")
    rows = []
    for a, tau, sb, se in product((2, 10, 25, 50), (1, 2, 4), SIGMAS, SIGMAS):
        scope = "UNCONFIRMED_2NM_EXPLORATORY" if a == 2 else "MAIN_144"
        for level, h in H_NM.items():
            ratio = tau * h / a
            rows.append({"candidate_id": candidate_id(a, tau, sb, se),
                         "candidate_scope": scope, "a_nm": a, "tau": tau,
                         "sigma_B_mC_m2": sb, "sigma_E_mC_m2": se,
                         "thickness_level": level, "nominal_h_nm": h,
                         "tau_h_over_a": ratio,
                         "geometry_status": "SHORT_PORE_CONDITIONAL" if ratio < 5 else "NOT_SHORT_AT_NOMINAL_H",
                         "measured_h_status": "NOT_MEASURED"})
    main = [r for r in rows if r["candidate_scope"] == "MAIN_144"]
    if len(rows) != 576 or len(main) != 432:
        raise ArithmeticError("Fixed candidate size changed")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    short = {r["candidate_id"] for r in main if r["geometry_status"] == "SHORT_PORE_CONDITIONAL"}
    print(f"candidates=192 rows={len(rows)} main_short_candidates={len(short)} "
          f"main_short_rows={sum(r['geometry_status'] == 'SHORT_PORE_CONDITIONAL' for r in main)}")


if __name__ == "__main__":
    main()
