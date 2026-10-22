"""v6 WP3A: algebraic delay/capacity diagnostics, never a fitted PE result."""

from __future__ import annotations

import csv
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "WP2/src"))
from transport import D_B_REF  # noqa: E402


OUTPUT = ROOT / "WP3/量级反算_20与60min.csv"
THICKNESSES_NM = (65, 120, 200, 350, 700, 2600)
TARGET_MIN = (20, 60)
DM_RATIOS = (1.0, 1e-8)


def main() -> None:
    fields = ["h_nm", "target_delay_min", "h2_over_Db_s", "Dreq_m2_s",
              "Dreq_over_Db", "Dm_over_Db_scenario", "chi_required_fast_exchange",
              "nonnegative_chi_possible", "epsilon_max_if_dead_end_water",
              "A_over_VR_min_formula", "interpretation"]
    rows = []
    for h_nm in THICKNESSES_NM:
        h = h_nm * 1e-9
        for time_min in TARGET_MIN:
            target_s = time_min * 60
            d_req = h * h / (6 * target_s)
            for dm_ratio in DM_RATIOS:
                dm = D_B_REF * dm_ratio
                chi_req = dm / d_req - 1
                feasible = chi_req >= 0
                rows.append({
                    "h_nm": h_nm,
                    "target_delay_min": time_min,
                    "h2_over_Db_s": f"{h*h/D_B_REF:.15g}",
                    "Dreq_m2_s": f"{d_req:.15g}",
                    "Dreq_over_Db": f"{d_req/D_B_REF:.15g}",
                    "Dm_over_Db_scenario": f"{dm_ratio:.15g}",
                    "chi_required_fast_exchange": f"{chi_req:.15g}",
                    "nonnegative_chi_possible": str(feasible).lower(),
                    "epsilon_max_if_dead_end_water": f"{1/(1+chi_req):.15g}" if feasible else "not_applicable",
                    "A_over_VR_min_formula": "gamma*(1+chi_required)" if feasible else "not_applicable",
                    "interpretation": "conditional algebraic diagnostic; not measured delay or diffusivity",
                })
    with OUTPUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    assert len(rows) == 24
    print(f"{OUTPUT}: {len(rows)} conditional rows; Db={D_B_REF:.15g} m2/s")


if __name__ == "__main__":
    main()
