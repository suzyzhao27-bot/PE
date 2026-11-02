"""Reconcile existing 16/31-node and post-freeze 35-node integral comparisons."""

from __future__ import annotations

import csv
from pathlib import Path

from mechanism_transport import PoreTable


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/salt_transport"
OUT = DATA / "22_T12_10mM三网格积分分解_20261029.csv"
FIELDS = ("lower_mM", "upper_mM", "integral_16_mol_m_s", "integral_31_mol_m_s",
          "integral_35_mol_m_s", "difference_31_minus_16_mol_m_s",
          "difference_35_minus_31_mol_m_s", "new_PB_points_in_interval")


def read(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    main_pb = read("22_pb_main_20261028.csv")
    extras = read("22_pb_target_refinement_20261028.csv")
    spots = read("22_pb_postfreeze_spotcheck_20261029.csv")
    if len(main_pb) != 372 or len(spots) != 4:
        raise ArithmeticError("Unexpected fixed table or new-point size")
    family = sorted((float(r["c_mol_m3"]), float(r["D_s_local_m2_s"])) for r in main_pb
                    if float(r["a_nm"]) == 25 and float(r["sigma_C_m2"]) == -.008)
    endpoints = sorted((float(r["c_mol_m3"]), float(r["D_s_local_m2_s"])) for r in extras
                       if float(r["a_nm"]) == 25 and float(r["sigma_C_m2"]) == -.008)
    new_points = sorted((float(r["c_mol_m3"]), float(r["D_s_local_m2_s"])) for r in spots)
    if len(family) != 31 or [c for c, _ in endpoints] != [1, 10]:
        raise ArithmeticError("Unexpected target PB family")
    fine, coarse, augmented = PoreTable(), PoreTable(), PoreTable()
    if fine.nodes != sorted([(0.0, 0.0)] + family + endpoints):
        raise ArithmeticError("Original transport table changed")
    coarse.nodes = sorted([(0.0, 0.0)] + family[::2] + endpoints)
    augmented.nodes = sorted(fine.nodes + new_points)
    bounds = [c for c, _ in coarse.nodes if c < 10] + [10.0]
    output = []
    for low, high in zip(bounds, bounds[1:]):
        i16, i31, i35 = (table.integral(low, high) for table in (coarse, fine, augmented))
        output.append({
            "lower_mM": low, "upper_mM": high,
            "integral_16_mol_m_s": i16, "integral_31_mol_m_s": i31,
            "integral_35_mol_m_s": i35,
            "difference_31_minus_16_mol_m_s": i31 - i16,
            "difference_35_minus_31_mol_m_s": i35 - i31,
            "new_PB_points_in_interval": sum(low < c < high for c, _ in new_points),
        })
    totals = [sum(float(r[key]) for r in output) for key in
              ("integral_16_mol_m_s", "integral_31_mol_m_s", "integral_35_mol_m_s")]
    direct = [table.integral(0, 10) for table in (coarse, fine, augmented)]
    if any(abs(a - b) > 1e-12 * b for a, b in zip(totals, direct)):
        raise ArithmeticError("Interval decomposition mismatch")
    if sum(r["new_PB_points_in_interval"] for r in output) != 4:
        raise ArithmeticError("New PB points not within audit intervals")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    print(f"intervals={len(output)} integral16={totals[0]:.10g} "
          f"integral31={totals[1]:.10g} integral35={totals[2]:.10g} "
          f"relative_16_31={(totals[1]-totals[0])/totals[1]:.8g} "
          f"relative_31_35={(totals[2]-totals[1])/totals[2]:.8g} "
          f"relative_16_35={(totals[2]-totals[0])/totals[2]:.8g}")


if __name__ == "__main__":
    main()
