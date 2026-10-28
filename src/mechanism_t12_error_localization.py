"""Localize 31/16 table integral differences without new PB solves."""

from __future__ import annotations

import csv
from pathlib import Path

from mechanism_transport import PoreTable


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/salt_transport"
OUT = DATA / "22_T12_10mM积分差区间定位_20261029.csv"
FIELDS = ("lower_mM", "upper_mM", "fine_integral_mol_m_s", "coarse_integral_mol_m_s",
          "signed_difference_mol_m_s", "absolute_difference_mol_m_s",
          "share_of_total_signed_difference", "fine_interior_main_nodes")


def read(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    main_rows = read("22_pb_main_20261028.csv")
    extras = read("22_pb_target_refinement_20261028.csv")
    family = sorted((float(r["c_mol_m3"]), float(r["D_s_local_m2_s"])) for r in main_rows
                    if float(r["a_nm"]) == 25 and float(r["sigma_C_m2"]) == -0.008)
    endpoint = sorted((float(r["c_mol_m3"]), float(r["D_s_local_m2_s"])) for r in extras
                      if float(r["a_nm"]) == 25 and float(r["sigma_C_m2"]) == -0.008)
    if len(family) != 31 or [c for c, _ in endpoint] != [1.0, 10.0]:
        raise ArithmeticError("Unexpected PB family")
    fine, coarse = PoreTable(), PoreTable()
    coarse.nodes = sorted([(0.0, 0.0)] + family[::2] + endpoint)
    bounds = [c for c, _ in coarse.nodes if 0 <= c < 10] + [10.0]
    rows = []
    for low, high in zip(bounds, bounds[1:]):
        a, b = fine.integral(low, high), coarse.integral(low, high)
        rows.append({"lower_mM": low, "upper_mM": high,
                     "fine_integral_mol_m_s": a, "coarse_integral_mol_m_s": b,
                     "signed_difference_mol_m_s": a - b,
                     "absolute_difference_mol_m_s": abs(a - b),
                     "fine_interior_main_nodes": sum(low < c < high for c, _ in family)})
    total = sum(r["signed_difference_mol_m_s"] for r in rows)
    direct = fine.integral(0, 10) - coarse.integral(0, 10)
    if total <= 0 or abs(total - direct) > 1e-12 * abs(direct):
        raise ArithmeticError("Interval decomposition does not reconcile")
    for row in rows:
        row["share_of_total_signed_difference"] = row["signed_difference_mol_m_s"] / total
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    largest = sorted(rows, key=lambda r: abs(r["signed_difference_mol_m_s"]), reverse=True)
    print(f"intervals={len(rows)} relative_0_10={total/fine.integral(0, 10):.6g} "
          f"largest={largest[0]['lower_mM']:.6g}..{largest[0]['upper_mM']:.6g} "
          f"top3_signed_share={sum(r['signed_difference_mol_m_s'] for r in largest[:3])/total:.6g}")


if __name__ == "__main__":
    main()
