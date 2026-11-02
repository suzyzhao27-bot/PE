"""Compare fixed 31/16 concentration-grid integrals without new PB solves."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from mechanism_pb_main_batch import OUT as MAIN_OUT, read_csv, write_csv
from mechanism_radial_pb import MATRIX, ROOT


OUT = ROOT / "data/processed/salt_transport/22_pb_integral_check_20261028.csv"
FIELDS = ["a_nm", "sigma_C_m2", "feed_mM", "integral_31_mol_m_s",
          "integral_16_mol_m_s", "relative_difference", "zero_to_min_lower_mol_m_s",
          "zero_to_min_upper_mol_m_s", "zero_interval_upper_fraction_of_integral",
          "status"]


def integrate(nodes: list[tuple[float, float]], upper: float) -> float:
    if upper < nodes[0][0] or upper > nodes[-1][0]:
        raise ValueError("Integration endpoint outside fixed PB table")
    total = 0.0
    for (c0, d0), (c1, d1) in zip(nodes, nodes[1:]):
        if upper <= c0:
            break
        end = min(upper, c1)
        d_end = d0 + (d1 - d0) * (end - c0) / (c1 - c0)
        total += (d0 + d_end) * (end - c0) / 2
        if end == upper:
            break
    return total


def main() -> None:
    with MATRIX.open(encoding="utf-8-sig") as stream:
        matrix = json.load(stream)
    table = {r["id"]: r for r in read_csv(MAIN_OUT)}
    assert len(matrix) == len(table) == 372
    rows = []
    for family_index in range(12):
        family = matrix[31 * family_index:31 * (family_index + 1)]
        nodes = [(float(p["c_mol_m3"]), float(table[p["id"]]["D_s_local_m2_s"])) for p in family]
        coarse = nodes[::2]
        assert len(nodes) == 31 and len(coarse) == 16
        c_min, d_min = nodes[0]
        for feed in (1.0, 10.0, 600.0):  # mM = mol/m^3
            fine_integral = integrate(nodes, feed) + 0.5 * c_min * d_min
            coarse_integral = integrate(coarse, feed) + 0.5 * c_min * d_min
            relative = abs(coarse_integral - fine_integral) / fine_integral
            upper_fraction = c_min * d_min / fine_integral
            rows.append({
                "a_nm": family[0]["a_nm"], "sigma_C_m2": family[0]["sigma_C_m2"],
                "feed_mM": feed, "integral_31_mol_m_s": fine_integral,
                "integral_16_mol_m_s": coarse_integral, "relative_difference": relative,
                "zero_to_min_lower_mol_m_s": 0.0,
                "zero_to_min_upper_mol_m_s": c_min * d_min,
                "zero_interval_upper_fraction_of_integral": upper_fraction,
                "status": "PASS" if relative <= 0.02 else "NUMERIC_RISK_OVER_2_PERCENT",
            })
    write_csv(OUT, rows, FIELDS)
    failed = [r for r in rows if r["status"] != "PASS"]
    print(f"targets={len(rows)} over_2pct={len(failed)} max_relative="
          f"{max(r['relative_difference'] for r in rows):.6g} "
          f"max_zero_upper_fraction={max(r['zero_interval_upper_fraction_of_integral'] for r in rows):.6g}")
    for row in failed:
        print(f"risk a={row['a_nm']} sigma={row['sigma_C_m2']} feed={row['feed_mM']} "
              f"difference={row['relative_difference']:.6g}")


if __name__ == "__main__":
    main()
