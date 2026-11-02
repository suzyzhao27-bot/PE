"""Recheck initial T12 risks after exact 1/10 mM endpoint PB supplements."""

from __future__ import annotations

import json

from mechanism_pb_integral_check import OUT as INITIAL_OUT, integrate
from mechanism_pb_main_batch import OUT as MAIN_OUT, read_csv, write_csv
from mechanism_pb_target_refinement import OUT as EXTRA_OUT
from mechanism_radial_pb import MATRIX, ROOT


OUT = ROOT / "data/processed/salt_transport/22_pb_target_recheck_20261028.csv"
FIELDS = ["a_nm", "sigma_C_m2", "feed_mM", "initial_31_16_relative_difference",
          "updated_31_16_relative_difference", "endpoint_change_to_31_relative",
          "updated_31_integral_mol_m_s", "updated_16_integral_mol_m_s", "status"]


def main() -> None:
    with MATRIX.open(encoding="utf-8-sig") as stream:
        matrix = json.load(stream)
    table = {r["id"]: r for r in read_csv(MAIN_OUT)}
    extras = read_csv(EXTRA_OUT)
    risks = [r for r in read_csv(INITIAL_OUT) if r["status"] != "PASS"]
    assert len(table) == 372 and len(extras) == len(risks) == 12
    rows = []
    for risk in risks:
        a, sigma, feed = (float(risk[k]) for k in ("a_nm", "sigma_C_m2", "feed_mM"))
        family = next(i for i in range(12) if float(matrix[31 * i]["a_nm"]) == a
                      and float(matrix[31 * i]["sigma_C_m2"]) == sigma)
        points = matrix[31 * family:31 * (family + 1)]
        nodes = [(float(p["c_mol_m3"]), float(table[p["id"]]["D_s_local_m2_s"])) for p in points]
        match = next(r for r in extras if float(r["a_nm"]) == a
                     and float(r["sigma_C_m2"]) == sigma and float(r["c_mM"]) == feed)
        endpoint = (feed, float(match["D_s_local_m2_s"]))
        c_min, d_min = nodes[0]
        fine = integrate(sorted(nodes + [endpoint]), feed) + 0.5 * c_min * d_min
        coarse = integrate(sorted(nodes[::2] + [endpoint]), feed) + 0.5 * c_min * d_min
        initial_fine = float(risk["integral_31_mol_m_s"])
        relative = abs(fine - coarse) / fine
        rows.append({
            "a_nm": a, "sigma_C_m2": sigma, "feed_mM": feed,
            "initial_31_16_relative_difference": risk["relative_difference"],
            "updated_31_16_relative_difference": relative,
            "endpoint_change_to_31_relative": abs(fine - initial_fine) / fine,
            "updated_31_integral_mol_m_s": fine,
            "updated_16_integral_mol_m_s": coarse,
            "status": "PASS_AFTER_ENDPOINT" if relative <= 0.02 else "NUMERIC_UNRESOLVED_OVER_2_PERCENT",
        })
    write_csv(OUT, rows, FIELDS)
    unresolved = [r for r in rows if r["status"].startswith("NUMERIC_UNRESOLVED")]
    print(f"rechecked={len(rows)} unresolved={len(unresolved)} "
          f"max_updated_difference={max(r['updated_31_16_relative_difference'] for r in rows):.6g} "
          f"max_endpoint_change={max(r['endpoint_change_to_31_relative'] for r in rows):.6g}")


if __name__ == "__main__":
    main()
