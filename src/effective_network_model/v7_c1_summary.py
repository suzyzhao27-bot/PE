"""Summarize the fixed C1 main matrix without fitting or seed selection."""

from __future__ import annotations

import csv
import gzip
import json
import math
import statistics
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def main():
    with (ROOT / "WP2/v7_c1_main_matrix.csv").open(encoding="utf-8-sig", newline="") as stream:
        raw = list(csv.DictReader(stream))
    assert len(raw) == 384 and all(row["status"] == "PASS" for row in raw)
    nz = {}
    with gzip.open(ROOT / "WP2/v7_network_tables.jsonl.gz", "rt", encoding="utf-8") as stream:
        for line in stream:
            graph = json.loads(line)
            nz[(graph["series"], graph["h_nm"], graph["seed"])] = len({node["z_m"] for node in graph["nodes"]})
    rows = []
    for series in ("H", "S"):
        for beta in (0.0, 1.0, 3.0, 6.0):
            for h in (65, 120, 200, 350, 700, 2600):
                selected = [r for r in raw if (r["series"], int(r["h_nm"]), float(r["beta"]))
                            == (series, h, beta)]
                assert len(selected) == 8
                for r in selected:
                    assert int(r["seed"]) in range(7001, 7009)
                rows.append({"series": series, "beta": beta, "h_nm": h,
                             "nz": nz[(series, h, 7001)],
                             "P_mean_m_s": statistics.mean(float(r["P_m_s"]) for r in selected),
                             "P_min_m_s": min(float(r["P_m_s"]) for r in selected),
                             "P_max_m_s": max(float(r["P_m_s"]) for r in selected),
                             "u_R_120min_mean": statistics.mean(float(r["u_R_120min"]) for r in selected),
                             "u_R_120min_min": min(float(r["u_R_120min"]) for r in selected),
                             "u_R_120min_max": max(float(r["u_R_120min"]) for r in selected),
                             "capacity_mean_m3": statistics.mean(float(r["salt_capacity_m3"]) for r in selected),
                             "port_dissipation_mean_fraction": statistics.mean(float(r["port_dissipation_fraction"]) for r in selected),
                             "P_ratio_previous": "", "u120_ratio_previous": "",
                             "delta_lnP_delta_lnh": "", "layer_count_change": ""})
    for series in ("H", "S"):
        for beta in (0.0, 1.0, 3.0, 6.0):
            selected = [r for r in rows if r["series"] == series and r["beta"] == beta]
            for previous, current in zip(selected, selected[1:]):
                current["P_ratio_previous"] = current["P_mean_m_s"] / previous["P_mean_m_s"]
                current["u120_ratio_previous"] = current["u_R_120min_mean"] / previous["u_R_120min_mean"]
                current["delta_lnP_delta_lnh"] = math.log(current["P_ratio_previous"]) / math.log(current["h_nm"] / previous["h_nm"])
                current["layer_count_change"] = current["nz"] != previous["nz"]
    assert len(rows) == 48
    with (ROOT / "WP2/v7_c1_main_summary.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print("C1 main summary: 48 rows, 8 fixed seeds per row; no fitted parameters")


if __name__ == "__main__":
    main()
