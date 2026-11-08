"""Frozen C2 domain-size comparison on three thicknesses and eight seeds."""

from __future__ import annotations

import copy
import csv
import gzip
import json
import math
import statistics
from pathlib import Path

import numpy as np

from v7_network_geometry import generate, morphology
from v7_network_transport import SPEC_PATH, closed_baths, fixed_bath_p, system


ROOT = Path(__file__).resolve().parents[2]
H_VALUES = (65, 350, 2600)
SIZES = (4, 6, 8)


def main():
    spec = json.loads(SPEC_PATH.read_text(encoding="utf-8-sig"))
    geom = spec["geometry"]
    max_nz = max(math.ceil(h / sum(morphology(h, series, geom["anchors"])[:2]))
                 for series in ("H", "S") for h in spec["matrices"]["C1"]["h_nm"])
    base = {}
    with gzip.open(ROOT / "WP2/v7_network_tables.jsonl.gz", "rt", encoding="utf-8") as stream:
        for line in stream:
            graph = json.loads(line)
            base[(graph["series"], graph["h_nm"], graph["seed"])] = graph
    rows, failures = [], []
    for size in SIZES:
        variant = copy.deepcopy(spec)
        variant["geometry"]["nx"] = size
        variant["geometry"]["ny"] = size
        for series in ("H", "S"):
            for h in H_VALUES:
                for seed in geom["main_seeds"]:
                    graph = base[(series, h, seed)] if size == 6 else generate(h, series, seed, variant, max_nz)[0]
                    for beta in (0.0, 3.0):
                        try:
                            capacities, edges, area = system(graph, beta, 1000.0, variant)
                            p = fixed_bath_p(capacities, edges, area)
                            run = closed_baths(capacities, edges, np.array([0.0, 7200.0]))
                            assert run["mass_relative_error"] <= 1e-8
                            assert run["minimum_u"] >= -1e-12
                            row = {"series": series, "h_nm": h, "seed": seed, "size": size,
                                   "beta": beta, "P_m_s": p,
                                   "u_R_120min": float(run["u_R"][-1]),
                                   "mass_relative_error": run["mass_relative_error"],
                                   "minimum_u": run["minimum_u"], "status": "PASS", "error": ""}
                        except Exception as exc:
                            failures.append((series,h,seed,size,beta,repr(exc)))
                            row = {"series": series, "h_nm": h, "seed": seed, "size": size,
                                   "beta": beta, "P_m_s": "", "u_R_120min": "",
                                   "mass_relative_error": "", "minimum_u": "",
                                   "status": "FAIL", "error": repr(exc)}
                        rows.append(row)
                        if len(rows) % 48 == 0:
                            print(f"C2 domain processed {len(rows)}/288; failures {len(failures)}", flush=True)
    assert len(rows) == 288
    with (ROOT / "WP2/v7_c2_domain_runs.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    if failures:
        raise RuntimeError(f"C2 domain has {len(failures)} failures; see CSV")
    summary = []
    for series in ("H", "S"):
        for h in H_VALUES:
            for beta in (0.0, 3.0):
                selected = [r for r in rows if r["series"] == series and r["h_nm"] == h and r["beta"] == beta]
                means = {size: {key: statistics.mean(r[key] for r in selected if r["size"] == size)
                                for key in ("P_m_s", "u_R_120min")} for size in SIZES}
                assert all(sum(r["size"] == size for r in selected) == 8 for size in SIZES)
                summary.append({"series": series, "h_nm": h, "beta": beta,
                                "P4_mean_m_s": means[4]["P_m_s"],
                                "P6_mean_m_s": means[6]["P_m_s"],
                                "P8_mean_m_s": means[8]["P_m_s"],
                                "u120_4_mean": means[4]["u_R_120min"],
                                "u120_6_mean": means[6]["u_R_120min"],
                                "u120_8_mean": means[8]["u_R_120min"],
                                "P_6v8_relative": abs(means[6]["P_m_s"]-means[8]["P_m_s"])/means[8]["P_m_s"],
                                "u120_6v8_relative": abs(means[6]["u_R_120min"]-means[8]["u_R_120min"])/means[8]["u_R_120min"],
                                "over_10pct": False})
    for row in summary:
        row["over_10pct"] = row["P_6v8_relative"] > 0.1 or row["u120_6v8_relative"] > 0.1
    with (ROOT / "WP2/v7_c2_domain_summary.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary[0]))
        writer.writeheader(); writer.writerows(summary)
    print(f"C2 domain rows={len(rows)}; failures=0; 6v8 over 10%: {sum(r['over_10pct'] for r in summary)}/12")


if __name__ == "__main__":
    main()
