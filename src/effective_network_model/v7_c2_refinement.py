"""One allowed 12x12 C2 refinement after the 6v8 threshold trigger."""

from __future__ import annotations

import csv
import json
import math
import statistics
from pathlib import Path

import numpy as np

from v7_network_geometry import generate, morphology
from v7_network_transport import SPEC_PATH, closed_baths, fixed_bath_p, system


ROOT = Path(__file__).resolve().parents[2]


def main():
    with (ROOT / "WP2/v7_c2_domain_summary.csv").open(encoding="utf-8-sig", newline="") as stream:
        previous = list(csv.DictReader(stream))
    assert len(previous) == 12 and any(r["over_10pct"] == "True" for r in previous)
    spec = json.loads(SPEC_PATH.read_text(encoding="utf-8-sig"))
    geom = spec["geometry"]
    max_nz = max(math.ceil(h / sum(morphology(h, series, geom["anchors"])[:2]))
                 for h in spec["matrices"]["C1"]["h_nm"] for series in ("H", "S"))
    geom["nx"] = geom["ny"] = 12
    rows, failures = [], []
    for series in ("H", "S"):
        for h in (65, 350, 2600):
            for seed in geom["main_seeds"]:
                graph, _audit = generate(h, series, seed, spec, max_nz)
                for beta in (0.0, 3.0):
                    try:
                        capacities, edges, area = system(graph, beta, 1000.0, spec)
                        p = fixed_bath_p(capacities, edges, area)
                        run = closed_baths(capacities, edges, np.array([0.0, 7200.0]))
                        assert run["mass_relative_error"] <= 1e-8
                        assert run["minimum_u"] >= -1e-12
                        row = {"series": series, "h_nm": h, "seed": seed, "size": 12,
                               "beta": beta, "P_m_s": p, "u_R_120min": float(run["u_R"][-1]),
                               "mass_relative_error": run["mass_relative_error"],
                               "minimum_u": run["minimum_u"], "status": "PASS", "error": ""}
                    except Exception as exc:
                        failures.append((series,h,seed,beta,repr(exc)))
                        row = {"series": series, "h_nm": h, "seed": seed, "size": 12,
                               "beta": beta, "P_m_s": "", "u_R_120min": "",
                               "mass_relative_error": "", "minimum_u": "",
                               "status": "FAIL", "error": repr(exc)}
                    rows.append(row)
                    if len(rows) % 16 == 0:
                        print(f"C2 12x12 processed {len(rows)}/96; failures {len(failures)}", flush=True)
    assert len(rows) == 96
    with (ROOT / "WP2/v7_c2_12x12_runs.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    if failures:
        raise RuntimeError(f"C2 12x12 has {len(failures)} failures; see CSV")
    summary = []
    for old in previous:
        selected = [r for r in rows if (r["series"], r["h_nm"], r["beta"])
                    == (old["series"], int(old["h_nm"]), float(old["beta"]))]
        assert len(selected) == 8
        p12 = statistics.mean(r["P_m_s"] for r in selected)
        u12 = statistics.mean(r["u_R_120min"] for r in selected)
        p8 = float(old["P8_mean_m_s"])
        u8 = float(old["u120_8_mean"])
        summary.append({"series": old["series"], "h_nm": old["h_nm"], "beta": old["beta"],
                        "P6_mean_m_s": old["P6_mean_m_s"], "P8_mean_m_s": p8,
                        "P12_mean_m_s": p12, "u120_6_mean": old["u120_6_mean"],
                        "u120_8_mean": u8, "u120_12_mean": u12,
                        "P_8v12_relative": abs(p8-p12)/p12,
                        "u120_8v12_relative": abs(u8-u12)/u12,
                        "over_10pct_8v12": abs(p8-p12)/p12 > 0.1 or abs(u8-u12)/u12 > 0.1})
    with (ROOT / "WP2/v7_c2_12x12_summary.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary[0]))
        writer.writeheader(); writer.writerows(summary)
    print(f"C2 12x12 rows={len(rows)}; failures=0; 8v12 over 10%: "
          f"{sum(r['over_10pct_8v12'] for r in summary)}/12")


if __name__ == "__main__":
    main()
