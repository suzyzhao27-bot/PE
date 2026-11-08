"""Frozen base/tight tolerance comparison at completed N0/NS B and L fits."""

from __future__ import annotations

import argparse
import csv
import json

import numpy as np

from v7_d1_l_fit_start import H_L, curves_l
from v7_d1_objective import H_VALUES, OUT, curves, graphs
from v7_network_transport import SPEC_PATH, closed_baths, system


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("model", choices=("N0", "NS"))
    args = parser.parse_args()
    full = json.loads((OUT / f"v7_{args.model}_fullB_summary.json").read_text(
        encoding="utf-8"))
    long = json.loads((OUT / f"v7_{args.model}_L_summary.json").read_text(
        encoding="utf-8"))
    if full["status"] != "multi_start_stable" or long["status"] != "multi_start_stable":
        raise RuntimeError("Formal B and L fits must be stable before T6")
    spec = json.loads(SPEC_PATH.read_text(encoding="utf-8-sig"))
    nets = graphs()
    rows = []
    for series, data, thicknesses, alpha, q in (
            ("B", curves(), H_VALUES, full["alpha_m_inv"],
             full["Q_uS_cm_per_unit_uR"]),
            ("L", curves_l(), H_L, long["alpha_L_m_inv"],
             long["Q_L_uS_cm_per_unit_uR"])):
        for h in thicknesses:
            times = np.array([t for t, _, _ in data[h]]) * 60.0
            for graph in nets[h]:
                capacities, edges, _ = system(graph, float(full["beta"]),
                                              alpha, spec)
                base = closed_baths(capacities, edges, times,
                                    rtol=1e-7, atol=1e-11)
                tight = closed_baths(capacities, edges, times,
                                     rtol=1e-9, atol=1e-13)
                y_base = q * base["u_R"]
                y_tight = q * tight["u_R"]
                scaled = np.abs(y_tight - y_base) / np.maximum(1.0, np.abs(y_tight))
                rows.append({"model": args.model, "series": series, "h_nm": h,
                             "seed": graph["seed"], "n_times": len(times),
                             "max_scaled_output_error": float(np.max(scaled)),
                             "max_absolute_output_error_uS_cm": float(np.max(
                                 np.abs(y_tight - y_base))),
                             "base_mass_relative_error": base["mass_relative_error"],
                             "tight_mass_relative_error": tight["mass_relative_error"],
                             "base_minimum_u": base["minimum_u"],
                             "tight_minimum_u": tight["minimum_u"]})
    assert len(rows) == (6 + 3) * 8
    with (OUT / f"v7_{args.model}_T6_runs.csv").open(
            "w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    summary = {"model": args.model, "n_runs": len(rows),
               "max_scaled_output_error": max(r["max_scaled_output_error"]
                                                  for r in rows),
               "max_absolute_output_error_uS_cm": max(
                   r["max_absolute_output_error_uS_cm"] for r in rows),
               "max_mass_relative_error": max(max(r["base_mass_relative_error"],
                                                   r["tight_mass_relative_error"])
                                              for r in rows),
               "minimum_u": min(min(r["base_minimum_u"], r["tight_minimum_u"])
                                for r in rows),
               "threshold_scaled_output": 0.01,
               "status": "PASS" if all(
                   r["max_scaled_output_error"] <= 0.01 and
                   r["base_mass_relative_error"] <= 1e-8 and
                   r["tight_mass_relative_error"] <= 1e-8 and
                   r["base_minimum_u"] >= -1e-12 and
                   r["tight_minimum_u"] >= -1e-12 for r in rows) else "FAIL",
               "scope": "B six and L three observed h values, eight fixed S seeds, formal fitted alpha/Q and beta"}
    (OUT / f"v7_{args.model}_T6_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"{args.model} T6: {summary['status']}; {len(rows)} runs; "
          f"max scaled={summary['max_scaled_output_error']:.3g}; "
          f"mass={summary['max_mass_relative_error']:.3g}")
    if summary["status"] != "PASS":
        raise RuntimeError("Formal T6 numerical gate failed; see saved rows")


if __name__ == "__main__":
    main()
