"""Frozen ND T8 checks at fitted B point and tau0 bounds."""

from __future__ import annotations

import csv
import json

import numpy as np

from v7_d1_objective import H_VALUES, OUT, curves, graphs
from v7_dynamic_network import (dynamic_closed_baths, dynamic_system,
                                effective_time_b0, free_energy)
from v7_network_transport import SPEC_PATH, closed_baths


def scaled_error(actual, reference):
    return float(np.max(np.abs(actual - reference)) /
                 max(1.0, float(np.max(np.abs(reference)))))


def check_run(capacities, run):
    energy = free_energy(capacities, run["all_u"])
    initial_mass = float(capacities[0])
    norm = max(abs(float(energy[0])), initial_mass)
    return {"max_mass_relative_error": float(run["mass_relative_error"]),
            "minimum_u": float(run["minimum_u"]),
            "max_negative_salt_fraction": float(max(
                0.0, -np.min(capacities[:, None] * run["all_u"])) / initial_mass),
            "free_energy_positive_drift_normalized": float(max(
                0.0, np.max(np.diff(energy))) / norm)}


def main():
    full = json.loads((OUT / "v7_ND_fullB_summary.json").read_text(encoding="utf-8"))
    if full["model"] != "ND":
        raise RuntimeError("ND full-B report required")
    spec = json.loads(SPEC_PATH.read_text(encoding="utf-8-sig"))
    data, nets = curves(), graphs()
    alpha, beta, q, b = (float(full[key]) for key in
                         ("alpha_m_inv", "beta", "Q_uS_cm_per_unit_uR", "b"))
    scenarios = (("fitted", float(full["tau0_s"])),
                 ("tau0_lower", 1.0), ("tau0_upper", 1e6))
    rows = []
    for scenario, tau0 in scenarios:
        for h in H_VALUES:
            times = np.array([t for t, _, _ in data[h]]) * 60.0
            for graph in nets[h]:
                capacities, edges, _, dynamic_edges = dynamic_system(
                    graph, beta, alpha, tau0, b, spec)
                run = dynamic_closed_baths(capacities, edges, dynamic_edges, times)
                rows.append({"scenario": scenario, "h_nm": h, "seed": graph["seed"],
                             "alpha_m_inv": alpha, "beta": beta, "tau0_s": tau0,
                             "b": b, "Q_uS_cm_per_unit_uR": q,
                             **check_run(capacities, run)})
    assert len(rows) == 3 * 6 * 8
    with (OUT / "v7_ND_T8_runs.csv").open("w", encoding="utf-8-sig",
                                           newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)

    graph = next(g for g in nets[200] if g["seed"] == 7001)
    times = np.array([t for t, _, _ in data[200]]) * 60.0
    limits = {}
    for name, tau0, b_case, f0 in (("f0_one", float(full["tau0_s"]), b, 1.0),
                                   ("tau0_zero", 0.0, b, 0.0),
                                   ("b_zero", float(full["tau0_s"]), 0.0, 0.0)):
        capacities, edges, _, dynamic_edges = dynamic_system(
            graph, beta, alpha, tau0, b_case, spec)
        run = dynamic_closed_baths(capacities, edges, dynamic_edges, times, f0=f0)
        reference_times = (effective_time_b0(times, tau0, f0)
                           if name == "b_zero" else times)
        reference = closed_baths(capacities, edges, reference_times)
        limits[name] = {"scaled_output_error": scaled_error(q * run["u_R"],
                                                             q * reference["u_R"]),
                        **check_run(capacities, run)}
    all_checks = [*rows, *limits.values()]
    core_passed = (all(r["max_mass_relative_error"] <= 1e-8 and
                  r["minimum_u"] >= -1e-12 and
                  r["max_negative_salt_fraction"] <= 1e-12 and
                  r["free_energy_positive_drift_normalized"] <= 1e-10
                  for r in all_checks) and
              all(v["scaled_output_error"] <= 0.01 for v in limits.values()))
    t6_path = OUT / "v7_ND_T6_summary.json"
    t6_status = (json.loads(t6_path.read_text(encoding="utf-8"))["status"]
                 if t6_path.exists() else "pending")
    status = ("FAIL" if not core_passed or t6_status == "FAIL" else
              "PASS" if t6_status == "PASS" else "core_pass_T6_pending")
    summary = {"model": "ND", "full_B_fit_status": full["status"],
               "n_fitted_and_tau_bound_network_runs": len(rows),
               "max_mass_relative_error": max(r["max_mass_relative_error"]
                                              for r in all_checks),
               "minimum_u": min(r["minimum_u"] for r in all_checks),
               "max_negative_salt_fraction": max(r["max_negative_salt_fraction"]
                                                  for r in all_checks),
               "max_free_energy_positive_drift_normalized": max(
                   r["free_energy_positive_drift_normalized"] for r in all_checks),
               "limit_cases_h200_seed7001": limits,
               "core_checks_passed": core_passed, "T6_status": t6_status,
               "status": status,
               "scope": "B six h x eight S seeds at fitted tau0 and both frozen tau0 bounds; NS and b=0 limits at h200 seed7001; ND formal T6 separate"}
    (OUT / "v7_ND_T8_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"ND T8: {summary['status']}; {len(rows)} full-network runs; "
          f"mass={summary['max_mass_relative_error']:.3g}; "
          f"energy={summary['max_free_energy_positive_drift_normalized']:.3g}")
    if status == "FAIL":
        raise RuntimeError("ND T8 formal check failed; inspect saved rows")


if __name__ == "__main__":
    main()
