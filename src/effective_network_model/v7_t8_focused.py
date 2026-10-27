"""Focused ND T8 checks before any five-parameter fit."""

from __future__ import annotations

import json

import numpy as np

from v7_dynamic_network import (dynamic_closed_baths, dynamic_system,
                                effective_time_b0, free_energy)
from v7_network_transport import ROOT, SPEC_PATH, closed_baths, load_graph


def scaled_error(a, b):
    return float(np.max(np.abs(a - b)) / max(1.0, float(np.max(np.abs(b)))))


def main():
    spec = json.loads(SPEC_PATH.read_text(encoding="utf-8-sig"))
    graph = load_graph("H", 200, 7001)
    times = np.array([0.0, 0.001, 0.01, 0.1, 1.0, 5.0, 20.0, 60.0, 120.0]) * 60.0
    results = []
    for beta, tau0, b, f0 in ((3.0, 60.0, 0.0, 0.0),
                              (3.0, 60.0, 2.0, 0.0),
                              (3.0, 1.0, 2.0, 0.0),
                              (3.0, 1e6, 2.0, 0.0),
                              (3.0, 60.0, 2.0, 1.0),
                              (3.0, 0.0, 2.0, 0.0)):
        capacities, edges, _, dynamic_edges = dynamic_system(
            graph, beta, 1e3, tau0, b, spec)
        run = dynamic_closed_baths(capacities, edges, dynamic_edges, times, f0=f0)
        energy = free_energy(capacities, run["all_u"])
        norm = max(abs(float(energy[0])), float(capacities[0]))
        drift = float(max(0.0, np.max(np.diff(energy))) / norm)
        negative_salt_fraction = float(max(
            0.0, -np.min(capacities[:, None] * run["all_u"])) / capacities[0])
        assert run["mass_relative_error"] <= 1e-8
        assert negative_salt_fraction <= 1e-12
        assert drift <= 1e-10
        record = {"series": "H", "h_nm": 200, "seed": 7001,
                  "alpha_m_inv": 1e3, "beta": beta, "tau0_s": tau0,
                  "b": b, "f0": f0,
                  "max_mass_relative_error": run["mass_relative_error"],
                  "minimum_u": run["minimum_u"],
                  "max_negative_salt_fraction": negative_salt_fraction,
                  "free_energy_positive_drift_normalized": drift,
                  "u_R": run["u_R"].tolist()}
        if b == 0 and tau0 > 0:
            transformed = effective_time_b0(times, tau0, f0)
            analytic = closed_baths(capacities, edges, transformed)
            record["b0_time_transform_scaled_error"] = scaled_error(
                run["u_R"], analytic["u_R"])
            assert record["b0_time_transform_scaled_error"] <= 0.01
        if f0 == 1 or tau0 == 0:
            static = closed_baths(capacities, edges, times)
            record["NS_limit_scaled_error"] = scaled_error(
                run["u_R"], static["u_R"])
            assert record["NS_limit_scaled_error"] <= 0.01
        results.append(record)
    path = ROOT / "WP2/v7_t8_focused.json"
    path.write_text(json.dumps({"cases": results,
                                "scope": "focused prefit T8; formal ND best point pending"},
                               ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"T8 focused: {len(results)} cases; max mass error "
          f"{max(r['max_mass_relative_error'] for r in results):.3g}; "
          f"max energy drift {max(r['free_energy_positive_drift_normalized'] for r in results):.3g}")


if __name__ == "__main__":
    main()
