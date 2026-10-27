"""Three prescribed C1 beta=3 single-factor ablations, eight fixed seeds."""

from __future__ import annotations

import copy
import csv
import gzip
import json
from pathlib import Path

import numpy as np

from v7_c1_main import connected_volume_fraction, steady_details
from v7_network_geometry import generate
from v7_network_transport import SPEC_PATH, closed_baths, fixed_bath_p, system


ROOT = Path(__file__).resolve().parents[2]


def main():
    spec = json.loads(SPEC_PATH.read_text(encoding="utf-8-sig"))
    graphs = []
    with gzip.open(ROOT / "WP2/v7_network_tables.jsonl.gz", "rt", encoding="utf-8") as stream:
        graphs = [json.loads(line) for line in stream]
    assert len(graphs) == 96
    max_nz = max(len({node["z_m"] for node in graph["nodes"]}) for graph in graphs)
    rows, failures = [], []
    for ablation in ("no_access", "q_1", "p_z_1"):
        variant = copy.deepcopy(spec)
        if ablation == "q_1":
            variant["geometry"]["q"] = 1.0
        elif ablation == "p_z_1":
            variant["geometry"]["p_z"] = 1.0
        for base_graph in graphs:
            identity = (ablation, base_graph["series"], base_graph["h_nm"], base_graph["seed"])
            try:
                graph = base_graph
                if ablation != "no_access":
                    graph, _audit = generate(base_graph["h_nm"], base_graph["series"],
                                             base_graph["seed"], variant, max_nz)
                capacities, edges, area = system(graph, 3.0, 1000.0, variant,
                                                 access=ablation != "no_access")
                p = fixed_bath_p(capacities, edges, area)
                p_detail, port_fraction, internal_fraction, energy_error, min_steady = steady_details(
                    capacities, edges, graph, area)
                run = closed_baths(capacities, edges, np.array([0.0, 7200.0]))
                assert abs(p-p_detail) <= 1e-8*max(abs(p),1e-300)
                assert energy_error <= 1e-8 and run["mass_relative_error"] <= 1e-8
                assert run["minimum_u"] >= -1e-12 and min_steady >= -1e-12
                row = {"ablation": ablation, "series": identity[1], "h_nm": identity[2],
                       "seed": identity[3], "beta": 3.0, "alpha_m_inv": 1000.0, "Q": 1.0,
                       "P_m_s": p, "u_R_120min": float(run["u_R"][-1]),
                       "salt_capacity_m3": float(sum(capacities[2:])),
                       "connected_volume_fraction": connected_volume_fraction(graph),
                       "port_dissipation_fraction": port_fraction,
                       "internal_dissipation_fraction": internal_fraction,
                       "energy_balance_relative_error": energy_error,
                       "mass_relative_error": run["mass_relative_error"],
                       "minimum_u": run["minimum_u"], "status": "PASS", "error": ""}
            except Exception as exc:
                failures.append((identity, repr(exc)))
                row = {"ablation": ablation, "series": identity[1], "h_nm": identity[2],
                       "seed": identity[3], "beta": 3.0, "alpha_m_inv": 1000.0, "Q": 1.0,
                       "P_m_s": "", "u_R_120min": "", "salt_capacity_m3": "",
                       "connected_volume_fraction": "", "port_dissipation_fraction": "",
                       "internal_dissipation_fraction": "", "energy_balance_relative_error": "",
                       "mass_relative_error": "", "minimum_u": "", "status": "FAIL",
                       "error": repr(exc)}
            rows.append(row)
            if len(rows) % 48 == 0:
                print(f"C1 ablations processed {len(rows)}/288; failures {len(failures)}", flush=True)
    assert len(rows) == 288
    with (ROOT / "WP2/v7_c1_ablations.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"C1 ablations rows={len(rows)}; failures={len(failures)}")
    if failures:
        raise RuntimeError(f"C1 ablations has {len(failures)} failures; see CSV")


if __name__ == "__main__":
    main()
