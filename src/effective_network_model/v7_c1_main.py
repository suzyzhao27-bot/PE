"""Frozen v7 C1 main matrix: no fitting, no geometry regeneration."""

from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import spsolve

from v7_network_transport import SPEC_PATH, closed_baths, fixed_bath_p, generator, system


ROOT = Path(__file__).resolve().parents[2]


def steady_details(capacities, edges, graph, area):
    operator = generator(capacities, edges)
    membrane = operator[2:, 2:]
    source = np.zeros(len(capacities) - 2)
    for i, j, g in edges:
        if i == 0 and j >= 2:
            source[j - 2] += g / capacities[j]
        elif j == 0 and i >= 2:
            source[i - 2] += g / capacities[i]
    states = np.r_[1.0, 0.0, spsolve(-membrane, source)]
    flux = sum(g * (states[i] - states[j]) for i, j, g in edges if i == 0 or j == 0)
    port_power = 0.0
    internal_power = 0.0
    for edge, (i, j, g) in zip(graph["edges"], edges, strict=True):
        power = g * (states[i] - states[j])**2
        if edge["port_flag"]:
            port_power += power
        else:
            internal_power += power
    total_power = port_power + internal_power
    return (float(flux / area), float(port_power / total_power),
            float(internal_power / total_power),
            float(abs(total_power - flux) / max(abs(flux), 1e-300)),
            float(states.min()))


def connected_volume_fraction(graph):
    n = len(graph["nodes"])
    parent = list(range(n + 2))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for edge in graph["edges"]:
        i = 0 if edge["i"] == "L" else 1 if edge["i"] == "R" else edge["i"] + 2
        j = 0 if edge["j"] == "L" else 1 if edge["j"] == "R" else edge["j"] + 2
        parent[find(i)] = find(j)
    both = find(0) == find(1)
    volumes = np.array([node["V_m3"] for node in graph["nodes"]])
    connected = sum(volumes[i] for i in range(n) if both and find(i + 2) == find(0))
    return float(connected / sum(volumes))


def main():
    spec = json.loads(SPEC_PATH.read_text(encoding="utf-8-sig"))
    rows = []
    failures = []
    path = ROOT / "WP2/v7_network_tables.jsonl.gz"
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            graph = json.loads(line)
            connected_fraction = connected_volume_fraction(graph)
            for beta in (0.0, 1.0, 3.0, 6.0):
                identity = (graph["series"], graph["h_nm"], graph["seed"], beta)
                try:
                    capacities, edges, area = system(graph, beta, 1000.0, spec)
                    p = fixed_bath_p(capacities, edges, area)
                    p_detail, port_fraction, internal_fraction, energy_error, min_steady = steady_details(
                        capacities, edges, graph, area)
                    run = closed_baths(capacities, edges, np.array([0.0, 7200.0]))
                    u_r = float(run["u_R"][-1])
                    mass_error = run["mass_relative_error"]
                    min_u = run["minimum_u"]
                    salt_capacity = float(sum(capacities[2:]))
                    assert abs(p - p_detail) <= 1e-8 * max(abs(p), 1e-300)
                    assert energy_error <= 1e-8
                    assert mass_error <= 1e-8
                    assert min_u >= -1e-12 and min_steady >= -1e-12
                    row = {"series": identity[0], "h_nm": identity[1], "seed": identity[2],
                           "beta": beta, "alpha_m_inv": 1000.0, "Q": 1.0,
                           "P_m_s": p, "u_R_120min": u_r,
                           "salt_capacity_m3": salt_capacity,
                           "connected_volume_fraction": connected_fraction,
                           "port_dissipation_fraction": port_fraction,
                           "internal_dissipation_fraction": internal_fraction,
                           "energy_balance_relative_error": energy_error,
                           "mass_relative_error": mass_error, "minimum_u": min_u,
                           "status": "PASS", "error": ""}
                except Exception as exc:
                    failures.append((identity, repr(exc)))
                    row = {"series": identity[0], "h_nm": identity[1], "seed": identity[2],
                           "beta": beta, "alpha_m_inv": 1000.0, "Q": 1.0,
                           "P_m_s": "", "u_R_120min": "", "salt_capacity_m3": "",
                           "connected_volume_fraction": connected_fraction,
                           "port_dissipation_fraction": "", "internal_dissipation_fraction": "",
                           "energy_balance_relative_error": "", "mass_relative_error": "",
                           "minimum_u": "", "status": "FAIL", "error": repr(exc)}
                rows.append(row)
                if len(rows) % 48 == 0:
                    print(f"C1 main processed {len(rows)}/384; failures {len(failures)}", flush=True)
    assert len(rows) == 384
    output = ROOT / "WP2/v7_c1_main_matrix.csv"
    with output.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"C1 main rows={len(rows)}; failures={len(failures)}; output={output}")
    if failures:
        raise RuntimeError(f"C1 main has {len(failures)} failures; see CSV")


if __name__ == "__main__":
    main()
