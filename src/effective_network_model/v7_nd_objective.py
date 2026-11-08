"""Frozen five-parameter ND B objective, activated only by saved D2 decision."""

from __future__ import annotations

import csv
import json
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from v7_d1_objective import H_VALUES, OUT, curves, graphs, residuals
from v7_dynamic_network import dynamic_closed_baths, dynamic_system
from v7_network_transport import SPEC_PATH


BOUNDS = ((0.0, 5.0), (0.0, 6.0), (0.0, 6.0), (0.0, 6.0), (0.0, 2.0))
_SPEC = None
_CURVES = None
_GRAPHS = None


def starts_nd():
    mid = [(lo + hi) / 2 for lo, hi in BOUNDS]
    rng = np.random.Generator(np.random.PCG64(20261022))
    jitter = np.column_stack([(rng.permutation(11) + rng.random(11)) / 11
                              for _ in BOUNDS])
    mapped = np.array([lo + (hi - lo) * jitter[:, k]
                       for k, (lo, hi) in enumerate(BOUNDS)]).T
    result = np.vstack([mid, mapped])
    assert result.shape == (12, 5)
    return result


def init_worker():
    global _SPEC, _CURVES, _GRAPHS
    _SPEC = json.loads(SPEC_PATH.read_text(encoding="utf-8-sig"))
    _CURVES = curves()
    _GRAPHS = {(h, g["seed"]): g for h, group in graphs().items()
               for g in group}


def forward_one(task):
    h, seed, alpha, beta, tau0, b = task
    times = np.array([t for t, _, _ in _CURVES[h]]) * 60.0
    capacities, edges, _, dynamic_edges = dynamic_system(
        _GRAPHS[h, seed], beta, alpha, tau0, b, _SPEC)
    run = dynamic_closed_baths(capacities, edges, dynamic_edges, times)
    return h, seed, run["u_R"], run["mass_relative_error"], run["minimum_u"]


def predict_nd(pool: ProcessPoolExecutor, alpha, q, beta, tau0, b):
    grouped = {h: [] for h in H_VALUES}
    tasks = [(h, seed, alpha, beta, tau0, b)
             for h in H_VALUES for seed in range(7001, 7009)]
    for h, seed, right, mass_error, minimum in pool.map(forward_one, tasks):
        grouped[h].append((seed, right, mass_error, minimum))
    predictions, diagnostics = {}, {}
    for h in H_VALUES:
        group = sorted(grouped[h], key=lambda x: x[0])
        assert [x[0] for x in group] == list(range(7001, 7009))
        right = np.vstack([x[1] for x in group])
        predictions[h] = q * np.mean(right, axis=0)
        diagnostics[h] = {"max_mass_relative_error": max(x[2] for x in group),
                          "minimum_u": min(x[3] for x in group),
                          "u_R_seed_min": np.min(right, axis=0).tolist(),
                          "u_R_seed_max": np.max(right, axis=0).tolist()}
    return predictions, diagnostics


def main():
    decision = json.loads((OUT / "v7_D2_decision.json").read_text(encoding="utf-8"))
    if not decision["ND_enabled_by_frozen_D2"]:
        raise RuntimeError("Frozen D2 did not activate ND")
    rows = []
    for index, vector in enumerate(starts_nd()):
        rows.append({"model": "ND", "start_index": index,
                     "log10_alpha": vector[0], "log10_Q": vector[1],
                     "beta": vector[2], "log10_tau0_s": vector[3], "b": vector[4]})
    with (OUT / "v7_ND_starts.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    print("ND 12 frozen starts saved; D2 activated; no fit run by this command")


if __name__ == "__main__":
    main()
