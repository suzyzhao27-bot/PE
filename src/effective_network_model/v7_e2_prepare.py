"""Prepare the frozen E2 paired sparse/dense synthetic observations and starts."""

from __future__ import annotations

import csv
import json

import numpy as np

from v7_d1_objective import H_VALUES, ROOT, curves, graphs
from v7_network_transport import SPEC_PATH, closed_baths, system


OUT = ROOT / "WP5"


def write_csv(path, rows):
    assert rows
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    frozen = json.loads((ROOT / "原始资料/赵珈绮/"
                         "14_论文补强学生任务书_v7.0.0/configs/frozen_spec.json"
                         ).read_text(encoding="utf-8-sig"))
    e2 = frozen["matrices"]["E2"]
    assert (e2["truth_sets"], e2["sparse_tasks"], e2["dense_tasks"],
            e2["dense_interval_min"], e2["noise_seed"],
            e2["startset_seeds"], e2["starts_per_set"]) == (
                2, 8, 8, 1, 20261023, [20261024, 20261025], 4)
    with (OUT / "v7_E1_equivalent_selection.csv").open(
            encoding="utf-8-sig", newline="") as stream:
        selection = {r["label"]: r for r in csv.DictReader(stream)}
    assert "fullB_best" in selection and "equivalent_1" in selection
    truth_labels = ("fullB_best", "equivalent_1")
    starts = []
    for set_index, seed in enumerate(e2["startset_seeds"]):
        rng = np.random.Generator(np.random.PCG64(seed))
        jitter = np.column_stack([
            (rng.permutation(4) + rng.random(4)) / 4
            for _ in range(3)])
        for i in range(4):
            x = jitter[i] * np.array([5.0, 6.0, 6.0])
            starts.append({
                "startset_index": set_index, "seed": seed,
                "start_index": i, "log10_alpha": x[0],
                "log10_Q": x[1], "beta": x[2]})
    assert len(starts) == 8

    observed = curves()
    network = graphs()
    spec = json.loads(SPEC_PATH.read_text(encoding="utf-8-sig"))
    noise_rng = np.random.Generator(np.random.PCG64(e2["noise_seed"]))
    rows = []
    diagnostics = []
    for truth_label in truth_labels:
        truth = selection[truth_label]
        alpha = 10 ** float(truth["log10_alpha"])
        q = 10 ** float(truth["log10_Q"])
        beta = float(truth["beta"])
        for h in H_VALUES:
            dense_times = np.arange(0.0, 121.0, 1.0)
            sparse_times = {t for t, _, _ in observed[h]}
            assert sparse_times <= set(dense_times.tolist())
            right = []
            mass = []
            minima = []
            for graph in network[h]:
                capacities, edges, _ = system(graph, beta, alpha, spec)
                run = closed_baths(capacities, edges, dense_times * 60.0)
                right.append(run["u_R"])
                mass.append(run["mass_relative_error"])
                minima.append(run["minimum_u"])
            assert len(right) == 8
            max_mass = max(mass)
            minimum_u = min(minima)
            assert max_mass <= 1e-8 and minimum_u >= -1e-12
            diagnostics.append({"truth_label": truth_label, "h_nm": h,
                                "max_mass_relative_error": max_mass,
                                "minimum_u": minimum_u})
            baseline = observed[h][0][1]
            raw_true = baseline + q * np.mean(np.vstack(right), axis=0)
            sigma = np.maximum(0.01, 0.01 * np.abs(raw_true))
            raw_noisy = raw_true + noise_rng.normal(0.0, sigma)
            for noise_eta, raw in ((0.0, raw_true), (0.01, raw_noisy)):
                for sampling in ("sparse_actual", "dense_1min"):
                    chosen = (set(dense_times.tolist()) if sampling == "dense_1min"
                              else sparse_times)
                    for i, t in enumerate(dense_times):
                        if t not in chosen:
                            continue
                        rows.append({
                            "truth_label": truth_label,
                            "noise_eta": noise_eta, "sampling": sampling,
                            "h_nm": h, "time_min": t,
                            "raw_kappa_uS_cm": raw[i],
                            "baseline_kappa_uS_cm": raw[0],
                            "increment_uS_cm": raw[i] - raw[0],
                            "noiseless_raw_kappa_uS_cm": raw_true[i],
                            "noise_sigma_uS_cm": 0.0 if noise_eta == 0
                                                 else sigma[i],
                            "noise_seed": "" if noise_eta == 0
                                          else e2["noise_seed"],
                            "max_mass_relative_error": max_mass,
                            "minimum_u": minimum_u,
                        })
    assert len(rows) == 2 * 2 * (6 * 7 + 6 * 121)
    assert len(diagnostics) == 12
    write_csv(OUT / "v7_E2_start_sets.csv", starts)
    write_csv(OUT / "v7_E2_synthetic_observations.csv", rows)
    write_csv(OUT / "v7_E2_truth_forward_diagnostics.csv", diagnostics)
    print(f"E2 prepared: truth_h={len(diagnostics)}; "
          f"observation_rows={len(rows)}; starts={len(starts)}")


if __name__ == "__main__":
    main()
