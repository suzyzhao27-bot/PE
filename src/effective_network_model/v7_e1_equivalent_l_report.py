"""Audit E1 equivalent-beta L early fits and late conditional predictions."""

from __future__ import annotations

import csv
import json
import math

import numpy as np

from v7_d1_l_fit_start import H_L, curves_l, residuals_early
from v7_d1_objective import OUT, ROOT, graphs
from v7_d1_report_l import predict_full
from v7_network_transport import SPEC_PATH


E1 = ROOT / "WP5"


def write_csv(path, rows):
    assert rows
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    with (E1 / "v7_E1_equivalent_selection.csv").open(
            encoding="utf-8-sig", newline="") as stream:
        selection = list(csv.DictReader(stream))
    assert [r["label"] for r in selection] == [
        "fullB_best", "equivalent_1", "equivalent_2", "equivalent_3"]
    runs = [json.loads(line) for line in
            (E1 / "v7_E1_equivalent_L_fit_runs.jsonl").read_text(
                encoding="utf-8").splitlines()]
    assert len(runs) == 12
    assert sorted((r["label"], r["start_index"]) for r in runs) == [
        (f"equivalent_{rank}", i) for rank in range(1, 4) for i in range(4)]
    baseline = json.loads((OUT / "v7_NS_L_summary.json").read_text(
        encoding="utf-8"))
    assert baseline["status"] == "multi_start_stable"
    data, nets = curves_l(), graphs()
    spec = json.loads(SPEC_PATH.read_text(encoding="utf-8-sig"))
    all_points, summaries = [], []
    for item in selection:
        label = item["label"]
        beta = float(item["beta"])
        if label == "fullB_best":
            x = baseline["best_transformed"]
            assert math.isclose(beta, baseline["fixed_beta_from_full_B"],
                                abs_tol=1e-12)
            recorded_L = baseline["L_early_recomputed"]
            status_counts = json.dumps(baseline["four_start_status_counts"],
                                       sort_keys=True)
            actual_calls = baseline["actual_objective_calls_total"]
        else:
            group = [r for r in runs if r["label"] == label]
            assert len(group) == 4
            assert all(math.isclose(beta, r["fixed_beta_from_selected_B"],
                                    abs_tol=1e-12) for r in group)
            good = [r for r in group if r["status"] == "solver_success"
                    and r["best_L"] is not None]
            assert good, f"{label}: no successful L early fit"
            best = min(good, key=lambda r: r["best_L"])
            assert any(r is not best and
                       abs(r["best_L"] - best["best_L"]) /
                       max(best["best_L"], 1e-12) <= 0.01
                       for r in good), f"{label}: best L not reproduced within 1%"
            assert not any(r["status"] in ("max_eval", "max_eval_or_solver_stop")
                           and r["best_L"] is not None
                           and r["best_L"] < best["best_L"] for r in group), (
                               f"{label}: lower unfinished L start requires review")
            x = best["best_transformed"]
            recorded_L = best["best_L"]
            status_counts = json.dumps({
                status: sum(r["status"] == status for r in group)
                for status in sorted({r["status"] for r in group})},
                sort_keys=True)
            actual_calls = sum(r["actual_objective_calls"] for r in group)
        alpha_l, q_l = 10 ** x[0], 10 ** x[1]
        predicted, diagnostic = predict_full(
            data, nets, spec, alpha_l, beta, q_l)
        early = {h: predicted[h][
            [i for i, (t, _, _) in enumerate(data[h]) if t <= 240]]
            for h in H_L}
        r = residuals_early(data, early)
        early_L = float(r @ r)
        assert math.isclose(early_L, recorded_L, rel_tol=1e-6, abs_tol=1e-8)
        max_mass = max(d["max_mass_relative_error"]
                       for d in diagnostic.values())
        min_state = min(d["minimum_u"] for d in diagnostic.values())
        assert max_mass <= 1e-8 and min_state >= -1e-12
        per_late = []
        for h in H_L:
            baseline_k = data[h][0][1]
            late_r = []
            for index, ((t, k, curve_id), y) in enumerate(
                    zip(data[h], predicted[h], strict=True)):
                observed = k - baseline_k
                error = math.asinh(y) - math.asinh(observed)
                role = ("baseline" if index == 0 else
                        "early_train" if t <= 240 else "late_prediction")
                if role == "late_prediction":
                    late_r.append(error)
                all_points.append({
                    "label": label, "h_nm": h, "curve_id": curve_id,
                    "point_index": index, "time_min": t, "role": role,
                    "observed_increment_uS_cm": observed,
                    "predicted_increment_uS_cm": y,
                    "prediction_seed_min_uS_cm":
                        diagnostic[h]["seed_min"][index],
                    "prediction_seed_max_uS_cm":
                        diagnostic[h]["seed_max"][index],
                    "asinh_residual": error,
                    "max_mass_relative_error":
                        diagnostic[h]["max_mass_relative_error"],
                    "minimum_u": diagnostic[h]["minimum_u"],
                })
            assert late_r
            per_late.append(float(np.sqrt(np.mean(np.square(late_r)))))
        summaries.append({
            "label": label, "fixed_beta_from_B": beta,
            "log10_alpha_L": x[0], "log10_Q_L": x[1],
            "L_early_recomputed": early_L,
            "E_Llate": float(np.sqrt(np.mean(np.square(per_late)))),
            "E_Llate_h200": per_late[0], "E_Llate_h350": per_late[1],
            "E_Llate_h700": per_late[2],
            "start_status_counts": status_counts,
            "actual_objective_calls_total": actual_calls,
            "max_mass_relative_error": max_mass,
            "minimum_u": min_state,
            "scope": "L t<=240 early fit; t>240 conditional prediction",
        })
    assert len(summaries) == 4 and len(all_points) == 4 * (19 + 15 + 36)
    write_csv(E1 / "v7_E1_equivalent_L_summary.csv", summaries)
    write_csv(E1 / "v7_E1_equivalent_L_points.csv", all_points)
    print("E1 equivalent L: four audited conditional transfers; "
          f"late points table rows={len(all_points)}")


if __name__ == "__main__":
    main()
