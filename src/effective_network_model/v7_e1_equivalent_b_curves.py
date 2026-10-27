"""Frozen E1: select distinct audited near-equivalent profile points and B curves."""

from __future__ import annotations

import csv
import json
import math
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from v7_d1_objective import OUT, ROOT, curves, init_forward_worker, predict_parallel, residuals


E1 = ROOT / "WP5"
CAP_MASS = 1e-8


def write_csv(path, rows):
    assert rows
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    full = json.loads((OUT / "v7_NS_fullB_summary.json").read_text(encoding="utf-8"))
    assert full["status"] == "multi_start_stable"
    best_L = float(full["L_recomputed"])
    ceiling = best_L + max(0.01, 0.05 * best_L)
    anchor = np.asarray(full["best_transformed"], dtype=float)
    assert anchor.shape == (3,)

    reports = [json.loads(line) for line in
               (E1 / "v7_E1_profile_point_reports.jsonl").read_text(
                   encoding="utf-8").splitlines()]
    assert len(reports) == 27
    assert sorted((r["parameter_index"], r["position_index"]) for r in reports) == [
        (p, i) for p in range(3) for i in range(9)]
    valid = []
    for r in reports:
        if (r["status"] != "multi_start_stable"
                or r["best_L"] is None
                or r["L_recomputed"] is None
                or r["max_mass_relative_error"] is None
                or r["minimum_u"] is None):
            continue
        if (r["best_L"] > ceiling or r["max_mass_relative_error"] > CAP_MASS
                or r["minimum_u"] < -1e-12
                or not math.isclose(r["best_L"], r["L_recomputed"],
                                    rel_tol=1e-6, abs_tol=1e-8)):
            continue
        x = np.asarray(r["best_transformed"], dtype=float)
        assert x.shape == (3,) and np.all(np.isfinite(x))
        valid.append((float(np.linalg.norm(x - anchor)), r))
    # Three optimizer starts at one fixed profile position validate one solution.
    # Use the audited point's selected best, never count start repeats as separate solutions.
    valid.sort(key=lambda item: (-item[0], item[1]["parameter_index"],
                                 item[1]["position_index"]))
    selected = valid[:3]
    assert len(selected) == 3
    entries = [("fullB_best", None, 0.0, anchor, best_L)]
    for rank, (distance, r) in enumerate(selected, 1):
        entries.append((f"equivalent_{rank}", r, distance,
                        np.asarray(r["best_transformed"], dtype=float),
                        float(r["best_L"])))

    selection = []
    b_points = []
    data = curves()
    with ProcessPoolExecutor(max_workers=6, initializer=init_forward_worker) as pool:
        for label, r, distance, x, recorded_L in entries:
            prediction, diagnostics = predict_parallel(
                pool, 10 ** x[0], x[2], 10 ** x[1])
            res = residuals(data, prediction)
            computed_L = float(res @ res)
            assert math.isclose(computed_L, recorded_L, rel_tol=1e-6,
                                abs_tol=1e-8)
            max_mass = max(d["max_mass_relative_error"] for d in diagnostics.values())
            minimum_u = min(d["minimum_u"] for d in diagnostics.values())
            assert max_mass <= CAP_MASS and minimum_u >= -1e-12
            selection.append({
                "label": label,
                "source": "fullB" if r is None else "E1_profile_point",
                "parameter_index": "" if r is None else r["parameter_index"],
                "position_index": "" if r is None else r["position_index"],
                "log10_alpha": x[0], "log10_Q": x[1], "beta": x[2],
                "distance_from_fullB_transformed": distance,
                "L_recorded": recorded_L, "L_recomputed": computed_L,
                "near_equivalent_ceiling_L": ceiling,
                "max_mass_relative_error": max_mass, "minimum_u": minimum_u,
            })
            for h, points in data.items():
                baseline = points[0][1]
                for i, ((t, observed, curve_id), y) in enumerate(
                        zip(points, prediction[h], strict=True)):
                    b_points.append({
                        "label": label, "h_nm": h, "curve_id": curve_id,
                        "point_index": i, "time_min": t,
                        "observed_increment_uS_cm": observed - baseline,
                        "predicted_increment_uS_cm": y,
                        "prediction_seed_min_uS_cm":
                            (10 ** x[1]) * diagnostics[h]["u_R_seed_min"][i],
                        "prediction_seed_max_uS_cm":
                            (10 ** x[1]) * diagnostics[h]["u_R_seed_max"][i],
                        "max_mass_relative_error":
                            diagnostics[h]["max_mass_relative_error"],
                        "minimum_u": diagnostics[h]["minimum_u"],
                    })
    assert len(selection) == 4 and len(b_points) == 4 * 6 * 7
    write_csv(E1 / "v7_E1_equivalent_selection.csv", selection)
    write_csv(E1 / "v7_E1_equivalent_B_curves.csv", b_points)
    print(f"E1 equivalent B: audited_profile_points={len(valid)}; "
          f"selected={len(selected)}; B_rows={len(b_points)}")


if __name__ == "__main__":
    main()
