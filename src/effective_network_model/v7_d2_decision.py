"""Apply the frozen D2 gate once NS full-B and time holdout are complete."""

from __future__ import annotations

import csv
import json
from collections import defaultdict

from v7_d1_objective import OUT, ROOT


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    full = json.loads((OUT / "v7_NS_fullB_summary.json").read_text(encoding="utf-8"))
    time = json.loads((OUT / "v7_NS_time_summary.json").read_text(encoding="utf-8"))
    if full["status"] != "multi_start_stable" or time["status"] != "multi_start_stable":
        raise RuntimeError("NS full/time multi-start status is not stable; D2 remains undecided")
    points = read_csv(OUT / "v7_NS_fullB_points.csv")
    curves = defaultdict(list)
    for row in points:
        if row["model"] == "NS" and row["series"] == "B" and int(row["point_index"]) > 0:
            curves[row["curve_id"]].append(row)
    if len(curves) != 6 or any(len(rows) != 6 for rows in curves.values()):
        raise RuntimeError("Expected six B curves with six nonzero points each")
    windows = []
    for curve_id, rows in curves.items():
        rows.sort(key=lambda r: int(r["point_index"]))
        for left, right in zip(rows, rows[1:]):
            r1, r2 = float(left["asinh_residual"]), float(right["asinh_residual"])
            if r1 * r2 > 0 and abs(r1) > 0.3 and abs(r2) > 0.3:
                windows.append({"curve_id": curve_id,
                                "h_nm": int(left["h_nm"]),
                                "point_indices": [int(left["point_index"]),
                                                  int(right["point_index"])],
                                "time_min": [float(left["time_min"]),
                                             float(right["time_min"])],
                                "asinh_residuals": [r1, r2]})
    residual_curves = sorted(set(row["curve_id"] for row in windows))
    shape_rows = read_csv(ROOT / "WP1/time_shape_fit_summary.csv")
    shape_pe = [r for r in shape_rows if r["role"] in ("B", "L") and
                (float(r["sse_reduction"]) >= 0.30 or
                 float(r["positive_interval_rate_ratio"]) >= 3.0)]
    shape_controls = [{"curve_id": r["curve_id"], "role": r["role"],
                       "sse_reduction": float(r["sse_reduction"]),
                       "positive_interval_rate_ratio": float(
                           r["positive_interval_rate_ratio"])}
                      for r in shape_rows if r["role"] not in ("B", "L")]
    first = len(residual_curves) >= 2 and float(time["E_time"]) > 0.3
    second = len(shape_pe) >= 2
    decision = {"first_condition": first, "second_condition": second,
                "ND_enabled_by_frozen_D2": first and second,
                "NS_full_B_status": full["status"], "NS_time_status": time["status"],
                "NS_time_E": float(time["E_time"]),
                "residual_window_threshold_strict": 0.3,
                "residual_windows": windows,
                "B_curves_with_residual_windows": residual_curves,
                "PE_curves_meeting_shape_condition": [r["curve_id"] for r in shape_pe],
                "control_shape_diagnostics": shape_controls,
                "interpretation_limit": "Curvature is also seen in controls; D2 does not prove PE-specific wetting."}
    (OUT / "v7_D2_decision.json").write_text(
        json.dumps(decision, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"D2: first={first}, second={second}, ND={first and second}; "
          f"B curves with residual windows={len(residual_curves)}, "
          f"NS E_time={time['E_time']:.6g}")


if __name__ == "__main__":
    main()
