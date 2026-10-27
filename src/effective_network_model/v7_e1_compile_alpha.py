"""Compile the nine audited NS log10(alpha) profile points without refitting."""

from __future__ import annotations

import csv
import json
import math

from v7_d1_objective import OUT, ROOT


def main() -> None:
    spec = json.loads((ROOT / "原始资料/赵珈绮/14_论文补强学生任务书_v7.0.0/configs/frozen_spec.json").read_text(encoding="utf-8-sig"))
    rule = spec["matrices"]["E1"]
    assert rule["profile_positions"] == 9
    full = json.loads((OUT / "v7_NS_fullB_summary.json").read_text(encoding="utf-8"))
    best_loss = full["L_recomputed"]
    ceiling = best_loss + max(rule["equivalent_loss_abs"], rule["equivalent_loss_rel"] * best_loss)
    raw = [json.loads(line) for line in (ROOT / "WP5/v7_E1_profile_runs.jsonl").read_text(encoding="utf-8").splitlines()]
    reports = [json.loads(line) for line in (ROOT / "WP5/v7_E1_profile_point_reports.jsonl").read_text(encoding="utf-8").splitlines()
               if json.loads(line)["parameter_index"] == 0]
    assert len(reports) == 9 and sorted(r["position_index"] for r in reports) == list(range(9))
    rows = []
    for report in sorted(reports, key=lambda r: r["position_index"]):
        pos = report["position_index"]
        starts = [r for r in raw if (r["parameter_index"], r["position_index"]) == (0, pos)]
        assert len(starts) == 3 and sorted(r["start_index"] for r in starts) == [0, 1, 2]
        fixed = 5.0 * pos / 8
        assert math.isclose(report["fixed_value_transformed"], fixed, abs_tol=1e-12)
        loss = report["best_L"]
        valid = report["status"] == "multi_start_stable" and loss is not None and math.isclose(
            report["L_recomputed"], loss, rel_tol=1e-6, abs_tol=1e-8) and report["max_mass_relative_error"] <= 1e-8 and report["minimum_u"] >= -1e-12
        x = report["best_transformed"]
        rows.append({
            "position_index": pos, "log10_alpha_fixed": fixed, "alpha_m_inv_fixed": 10 ** fixed,
            "best_L": loss, "delta_L_from_full_best": loss - best_loss if loss is not None else None,
            "near_equivalent_ceiling_L": ceiling, "near_equivalent_technically_valid": valid and loss <= ceiling,
            "point_status": report["status"], "start_statuses": ";".join(r["status"] for r in sorted(starts, key=lambda r: r["start_index"])),
            "actual_objective_calls_total": report["actual_objective_calls_total"],
            "log10_Q_fitted": x[1] if x else None, "beta_fitted": x[2] if x else None,
            "log10_Q_near_bound_1pct": bool(x and min(x[1], 6-x[1]) <= .06),
            "beta_near_bound_1pct": bool(x and min(x[2], 6-x[2]) <= .06),
            "max_mass_relative_error": report["max_mass_relative_error"], "minimum_u": report["minimum_u"],
        })
    target = ROOT / "WP5/v7_E1_alpha_profile.csv"
    with target.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    print(f"{target}: 9 audited points; near-equivalent={sum(r['near_equivalent_technically_valid'] for r in rows)}; ceiling={ceiling:.9f}")


if __name__ == "__main__":
    main()
