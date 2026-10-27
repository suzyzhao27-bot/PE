"""Reaggregate archived v6 held-out residuals under the v7 six-film RMS rule.

This reads frozen v6 outputs; it neither refits old models nor certifies M2.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OLD = ROOT / "原始资料/赵珈绮/14_论文补强学生任务书_v7.0.0/references/v6/WP5"
OUT = ROOT / "WP4/v6_metric_harmonization.csv"


def main() -> None:
    archived = json.loads((OLD / "09_留出门槛核对.json").read_text(encoding="utf-8"))["metrics"]
    rows = []
    for model in ("M0", "M1", "M2"):
        for split, stem, key in (
            ("time", "早段预测逐点", "E_time_macro"),
            ("film", "厚度留出逐点", "E_h_macro"),
        ):
            with (OLD / f"{model}_{stem}.csv").open(encoding="utf-8-sig", newline="") as stream:
                points = list(csv.DictReader(stream))
            per_film = {}
            for h in (65, 120, 200, 350, 700, 2600):
                residuals = [float(p["residual_asinh"]) for p in points
                             if int(p["thickness_nm"]) == h and p["used_for_test"] == "True"]
                assert len(residuals) == (3 if split == "time" else 6), (model, split, h)
                value = math.sqrt(sum(r * r for r in residuals) / len(residuals))
                assert math.isclose(value, archived[model][f"per_film_{'time' if split == 'time' else 'h'}"][str(h)],
                                    rel_tol=1e-9, abs_tol=1e-10), (model, split, h)
                per_film[h] = value
            macro = sum(per_film.values()) / 6
            assert math.isclose(macro, archived[model][key], rel_tol=1e-9, abs_tol=1e-10)
            rows.append({
                "model": model, "split": split,
                "v6_reported_macro_E_asinh": macro,
                "v7_definition_six_film_RMS_E_asinh": math.sqrt(sum(e * e for e in per_film.values()) / 6),
                "n_film": 6, "n_test_points": 18 if split == "time" else 36,
                "source": f"references/v6/WP5/{model}_{stem}.csv",
                "qualification": "historical_M2_technical_conditional" if model == "M2" else "historical_v6",
            })
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print(f"{row['model']} {row['split']}: v6 macro={row['v6_reported_macro_E_asinh']:.6f}; "
              f"v7 RMS={row['v7_definition_six_film_RMS_E_asinh']:.6f}")


if __name__ == "__main__":
    main()
