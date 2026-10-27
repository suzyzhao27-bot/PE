"""Compile frozen L early-fit/late-prediction records without rerunning models."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "WP4"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    rows = []
    observed_identity = None
    for model in ("N0", "NS", "ND"):
        summary = json.loads((OUT / f"v7_{model}_L_summary.json").read_text(encoding="utf-8"))
        full_b = json.loads((OUT / f"v7_{model}_fullB_summary.json").read_text(encoding="utf-8"))
        points = read_csv(OUT / f"v7_{model}_L_points.csv")
        assert summary["status"] == "multi_start_stable" and len(points) == 70, model
        assert len(summary["per_membrane"]) == 3, model
        assert full_b["status"] == "multi_start_stable", model
        fixed = [full_b["beta"]]
        if model == "ND":
            fixed.extend((full_b["tau0_s"], full_b["b"]))
            assert all(math.isclose(a, b, rel_tol=1e-12) for a, b in
                       zip(fixed, summary["fixed_beta_tau0_b_from_full_B"], strict=True))
        else:
            assert math.isclose(fixed[0], summary["fixed_beta_from_full_B"], rel_tol=1e-12)
        identity = [(p["h_nm"], p["curve_id"], p["point_index"], p["time_min"],
                     p["role"], p["conductivity_uS_cm"], p["baseline_uS_cm"],
                     p["observed_increment_uS_cm"]) for p in points]
        if observed_identity is None:
            observed_identity = identity
        else:
            assert identity == observed_identity, f"L observations/roles changed for {model}"
        per_film = []
        for item in summary["per_membrane"]:
            h = int(item["h_nm"])
            early = [float(p["asinh_residual"]) for p in points
                     if int(p["h_nm"]) == h and p["role"] == "early_train"]
            late = [float(p["asinh_residual"]) for p in points
                    if int(p["h_nm"]) == h and p["role"] == "late_prediction"]
            assert len(early) == item["n_early_train"] and len(late) == item["n_late_prediction"]
            assert all(float(p["time_min"]) <= 240 for p in points
                       if int(p["h_nm"]) == h and p["role"] == "early_train")
            assert all(float(p["time_min"]) > 240 for p in points
                       if int(p["h_nm"]) == h and p["role"] == "late_prediction")
            assert math.isclose(math.sqrt(sum(r*r for r in early)/len(early)), item["E_early"], rel_tol=1e-10)
            assert math.isclose(math.sqrt(sum(r*r for r in late)/len(late)), item["E_late"], rel_tol=1e-10)
            per_film.append(item["E_late"])
        assert math.isclose(math.sqrt(sum(e*e for e in per_film)/3), summary["E_Llate"], rel_tol=1e-10)
        for p in points:
            rows.append({
                "model": model, "h_nm": p["h_nm"], "curve_id": p["curve_id"],
                "point_index": p["point_index"], "time_min": p["time_min"],
                "role": p["role"], "conductivity_uS_cm": p["conductivity_uS_cm"],
                "baseline_uS_cm": p["baseline_uS_cm"],
                "observed_increment_uS_cm": p["observed_increment_uS_cm"],
                "predicted_increment_uS_cm": p["predicted_increment_uS_cm"],
                "original_residual_uS_cm": float(p["predicted_increment_uS_cm"]) -
                                           float(p["observed_increment_uS_cm"]),
                "asinh_residual": p["asinh_residual"],
                "alpha_L_m_inv": summary["alpha_L_m_inv"],
                "Q_L_uS_cm_per_unit_uR": summary["Q_L_uS_cm_per_unit_uR"],
                "B_fixed_beta": fixed[0],
                "B_fixed_tau0_s": fixed[1] if model == "ND" else "",
                "B_fixed_b": fixed[2] if model == "ND" else "",
                "E_Llate_asinh": summary["E_Llate"],
                "fit_status": summary["status"],
                "source_points": f"WP4/v7_{model}_L_points.csv",
            })
    assert len(rows) == 210
    with (OUT / "cross_series_transfer.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    print(f"L conditional transfer: {len(rows)} rows; 3 models; early<=240 min fit, late>240 min prediction")
    for model in ("N0", "NS", "ND"):
        print(f"{model}: E_Llate={next(r['E_Llate_asinh'] for r in rows if r['model'] == model):.9f}")


if __name__ == "__main__":
    main()
