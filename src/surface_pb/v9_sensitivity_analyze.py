"""Compare the frozen MODEL_ONLY sensitivity matrix without refitting."""

import csv
import json
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "WP5/v9_MODEL_ONLY_八项敏感性.csv"
DELTA = ROOT / "WP5/v9_MODEL_ONLY_八项敏感性差值.csv"
SUMMARY = ROOT / "WP5/v9_MODEL_ONLY_八项敏感性判定.json"


def sign(value, tolerance=1e-12):
    return 1 if value > tolerance else -1 if value < -tolerance else 0


def main():
    with RAW.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 144:
        raise AssertionError("frozen 9-condition eight-sensitivity matrix incomplete")
    lookup = {(r["case_id"], r["kind"], r["variation"]): r for r in rows}
    if len(lookup) != len(rows):
        raise AssertionError("duplicate sensitivity state")
    if any(r["status"] != "VALID" for r in rows):
        raise AssertionError("unresolved states require disposition before comparison")
    variations = ["g_x0.5", "g_x2", "H_x0.9", "H_x1.1", "cs_x0.8",
                  "cs_x1.2", "PE_eps_2.3", "aligned_layers"]
    deltas = []
    for r in rows:
        if r["variation"] == "baseline":
            continue
        base = lookup[(r["case_id"], r["kind"], "baseline")]
        base_sigma, sigma = float(base["sigma_mean_C_m2"]), float(r["sigma_mean_C_m2"])
        zeta_delta = (float(r["zeta_app_mV"]) - float(base["zeta_app_mV"])
                      if r["kind"] == "backed" else None)
        deltas.append({"case_id": r["case_id"], "kind": r["kind"],
                       "variation": r["variation"],
                       "sigma_base_C_m2": base_sigma, "sigma_varied_C_m2": sigma,
                       "delta_sigma_C_m2": sigma - base_sigma,
                       "abs_delta_sigma_relative": abs(sigma - base_sigma) / max(abs(base_sigma), 1e-30),
                       "zeta_base_mV": float(base["zeta_app_mV"]) if r["kind"] == "backed" else "NA",
                       "zeta_varied_mV": float(r["zeta_app_mV"]) if r["kind"] == "backed" else "NA",
                       "delta_zeta_mV": zeta_delta if zeta_delta is not None else "NA",
                       "N_base": int(base["N"]), "N_varied": int(r["N"]),
                       "g_base_nm": float(base["g_nm"]), "g_varied_nm": float(r["g_nm"]),
                       "salt_base_mM": float(base["salt_mM"]), "salt_varied_mM": float(r["salt_mM"])})
    with DELTA.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(deltas[0]))
        writer.writeheader()
        writer.writerows(deltas)
    aggregate = {}
    for variation in variations:
        group = [r for r in deltas if r["variation"] == variation]
        aggregate[variation] = {
            "rows": len(group),
            "max_abs_delta_sigma_relative": max(r["abs_delta_sigma_relative"] for r in group),
            "max_abs_delta_zeta_mV": max(abs(r["delta_zeta_mV"]) for r in group
                                           if r["delta_zeta_mV"] != "NA"),
            "changed_layer_count": sum(r["N_base"] != r["N_varied"] for r in group),
        }
    if any(r["abs_delta_sigma_relative"] != 0 for r in deltas
           if r["variation"].startswith("H_")):
        raise AssertionError("H-only variation changed electrostatic field")

    trends = {}
    for variation in ("baseline", "g_x0.5", "g_x2", "cs_x0.8", "cs_x1.2",
                      "PE_eps_2.3", "aligned_layers"):
        neutral182_free = lookup[("182_pH7", "free", variation)]
        neutral604_free = lookup[("604_pH7", "free", variation)]
        neutral182_backed = lookup[("182_pH7", "backed", variation)]
        neutral604_backed = lookup[("604_pH7", "backed", variation)]
        sigma_contrast = (abs(float(neutral604_free["sigma_mean_C_m2"]))
                          - abs(float(neutral182_free["sigma_mean_C_m2"])))
        zeta_contrast = (float(neutral604_backed["zeta_app_mV"])
                         - float(neutral182_backed["zeta_app_mV"]))
        trends[variation] = {"free_abs_sigma_604_minus182_C_m2": sigma_contrast,
                             "backed_zeta_604_minus182_mV": zeta_contrast,
                             "free_abs_sigma_sign": sign(sigma_contrast),
                             "backed_zeta_sign": sign(zeta_contrast)}
    base_trend = trends["baseline"]
    for variation, record in trends.items():
        record["free_sigma_direction_changed_from_baseline"] = (
            record["free_abs_sigma_sign"] != base_trend["free_abs_sigma_sign"])
        record["backed_zeta_direction_changed_from_baseline"] = (
            record["backed_zeta_sign"] != base_trend["backed_zeta_sign"])
    dielectric = aggregate["PE_eps_2.3"]
    dielectric["exceeds_frozen_5mV_or_10percent"] = (
        dielectric["max_abs_delta_zeta_mV"] > 5 or
        dielectric["max_abs_delta_sigma_relative"] > 0.1)
    outcome = {"parameter_role": "MODEL_ONLY_LEGACY_REFERENCE_NOT_PE_FIT",
               "raw_rows": len(rows), "delta_rows": len(deltas),
               "all_numeric_valid": True, "aggregates": aggregate,
               "neutral_604_minus182_proxy_trends": trends,
               "dielectric_physics_unresolved_by_frozen_threshold":
                   dielectric["exceeds_frozen_5mV_or_10percent"],
               "material_specific_support": False}
    SUMMARY.write_text(json.dumps(outcome, ensure_ascii=False, indent=2) + "\n",
                       encoding="utf-8")
    print(json.dumps({"dielectric": dielectric, "trends": trends,
                      "raw_rows": len(rows), "delta_rows": len(deltas)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
