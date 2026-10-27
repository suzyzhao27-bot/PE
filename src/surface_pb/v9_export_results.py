"""Export taskbook-required field/observation/layer indexes from saved rows."""

import csv
import hashlib
import json
import math
from dataclasses import asdict
from pathlib import Path

from v9_geometry import make_geometry
from v9_pb_1d import EPS0, F, R, T, Z, reservoir

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / "原始资料/赵珈绮/18_薄化荷电主线任务书_v9.0.0"
if not PACKAGE.is_dir():
    PACKAGE = ROOT / "taskbook"
OUT = ROOT / "WP6"


def read_csv(path):
    with (ROOT / path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path, records):
    with (OUT / path).open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def value(row, key, fallback="NA"):
    item = row.get(key, fallback)
    return fallback if item in (None, "") else item


def main():
    with (PACKAGE / "file_sha256.csv").open(encoding="utf-8-sig", newline="") as stream:
        manifest = list(csv.DictReader(stream))
    spec_hash = next(row["sha256"] for row in manifest if row["path"] == "configs/frozen_spec_v9.json")
    observations = {r["case_id"]: r for r in read_csv("WP1/v9_九组观测与条件盐度.csv")}
    sources = [
        ("WP5/v9_MODEL_ONLY_精简厚度图.csv", "pure_regulated", read_csv("WP5/v9_MODEL_ONLY_精简厚度图.csv")),
        ("WP5/v9_MODEL_ONLY_固定电荷对照.csv", "pure_fixed", read_csv("WP5/v9_MODEL_ONLY_固定电荷对照.csv")),
        ("WP5/v9_MODEL_ONLY_八项敏感性.csv", "family_sensitivity", read_csv("WP5/v9_MODEL_ONLY_八项敏感性.csv")),
    ]
    old = json.loads((PACKAGE / "configs/frozen_spec_v9.json").read_text(encoding="utf-8"))["chemistry"]["legacy_fallback"]
    fields, observed, layers = [], [], []
    for source_file, category, source_rows in sources:
        for row in source_rows:
            if row["status"] != "VALID":
                raise AssertionError(f"cannot export unresolved state {source_file}")
            kind = row["kind"]
            N = int(row["N"])
            h, d, g = float(row["h_nm"]), float(row["d_nm"]), float(row["g_nm"])
            aligned = row.get("aligned_layers", "False") == "True"
            geometry = make_geometry(kind, d, g, N, h, aligned)
            payload = json.dumps(asdict(geometry), ensure_ascii=False, sort_keys=True,
                                 separators=(",", ":"))
            geometry_hash = hashlib.sha256(payload.encode("utf-8")).hexdigest()
            salt = float(row["salt_mM"])
            pH = float(row["pH_measured"] if category == "family_sensitivity" else row["pH"])
            bath = reservoir(salt, pH)
            lam = math.sqrt(EPS0 * 78.3 * R * T / (F * F * float((Z * Z) @ bath))) * 1e9
            mesh_step = min(d / 12, g / 20, lam / 8)
            L = max(5 * geometry.period_nm, 12 * lam)
            case_id = (f"pure_N{N}_c{salt:g}_{kind}"
                       if category == "pure_regulated" else
                       f"fixed_N{N}_c{salt:g}_{kind}"
                       if category == "pure_fixed" else
                       f"{row['case_id']}_{kind}_{row['variation']}")
            parameter_set_id = ("LEGACY_FIXED_FROM_FREE_N20" if category == "pure_fixed"
                                else "LEGACY_REGULATED_MODEL_ONLY")
            field = {"case_id": case_id, "source_csv": source_file,
                     "geometry_hash": geometry_hash, "spec_hash": spec_hash,
                     "free_or_backed": kind, "h_nm": h, "N": N, "d_nm": d, "g_nm": g,
                     "pH": pH, "c_s_mM": salt,
                     "Gamma_nm2": old["Gamma_nm2"] if category != "pure_fixed" else "NA",
                     "pI": old["pI"] if category != "pure_fixed" else "NA",
                     "b": old["b"] if category != "pure_fixed" else "NA",
                     "fixed_sigma_source_C_m2": value(row, "fixed_sigma_source_C_m2"),
                     "parameter_set_id": parameter_set_id,
                     "sigma_mean_C_m2": row["sigma_mean_C_m2"],
                     "sigma_outer_C_m2": row["sigma_outer_C_m2"],
                     "sigma_core_C_m2": row["sigma_core_C_m2"],
                     "Q_proj_C_m2": row["Q_proj_C_m2"],
                     "max_c_M": row["max_species_M"],
                     "Gauss_absolute_C_m": row["Gauss_absolute_C_m"],
                     "Gauss_relative": row["Gauss_relative"],
                     "mesh_target_nm": value(row, "mesh_step_nm", mesh_step),
                     "mesh_L_nm": value(row, "L_nm", L),
                     "convergence_status": row["status"],
                     "physics_validity": "MODEL_ONLY_NUMERICALLY_VALID"}
            fields.append(field)
            parsed_layers = json.loads(row["layer_sigma_json"])
            for layer in geometry.layers:
                layers.append({"case_id": case_id, "geometry_hash": geometry_hash,
                               "layer_id": layer.layer_id,
                               "z_nm": (layer.z_bottom_nm + layer.z_top_nm) / 2,
                               "wet_area_m2_per_m": layer.wet_perimeter_nm * 1e-9,
                               "sigma_mean_C_m2": parsed_layers[str(layer.layer_id)],
                               "is_outer": layer.is_outer, "is_core": layer.is_core})
            if kind == "backed":
                H = float(row["H_um"] if category == "family_sensitivity" else 100)
                zeta = float(row["zeta_app_mV"])
                raw_case_id = row["case_id"] if category == "family_sensitivity" else None
                zeta_obs = (float(observations[raw_case_id]["zeta_reported_mean_mV"])
                            if raw_case_id else None)
                scale = (float(observations[raw_case_id]["comparison_scale_mV"])
                         if raw_case_id else None)
                I_over_Q = -(zeta / 1000) * 12 * EPS0 * 78.3 / (H * 1e-6) ** 2
                observed.append({"case_id": case_id, "H_um": H,
                                 "I_over_Q_C_m3": I_over_Q,
                                 "I_over_Q_source": "algebraically derived from saved zeta_app; not independent fit target",
                                 "zeta_pred_mV": zeta,
                                 "zeta_observed_mV": zeta_obs if zeta_obs is not None else "NA",
                                 "comparison_scale_mV": scale if scale is not None else "NA",
                                 "residual_mV": zeta - zeta_obs if zeta_obs is not None else "NA",
                                 "train_or_predict": "MODEL_ONLY_NOT_VALIDATION",
                                 "calibration_set_id": "NO_COMPATIBLE_SET_LEGACY_ONLY",
                                 "assembly_assumption": "two identical backed layers; full gap H; unverified"})
    if len(fields) != 34 + 24 + 144 or len(layers) <= len(fields):
        raise AssertionError("saved-field export incomplete")
    write_csv("v9_平衡场索引.csv", fields)
    write_csv("v9_观察量索引.csv", observed)
    write_csv("v9_层荷电索引.csv", layers)
    print(f"exported {len(fields)} field rows, {len(observed)} observation rows, {len(layers)} layer rows; no PB solves")


if __name__ == "__main__":
    main()
