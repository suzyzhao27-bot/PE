"""Cost profiles and local dimensionless sensitivities for WP5."""

from __future__ import annotations

import json
import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "effective_network_model"))
from v6_fit_wp4 import bounds, csv_write, load_films, residual, unpack  # noqa: E402
from v6_wp5_holdouts import MANIFEST_SHA256  # noqa: E402

OUT = ROOT / "WP5"
WP4 = ROOT / "WP4"
NAMES = ["ws_ref", "log10_ell_nm", "log10_dm_ratio", "log10_gamma_per_m",
         "log10_q_conductivity", "log10_chi", "log10_tau_s"]


def best_vector(model: str) -> np.ndarray:
    summary = json.loads((WP4 / f"{model}_拟合摘要.json").read_text(encoding="utf-8"))
    searches = json.loads((WP4 / f"{model}_搜索记录.json").read_text(encoding="utf-8"))
    return np.array(next(r["vector"] for r in searches["runs"]
                         if r["start_id"] == summary["best_start_id"]))


def profile(model: str, indices: list[int], films: list) -> list[dict]:
    z0 = best_vector(model)
    lo, hi = bounds(model)
    rows = []
    for fixed in indices:
        free = np.array([j for j in range(len(lo)) if j != fixed])
        for point_id, value in enumerate(np.linspace(lo[fixed], hi[fixed], 15), 1):
            z = z0.copy()
            z[fixed] = value
            def fun(x: np.ndarray) -> np.ndarray:
                full = z.copy()
                full[free] = x
                return residual(full, films, model)
            try:
                opt = least_squares(fun, z0[free], bounds=(lo[free], hi[free]),
                                    max_nfev=100, diff_step=1e-3,
                                    ftol=1e-6, xtol=1e-6, gtol=1e-6)
                z[free] = opt.x
                rows.append({"model_id": model, "scenario_id": "profile", "split_id": "full_data",
                             "data_manifest_hash": MANIFEST_SHA256,
                             "fixed_parameter": NAMES[fixed],
                             "point_id": point_id, "fixed_value_transformed": value,
                             "loss": float(np.sum(opt.fun ** 2)),
                             "status": str(opt.message), "success": bool(opt.success),
                             "nfev": int(opt.nfev),
                             "reoptimized_parameters": json.dumps(unpack(z, model), ensure_ascii=False),
                             "boundary_hits": json.dumps([j for j in range(len(lo)) if
                                                          min((z[j]-lo[j])/(hi[j]-lo[j]),
                                                              (hi[j]-z[j])/(hi[j]-lo[j])) < 1e-3])})
            except Exception as exc:
                rows.append({"model_id": model, "scenario_id": "profile", "split_id": "full_data",
                             "data_manifest_hash": MANIFEST_SHA256,
                             "fixed_parameter": NAMES[fixed],
                             "point_id": point_id, "fixed_value_transformed": value,
                             "loss": "", "status": f"{type(exc).__name__}: {exc}",
                             "success": False, "nfev": "",
                             "reoptimized_parameters": "", "boundary_hits": ""})
        print(model, NAMES[fixed], "profile completed", flush=True)
    return rows


def sensitivity(model: str, films: list) -> dict:
    z0 = best_vector(model)
    lo, hi = bounds(model)
    step = 1e-4
    columns = []
    for i in range(len(z0)):
        plus, minus = z0.copy(), z0.copy()
        plus[i] = min(hi[i], z0[i] + step)
        minus[i] = max(lo[i], z0[i] - step)
        if plus[i] == minus[i]:
            columns.append(np.zeros(len(residual(z0, films, model))))
        else:
            columns.append((residual(plus, films, model) - residual(minus, films, model)) /
                           (plus[i] - minus[i]))
    jac = np.column_stack(columns)
    u, singular, vt = np.linalg.svd(jac, full_matrices=False)
    # Correlation of sensitivities, not a statistical parameter covariance.
    norm = np.linalg.norm(jac, axis=0)
    corr = jac.T @ jac / np.maximum(norm[:, None] * norm[None, :], 1e-300)
    return {"model": model, "parameter_coordinates": NAMES[:len(z0)],
            "coordinate_definition": "linear ws; log10 physical positive parameters",
            "singular_values": singular.tolist(),
            "condition_ratio_largest_to_smallest": float(singular[0] / max(singular[-1], 1e-300)),
            "weakest_direction": vt[-1].tolist(),
            "sensitivity_column_norms": norm.tolist(),
            "sensitivity_correlation": corr.tolist(),
            "scope": "local transformed objective derivatives; not statistical confidence intervals"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", choices=("all", "M1", "M2"), default="all")
    args = parser.parse_args()
    films = load_films()
    started = datetime.now(timezone.utc).isoformat()
    if args.model in ("all", "M1"):
        rows = profile("M1", list(range(5)), films)
        csv_write(OUT / "M1_参数剖面.csv", rows)
    else:
        rows = list(__import__("csv").DictReader((OUT / "M1_参数剖面.csv").open(encoding="utf-8-sig")))
    if args.model in ("all", "M2"):
        rows += profile("M2", [5, 6], films)
    csv_write(OUT / "03_参数辨识与稳健性结果.csv", rows)
    evidence = {"started_utc": started,
                "finished_utc": datetime.now(timezone.utc).isoformat(),
                "M1": sensitivity("M1", films),
                "M2": sensitivity("M2", films),
                "profile_rule": "M1 provisional minimal candidate: all five shared coordinates, 15 points each; M2 capacity/exchange: 15 points each, reoptimize all other shared coordinates from full-fit best as local cost profiles"}
    (OUT / "06_局部灵敏度与奇异值.json").write_text(
        json.dumps(evidence, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({m: evidence[m]["singular_values"] for m in ("M1", "M2")}))


if __name__ == "__main__":
    main()
