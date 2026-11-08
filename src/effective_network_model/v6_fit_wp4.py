"""Frozen v6 WP4 equal-thickness full-curve M0/M1/M2 fitting.

Run from the repository root. All data remain in their original units;
conductivity is a shared conditional scale Q times the right-pool fraction.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import platform
import sys
from collections import defaultdict
from dataclasses import asdict, dataclass, replace
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy
from scipy.linalg import eigh, expm
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "WP2" / "src"))
from transport import Config, D_B_REF, build_operator, m0_right_fraction, solve  # noqa: E402

MANIFEST = ROOT / "WP1" / "03_数据与参数清单.csv"
OUT = ROOT / "WP4"
THICKNESSES = (65, 120, 200, 350, 700, 2600)
SOURCE_SHA256 = "c9446a965bac8b875bff135676b59352ae3e57a35d6f0be51ad860c20b5a596d"
SEED = 20261021
N_STARTS = 12
N_CELLS = 50
RTOL = 1e-7
ATOL = 1e-10


@dataclass(frozen=True)
class Film:
    h_nm: int
    dataset_id: str
    source_file: str
    time_column: str
    value_column: str
    unit_status: str
    times_min: np.ndarray
    raw: np.ndarray

    @property
    def increment(self) -> np.ndarray:
        return self.raw - self.raw[0]


def load_films() -> list[Film]:
    groups: dict[int, list[dict[str, str]]] = defaultdict(list)
    with MANIFEST.open(encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh):
            if row["record_type"] == "observation" and row["use_role"] == "primary_conditional_fit":
                groups[int(float(row["thickness_nm"]))].append(row)
    if tuple(sorted(groups)) != THICKNESSES:
        raise ValueError(f"unexpected thickness set: {sorted(groups)}")
    films = []
    for h in THICKNESSES:
        rows = sorted(groups[h], key=lambda r: float(r["time_min"]))
        times = np.array([float(r["time_min"]) for r in rows])
        raw = np.array([float(r["conductivity_raw"]) for r in rows])
        if len(rows) != 7 or times[0] != 0 or times[-1] != 120 or not np.all(np.diff(times) > 0):
            raise ValueError(f"bad time sequence at {h} nm")
        if times[4] != (85 if h == 65 else 80):
            raise ValueError(f"wrong fifth sample time at {h} nm")
        keys = ("dataset_id", "source_file", "time_column", "value_column", "conductivity_unit_status")
        if any(len({r[k] for r in rows}) != 1 for k in keys):
            raise ValueError(f"inconsistent provenance at {h} nm")
        r = rows[0]
        films.append(Film(h, r["dataset_id"], r["source_file"], r["time_column"],
                          r["value_column"], r["conductivity_unit_status"], times, raw))
    return films


def bounds(model: str) -> tuple[np.ndarray, np.ndarray]:
    lo = [0, 1, -12, -2, 0]
    hi = [20, 4, 0, 6, 8]
    if model == "M2":
        lo += [-4, -3]
        hi += [6, 6]
    return np.asarray(lo, float), np.asarray(hi, float)


def unpack(z: np.ndarray, model: str) -> dict[str, float]:
    out = {"ws_ref": float(z[0]), "ell_nm": float(10 ** z[1]),
           "dm_ratio": float(10 ** z[2]), "gamma_per_m": float(10 ** z[3]),
           "q_conductivity": float(10 ** z[4])}
    if model == "M2":
        out.update(chi=float(10 ** z[5]), tau_s=float(10 ** z[6]))
    else:
        out.update(chi=0.0, tau_s=1.0)
    return out


def config(film: Film, p: dict[str, float]) -> Config:
    return Config(h_m=film.h_nm * 1e-9, a_m=25e-9, ell_m=p["ell_nm"] * 1e-9,
                  ws_ref=p["ws_ref"], robin_r=0.0,
                  dm_m2_s=p["dm_ratio"] * D_B_REF, db_m2_s=D_B_REF,
                  gamma_per_m=p["gamma_per_m"], volume_ratio=1.0,
                  chi=p["chi"], tau_s=p["tau_s"], q_conductivity=p["q_conductivity"])


def predict(film: Film, p: dict[str, float], model: str,
            n_cells: int = N_CELLS, *, cfg_override: Config | None = None,
            left: float = 1.0, mobile: np.ndarray | None = None,
            storage: np.ndarray | None = None) -> tuple[np.ndarray, dict]:
    cfg = cfg_override or config(film, p)
    times = film.times_min * 60
    if model == "M0":
        if left != 1.0 or mobile is not None or storage is not None:
            raise ValueError("M0 has no preloaded mobile or storage state")
        y = p["q_conductivity"] * m0_right_fraction(times, cfg)
        return y, {"mass_error": 0.0, "minimum_concentration": 0.0,
                   "left_final": float(1 - y[-1] / p["q_conductivity"]),
                   "right_final": float(y[-1] / p["q_conductivity"]),
                   "mobile_final": None, "storage_final": None, "nfev": 0}
    traj = solve(cfg, times, n_cells=n_cells, left=left, mobile=mobile,
                 storage=storage, rtol=RTOL, atol=ATOL)
    op = traj.operator
    n = op.n_cells
    last = traj.states[:, -1]
    return traj.conductivity_increment, {
        "mass_error": traj.mass_error,
        "minimum_concentration": float(np.min(traj.states)),
        "left_final": float(last[0]), "right_final": float(last[n + 1]),
        "mobile_final": float(cfg.gamma_per_m * op.dx_m * np.sum(last[1:n + 1])),
        "storage_final": float(cfg.gamma_per_m * op.dx_m * np.sum(last[n + 2:])) if op.has_storage else 0.0,
        "nfev": traj.solver_nfev,
    }


def predict_for_fit(film: Film, p: dict[str, float], model: str,
                    *, cfg_override: Config | None = None,
                    left: float = 1.0, mobile: np.ndarray | None = None,
                    storage: np.ndarray | None = None,
                    n_cells: int = N_CELLS) -> np.ndarray:
    """Reversible similarity transform and symmetric eigensolution for fitting.

    The discrete chain with storage leaves has detailed balance on every edge.
    Final exported trajectories independently use the validated BDF solver.
    """
    if model == "M0":
        return predict(film, p, model, cfg_override=cfg_override,
                       left=left, mobile=mobile, storage=storage,
                       n_cells=n_cells)[0]
    cfg = cfg_override or config(film, p)
    rapid_limit = (model == "M2" and cfg.chi / cfg.tau_s > 1e6 and
                   cfg.tau_s <= 0.1 and mobile is None and storage is None)
    if rapid_limit:
        # The specified M2 fast-exchange limit: (1+chi) dc/dt = -div(j).
        # Used only to screen very stiff optimizer points with empty membranes;
        # accepted outputs are recomputed with the full M2 BDF solver.
        op = build_operator(replace(cfg, chi=0.0), n_cells)
        a = op.matrix.toarray()
        a[1:op.n_cells + 1, :] /= 1 + cfg.chi
    else:
        op = build_operator(cfg, n_cells)
        a = op.matrix.toarray()
    if model == "M2" and cfg.chi / cfg.tau_s > 1e6 and not rapid_limit:
        # Extreme exchange separates eigenvalues by many orders of magnitude;
        # dense matrix exponentials avoid cancellation in small early signals.
        y0 = op.initial(left=left, mobile=mobile, storage=storage)
        right_idx = op.n_cells + 1
        out = np.array([(expm(a * (t * 60)) @ y0)[right_idx] - y0[right_idx]
                        for t in film.times_min])
        out[0] = 0.0
        return p["q_conductivity"] * out
    d = np.ones(a.shape[0])
    for j in range(1, op.n_cells + 2):
        d[j] = d[j - 1] * math.sqrt(a[j, j - 1] / a[j - 1, j])
    if op.has_storage:
        for j in range(op.n_cells):
            c, s = 1 + j, op.n_cells + 2 + j
            d[s] = d[c] * math.sqrt(a[s, c] / a[c, s])
    symmetric = (a / d[:, None]) * d[None, :]
    values, vectors = eigh(symmetric, check_finite=False)
    t = film.times_min * 60
    y0 = op.initial(left=left, mobile=mobile, storage=storage)
    coeff = vectors.T @ (y0 / d)
    right = d[op.n_cells + 1] * ((vectors[op.n_cells + 1] * coeff) @
                                   np.exp(values[:, None] * t[None, :]))
    right -= y0[op.n_cells + 1]
    # Initial right-pool increment is exactly zero; remove eigensum cancellation.
    right[0] = 0.0
    if not np.all(np.isfinite(right)) or np.min(right) < -1e-8:
        return predict(film, p, model, cfg_override=cfg,
                       left=left, mobile=mobile, storage=storage,
                       n_cells=n_cells)[0]
    return p["q_conductivity"] * right


def residual(z: np.ndarray, films: list[Film], model: str, n_cells: int = N_CELLS) -> np.ndarray:
    p = unpack(z, model)
    pieces = []
    for film in films:
        pred = predict_for_fit(film, p, model, n_cells=n_cells)
        if not np.all(np.isfinite(pred)) or np.min(pred) < -1e-7:
            raise FloatingPointError("nonfinite or negative prediction")
        pieces.append((np.arcsinh(pred[1:]) - np.arcsinh(film.increment[1:])) /
                      math.sqrt(len(films) * (len(film.times_min) - 1)))
    return np.concatenate(pieces)


def initial_points(model: str) -> list[np.ndarray]:
    lo, hi = bounds(model)
    rng = np.random.default_rng(SEED)
    # Fixed, reusable Latin-hypercube starts, independent of all fitted results.
    design = np.empty((N_STARTS, len(lo)))
    for j in range(len(lo)):
        design[:, j] = (rng.permutation(N_STARTS) + 0.5) / N_STARTS
    points = lo + design * (hi - lo)
    return [points[i] for i in range(N_STARTS)]


def csv_write(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("model", choices=("M0", "M1", "M2"))
    parser.add_argument("--max-nfev", type=int, default=150)
    args = parser.parse_args()
    films = load_films()
    lo, hi = bounds(args.model)
    started = datetime.now(timezone.utc).isoformat()
    runs = []
    best = None
    for i, start in enumerate(initial_points(args.model), 1):
        record = {"start_id": i, "initial_vector": start.tolist()}
        try:
            fit = least_squares(residual, start, bounds=(lo, hi),
                                args=(films, args.model), method="trf",
                                max_nfev=args.max_nfev, diff_step=1e-3,
                                ftol=1e-6, xtol=1e-6, gtol=1e-6)
            record.update(success=bool(fit.success), status=int(fit.status),
                          message=str(fit.message), nfev=int(fit.nfev),
                          loss=float(np.sum(fit.fun ** 2)), vector=fit.x.tolist(),
                          boundary_hits=[j for j in range(len(lo)) if
                                         min((fit.x[j] - lo[j]) / (hi[j] - lo[j]),
                                             (hi[j] - fit.x[j]) / (hi[j] - lo[j])) < 1e-3])
            if best is None or record["loss"] < best["loss"]:
                best = record
        except Exception as exc:
            record.update(success=False, status=-1, message=f"{type(exc).__name__}: {exc}")
        runs.append(record)
        print(f"{args.model} start {i}/{N_STARTS}: {record.get('loss', 'FAIL')} {record['message']}", flush=True)
        (OUT / f"{args.model}_搜索记录.json").write_text(
            json.dumps({"model": args.model, "started_utc": started,
                        "updated_utc": datetime.now(timezone.utc).isoformat(),
                        "settings": {"seed": SEED, "n_starts": N_STARTS, "n_cells": N_CELLS,
                                     "rtol": RTOL, "atol": ATOL, "max_nfev": args.max_nfev},
                        "runs": runs}, ensure_ascii=False, indent=2), encoding="utf-8")
    if best is None:
        raise RuntimeError(f"all {args.model} starts failed")
    p = unpack(np.array(best["vector"]), args.model)
    obs_rows = []
    film_rows = []
    for film in films:
        pred, info = predict(film, p, args.model)
        obs = film.increment
        transformed = np.arcsinh(pred) - np.arcsinh(obs)
        for j, t in enumerate(film.times_min):
            obs_rows.append({"model_id": args.model, "scenario_id": "main_r0_a25_rv1_dbref",
                             "dataset_id": film.dataset_id, "source_file": film.source_file,
                             "thickness_nm": film.h_nm, "time_column": film.time_column,
                             "value_column": film.value_column, "time_min": t,
                             "conductivity_unit_status": film.unit_status,
                             "conductivity_raw": film.raw[j], "baseline_raw": film.raw[0],
                             "increment_observed": obs[j], "increment_predicted": pred[j],
                             "residual_raw_unit": pred[j] - obs[j],
                             "residual_asinh": transformed[j], "used_for_fit": j > 0})
        film_rows.append({"model_id": args.model, "thickness_nm": film.h_nm,
                          "n_fit_points": 6,
                          "rmse_asinh": float(np.sqrt(np.mean(transformed[1:] ** 2))),
                          "rmse_raw_unit": float(np.sqrt(np.mean((pred[1:] - obs[1:]) ** 2))),
                          "endpoint_observed": float(obs[-1]),
                          "endpoint_predicted": float(pred[-1]),
                          "endpoint_pred_over_obs": float(pred[-1] / obs[-1]),
                          **info})
    csv_write(OUT / f"{args.model}_观测预测.csv", obs_rows)
    csv_write(OUT / f"{args.model}_逐膜残差.csv", film_rows)
    summary = {"model_id": args.model, "started_utc": started,
               "finished_utc": datetime.now(timezone.utc).isoformat(),
               "actual_executor": "Codex on local workspace; student independent execution not recorded",
               "python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
               "command": f".venv/Scripts/python.exe src/effective_network_model/v6_fit_wp4.py {args.model} --max-nfev {args.max_nfev}",
               "manifest": str(MANIFEST.relative_to(ROOT)),
               "raw_source_sha256_reused_from_wp2": SOURCE_SHA256,
               "data_scope": "Book1 six PE curves, 42 raw points, 36 nonzero fitted points",
               "objective": "equal-thickness mean of six-point asinh transformed squared residuals; s0=1 raw unit",
               "best_start_id": best["start_id"], "best_loss": best["loss"],
               "best_success": best["success"], "best_message": best["message"],
               "boundary_hits": best["boundary_hits"], "parameters": p,
               "parameter_units": {"ws_ref": "1", "ell_nm": "nm", "dm_ratio": "Dm/Db",
                                   "gamma_per_m": "m^-1", "q_conductivity": "raw conductivity unit",
                                   "chi": "1", "tau_s": "s"},
               "film_results": film_rows,
               "start_status": [{k: v for k, v in r.items() if k not in ("vector", "initial_vector")}
                                for r in runs]}
    (OUT / f"{args.model}_拟合摘要.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"model": args.model, "loss": best["loss"], "parameters": p,
                      "boundary_hits": best["boundary_hits"]}, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
