"""Prespecified one-at-a-time WP5 sensitivity matrix for provisional M1."""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "effective_network_model"))
sys.path.insert(0, str(ROOT / "WP2" / "src"))
from transport import D_B_REF, build_operator  # noqa: E402
from v6_fit_wp4 import (bounds, config, csv_write, initial_points,
                        load_films, predict, predict_for_fit, unpack)  # noqa: E402
from v6_wp5_holdouts import MANIFEST_SHA256  # noqa: E402

OUT = ROOT / "WP5"
SEED = 20261021
MODEL = "M1"
SCENARIOS = [
    ("robin_r_0.75", {"robin_r": 0.75}),
    ("a_20nm", {"a_m": 20e-9}), ("a_30nm", {"a_m": 30e-9}),
    ("rv_0.5", {"volume_ratio": 0.5}), ("rv_2", {"volume_ratio": 2.0}),
    ("db_x0.7", {"db_m2_s": 0.7 * D_B_REF}),
    ("db_x1.3", {"db_m2_s": 1.3 * D_B_REF}),
    ("layer_10um", {"boundary_layer_m": 10e-6}),
    ("layer_100um", {"boundary_layer_m": 100e-6}),
    ("preload_cstar_0.25", {"preload_cstar": 0.25}),
]


def initial(film, cfg, scenario):
    if "preload_cstar" not in scenario:
        return 1.0, None, None
    op = build_operator(cfg, 50)
    mobile = scenario["preload_cstar"] * np.exp(-op.field_centers)
    stored = cfg.chi * mobile if cfg.chi > 0 else None
    m_mem = cfg.gamma_per_m * op.dx_m * (np.sum(mobile) + (np.sum(stored) if stored is not None else 0))
    left = 1 - m_mem / cfg.volume_ratio
    if left < 0:
        raise ValueError("preloaded membrane exceeds frozen total salt budget")
    return float(left), mobile, stored


def scenario_prediction(film, p, scenario, fast):
    overrides = {k: v for k, v in scenario.items() if k != "preload_cstar"}
    cfg = replace(config(film, p), **overrides)
    left, mobile, stored = initial(film, cfg, scenario)
    if fast:
        return predict_for_fit(film, p, MODEL, cfg_override=cfg,
                               left=left, mobile=mobile, storage=stored), None
    return predict(film, p, MODEL, cfg_override=cfg,
                   left=left, mobile=mobile, storage=stored)


def residual(z, films, scenario):
    p = unpack(z, MODEL)
    items = []
    for film in films:
        y, _ = scenario_prediction(film, p, scenario, fast=True)
        if not np.all(np.isfinite(y)) or np.min(y) < -1e-7:
            raise FloatingPointError("invalid prediction")
        items.append((np.arcsinh(y[1:]) - np.arcsinh(film.increment[1:])) /
                     math.sqrt(len(films) * (len(film.times_min) - 1)))
    return np.concatenate(items)


def fit(films, scenario, max_nfev):
    lo, hi = bounds(MODEL)
    best = None
    runs = []
    for i, start in enumerate(initial_points(MODEL), 1):
        item = {"start_id": i, "initial_vector": start.tolist()}
        try:
            opt = least_squares(residual, start, bounds=(lo, hi),
                                args=(films, scenario), max_nfev=max_nfev,
                                diff_step=1e-3, ftol=1e-6, xtol=1e-6, gtol=1e-6)
            item.update(success=bool(opt.success), status=int(opt.status),
                        message=str(opt.message), nfev=int(opt.nfev),
                        loss=float(np.sum(opt.fun ** 2)), vector=opt.x.tolist())
            if best is None or item["loss"] < best["loss"]:
                best = item
        except Exception as exc:
            item.update(success=False, status=-1, message=f"{type(exc).__name__}: {exc}")
        runs.append(item)
    return best, runs


def perturbed_films(films, epsilon, sample_id):
    rng = np.random.default_rng(SEED + 1000 * (1 if epsilon < 0.2 else 2) + sample_id)
    return [replace(f, raw=f.raw + rng.uniform(-epsilon, epsilon, len(f.raw))) for f in films]


def directed_films(films, epsilon):
    updated = []
    for f in films:
        if f.h_nm != 200:
            updated.append(f)
            continue
        raw = f.raw.copy()
        for t, sign in ((0, 1), (20, -1), (40, 1), (120, -1)):
            raw[np.where(f.times_min == t)[0][0]] += sign * epsilon
        updated.append(replace(f, raw=raw))
    return updated


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-nfev", type=int, default=120)
    args = parser.parse_args()
    films = load_films()
    tasks = [(name, config_values, films) for name, config_values in SCENARIOS]
    for epsilon in (0.111, 0.25):
        for sample_id in range(16):
            tasks.append((f"reading_{epsilon}_seed{sample_id:02d}", {},
                          perturbed_films(films, epsilon, sample_id)))
        tasks.append((f"reading_{epsilon}_directed200", {},
                      directed_films(films, epsilon)))
    started = datetime.now(timezone.utc).isoformat()
    rows = []
    details = {}
    for task_i, (name, scenario, task_films) in enumerate(tasks, 1):
        best, runs = fit(task_films, scenario, args.max_nfev)
        details[name] = {"scenario": scenario, "best": best, "runs": runs}
        if best is None:
            rows.append({"model_id": MODEL, "scenario_id": name, "split_id": "full_data",
                         "data_manifest_hash": MANIFEST_SHA256,
                         "status": "all starts failed",
                         "loss_bdf": "", "parameters": "", "max_mass_error": "",
                         "minimum_concentration": "", "success_starts": 0,
                         "failed_starts": len(runs), "per_film_rmse_asinh": ""})
        else:
            p = unpack(np.array(best["vector"]), MODEL)
            per_film = {}
            mass = []
            minimum = []
            for film in task_films:
                y, diag = scenario_prediction(film, p, scenario, fast=False)
                per_film[film.h_nm] = float(np.sqrt(np.mean((np.arcsinh(y[1:]) -
                                                              np.arcsinh(film.increment[1:])) ** 2)))
                mass.append(diag["mass_error"])
                minimum.append(diag["minimum_concentration"])
            rows.append({"model_id": MODEL, "scenario_id": name, "split_id": "full_data",
                         "data_manifest_hash": MANIFEST_SHA256,
                         "status": best["message"], "loss_bdf": float(np.mean(np.square(list(per_film.values())))),
                         "parameters": json.dumps(p, ensure_ascii=False),
                         "max_mass_error": max(mass), "minimum_concentration": min(minimum),
                         "success_starts": sum(r["success"] for r in runs),
                         "failed_starts": sum(not r["success"] for r in runs),
                         "per_film_rmse_asinh": json.dumps(per_film)})
        csv_write(OUT / "07_固定情景与读数扰动.csv", rows)
        (OUT / "08_敏感性搜索记录.json").write_text(json.dumps(
            {"started_utc": started, "updated_utc": datetime.now(timezone.utc).isoformat(),
             "candidate": MODEL, "seed": SEED, "starts_per_fit": 12,
             "max_nfev": args.max_nfev, "preload_cstar": 0.25,
             "tasks": details}, ensure_ascii=False, indent=2), encoding="utf-8")
        print(task_i, "/", len(tasks), name, rows[-1]["loss_bdf"], flush=True)


if __name__ == "__main__":
    main()
