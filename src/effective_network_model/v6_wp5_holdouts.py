"""Frozen v6 retrospective time and six-fold thickness holdouts."""

from __future__ import annotations

import argparse
import csv
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
from v6_fit_wp4 import (Film, bounds, csv_write, initial_points, load_films,
                        predict, residual, unpack)  # noqa: E402

OUT = ROOT / "WP5"
PARAM_UNITS = {"ws_ref": "1", "ell_nm": "nm", "dm_ratio": "Dm/Db",
               "gamma_per_m": "m^-1", "q_conductivity": "source conductivity unit",
               "chi": "1", "tau_s": "s"}
MANIFEST_SHA256 = "550cca87c80fda62fefa0b2ac9b14cefc5e8ea71d7bf09c946af00676512f2ce"


def fit(model: str, films: list[Film], max_nfev: int) -> tuple[dict | None, list[dict]]:
    lo, hi = bounds(model)
    runs = []
    best = None
    for i, start in enumerate(initial_points(model), 1):
        item = {"start_id": i, "initial_vector": start.tolist()}
        try:
            opt = least_squares(residual, start, bounds=(lo, hi),
                                args=(films, model), method="trf", max_nfev=max_nfev,
                                diff_step=1e-3, ftol=1e-6, xtol=1e-6, gtol=1e-6)
            item.update(success=bool(opt.success), status=int(opt.status),
                        message=str(opt.message), nfev=int(opt.nfev),
                        loss=float(np.sum(opt.fun ** 2)), vector=opt.x.tolist(),
                        boundary_hits=[j for j in range(len(lo)) if
                                       min((opt.x[j] - lo[j]) / (hi[j] - lo[j]),
                                           (hi[j] - opt.x[j]) / (hi[j] - lo[j])) < 1e-3])
            if best is None or item["loss"] < best["loss"]:
                best = item
        except Exception as exc:
            item.update(success=False, status=-1, message=f"{type(exc).__name__}: {exc}")
        runs.append(item)
    return best, runs


def trimmed(film: Film) -> Film:
    mask = film.times_min <= 60
    return replace(film, times_min=film.times_min[mask], raw=film.raw[mask])


def rmse(film: Film, predicted: np.ndarray, mask: np.ndarray) -> float:
    return float(np.sqrt(np.mean((np.arcsinh(predicted[mask]) -
                                  np.arcsinh(film.increment[mask])) ** 2)))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("model", choices=("M0", "M1", "M2"))
    parser.add_argument("--max-nfev", type=int, default=150)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--rebuild-only", action="store_true",
                        help="recreate CSVs from completed saved split searches without refitting")
    args = parser.parse_args()
    films = load_films()
    started = datetime.now(timezone.utc).isoformat()
    def prior_rows(name):
        path = OUT / name
        if not args.resume or args.rebuild_only or not path.exists():
            return []
        with path.open(encoding="utf-8-sig", newline="") as fh:
            return list(csv.DictReader(fh))
    thickness_rows = prior_rows(f"{args.model}_厚度留出.csv")
    thickness_curves = prior_rows(f"{args.model}_厚度留出逐点.csv")
    time_rows = prior_rows(f"{args.model}_早段预测.csv")
    time_curves = prior_rows(f"{args.model}_早段预测逐点.csv")
    search_path = OUT / f"{args.model}_留出搜索记录.json"
    search_details = (json.loads(search_path.read_text(encoding="utf-8"))["details"]
                      if (args.resume or args.rebuild_only) and search_path.exists() else {})
    for split_kind, hold in [("time", None)] + [("thickness", f.h_nm) for f in films]:
        split_id = "train_0_60" if hold is None else f"hold_{hold}nm"
        if args.resume and not args.rebuild_only and split_id in search_details:
            print(args.model, split_id, "resumed from saved run", flush=True)
            continue
        train = [trimmed(f) for f in films] if hold is None else [f for f in films if f.h_nm != hold]
        if args.rebuild_only:
            saved = search_details[split_id]
            best, runs = saved["best"], saved["runs"]
        else:
            best, runs = fit(args.model, train, args.max_nfev)
            search_details[split_id] = {"best": best, "runs": runs}
        if best is None:
            print(args.model, split_id, "all starts failed", flush=True)
            continue
        p = unpack(np.array(best["vector"]), args.model)
        target = films if hold is None else [next(f for f in films if f.h_nm == hold)]
        for film in target:
            pred, diag = predict(film, p, args.model)
            mask = film.times_min > 60 if hold is None else film.times_min > 0
            e = rmse(film, pred, mask)
            row = {"model_id": args.model, "scenario_id": "main_r0_a25_rv1_dbref",
                   "split_id": split_id, "thickness_nm": film.h_nm,
                   "split_role": "retrospective_future_time" if hold is None else
                                 ("boundary_extrapolation" if hold in (65, 2600) else "interior_interpolation"),
                   "data_manifest_hash": MANIFEST_SHA256,
                   "source_opju_sha256_from_wp2": "c9446a965bac8b875bff135676b59352ae3e57a35d6f0be51ad860c20b5a596d",
                   "parameter_values": json.dumps(p, ensure_ascii=False),
                   "parameter_units": json.dumps(PARAM_UNITS, ensure_ascii=False),
                   "boundary_hits": json.dumps(best["boundary_hits"]),
                   "optimizer_status": best["message"],
                   "numerical_status": "BDF success",
                   "E_train": math.sqrt(best["loss"]),
                   "E_time": e if hold is None else "",
                   "E_h": e if hold is not None else "",
                   "per_film_errors": json.dumps({film.h_nm: e}),
                   "mass_error": diag["mass_error"],
                   "minimum_concentration": diag["minimum_concentration"],
                   "scientific_claim_status": "retrospective internal holdout; not independent validation",
                   "successful_starts": sum(r["success"] for r in runs),
                   "failed_starts": sum(not r["success"] for r in runs),
                   "endpoint_observed": film.increment[-1],
                   "endpoint_predicted": pred[-1],
                   "endpoint_error_factor": max(pred[-1] / film.increment[-1],
                                                film.increment[-1] / max(pred[-1], 1e-300))}
            if hold is None:
                time_rows.append(row)
            else:
                thickness_rows.append(row)
            for j, t in enumerate(film.times_min):
                detail = {"model_id": args.model, "split_id": split_id,
                          "thickness_nm": film.h_nm, "dataset_id": film.dataset_id,
                          "time_min": t, "conductivity_raw": film.raw[j],
                          "baseline_raw": film.raw[0], "increment_observed": film.increment[j],
                          "increment_predicted": pred[j], "residual_raw_unit": pred[j] - film.increment[j],
                          "residual_asinh": np.arcsinh(pred[j]) - np.arcsinh(film.increment[j]),
                          "used_for_train": (t > 0 and (hold is not None and film.h_nm != hold or
                                                             hold is None and t <= 60)),
                          "used_for_test": bool(mask[j])}
                if hold is None:
                    time_curves.append(detail)
                else:
                    thickness_curves.append(detail)
        print(args.model, split_id, "loss", best["loss"],
              "success", sum(r["success"] for r in runs), flush=True)
        (OUT / f"{args.model}_留出搜索记录.json").write_text(
            json.dumps({"model": args.model, "started_utc": started,
                        "updated_utc": datetime.now(timezone.utc).isoformat(),
                        "seed": 20261021, "fixed_starts_per_split": 12,
                        "max_nfev": args.max_nfev, "details": search_details},
                       ensure_ascii=False, indent=2), encoding="utf-8")
    if thickness_rows:
        csv_write(OUT / f"{args.model}_厚度留出.csv", thickness_rows)
        csv_write(OUT / f"{args.model}_厚度留出逐点.csv", thickness_curves)
    if time_rows:
        csv_write(OUT / f"{args.model}_早段预测.csv", time_rows)
        csv_write(OUT / f"{args.model}_早段预测逐点.csv", time_curves)


if __name__ == "__main__":
    main()
