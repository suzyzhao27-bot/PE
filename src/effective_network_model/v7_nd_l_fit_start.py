"""Frozen ND L conditional transfer: alpha/Q fit, B beta/tau0/b fixed."""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.optimize import least_squares

from v7_d1_objective import OUT, TASK, graphs
from v7_dynamic_network import dynamic_closed_baths, dynamic_system
from v7_nd_objective import starts_nd
from v7_network_transport import SPEC_PATH


H_L = (200, 350, 700)
_SPEC = None
_DATA = None
_GRAPHS = None
_MATERIAL = None


class BudgetReached(Exception):
    pass


def curves_l():
    with (TASK / "inputs/salt/curve_roles_v7.csv").open(
            encoding="utf-8-sig", newline="") as stream:
        roles = {r["curve_id"]: r["role"] for r in csv.DictReader(stream)}
    by_h = {h: [] for h in H_L}
    with (TASK / "inputs/salt/salt_transport_tidy.csv").open(
            encoding="utf-8-sig", newline="") as stream:
        for row in csv.DictReader(stream):
            if roles[row["curve_id"]] != "L":
                continue
            assert row["time_unit"] == "min" and row["conductivity_unit"] == "μS/cm"
            h = int(row["sample_id"].split()[0])
            by_h[h].append((float(row["time"]), float(row["conductivity"]),
                            row["curve_id"]))
    for h, points in by_h.items():
        points.sort()
        assert points[0][0] == 0 and len({p[0] for p in points}) == len(points)
        assert sum(0 < t <= 240 for t, _, _ in points) >= 3
        assert sum(t > 240 for t, _, _ in points) >= 1
    assert {h: len(p) for h, p in by_h.items()} == {200: 19, 350: 15, 700: 36}
    return by_h


def init_worker(material):
    global _SPEC, _DATA, _GRAPHS, _MATERIAL
    _SPEC = json.loads(SPEC_PATH.read_text(encoding="utf-8-sig"))
    _DATA = curves_l()
    _GRAPHS = {(h, g["seed"]): g for h, group in graphs().items()
               for g in group if h in H_L}
    _MATERIAL = material


def forward_one(task):
    h, seed, alpha = task
    times = np.array([t for t, _, _ in _DATA[h] if t <= 240]) * 60.0
    beta, tau0, b = _MATERIAL
    capacities, edges, _, dynamic_edges = dynamic_system(
        _GRAPHS[h, seed], beta, alpha, tau0, b, _SPEC)
    run = dynamic_closed_baths(capacities, edges, dynamic_edges, times)
    return h, seed, run["u_R"], run["mass_relative_error"], run["minimum_u"]


def predict_early(pool, alpha, q):
    grouped = {h: [] for h in H_L}
    tasks = [(h, seed, alpha) for h in H_L for seed in range(7001, 7009)]
    for h, seed, right, mass_error, minimum in pool.map(forward_one, tasks):
        grouped[h].append((seed, right, mass_error, minimum))
    predictions, diagnostics = {}, {}
    for h in H_L:
        group = sorted(grouped[h], key=lambda x: x[0])
        assert [r[0] for r in group] == list(range(7001, 7009))
        predictions[h] = q * np.mean(np.vstack([r[1] for r in group]), axis=0)
        diagnostics[h] = (max(r[2] for r in group), min(r[3] for r in group))
    return predictions, diagnostics


def residuals_early(data, predictions):
    residual = []
    for h in H_L:
        points = [(t, k) for t, k, _ in data[h] if t <= 240]
        baseline = points[0][1]
        obs = np.array([k - baseline for _, k in points[1:]])
        pred = predictions[h][1:]
        residual.extend((np.arcsinh(pred) - np.arcsinh(obs)) /
                        math.sqrt(len(obs) * len(H_L)))
    result = np.array(residual)
    assert len(result) == 27 and np.all(np.isfinite(result))
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("start_index", type=int)
    args = parser.parse_args()
    if not 0 <= args.start_index < 4:
        raise ValueError("L start_index must be 0..3")
    full = json.loads((OUT / "v7_ND_fullB_summary.json").read_text(
        encoding="utf-8"))
    if not all(key in full for key in ("beta", "tau0_s", "b")):
        raise RuntimeError("ND full-B report has no usable material point")
    material = (float(full["beta"]), float(full["tau0_s"]), float(full["b"]))
    out = OUT / "v7_d1_fit_ND_L_runs.jsonl"
    if out.exists():
        for line in out.read_text(encoding="utf-8").splitlines():
            saved = json.loads(line)
            if (saved["model"], saved["start_index"]) == (
                    "ND", args.start_index):
                raise RuntimeError("This frozen L start already has a recorded result")
    data = curves_l()
    x0 = starts_nd()[args.start_index, :2]
    lower, upper = np.array([0.0, 0.0]), np.array([5.0, 6.0])
    count = 0
    best = {"L": float("inf"), "x": None}
    started = time.perf_counter()
    pool = ProcessPoolExecutor(max_workers=6, initializer=init_worker,
                               initargs=(material,))

    def objective(x):
        nonlocal count
        if count >= 200:
            raise BudgetReached()
        count += 1
        predicted, diagnostics = predict_early(pool, 10 ** x[0], 10 ** x[1])
        if max(v[0] for v in diagnostics.values()) > 1e-8 or \
           min(v[1] for v in diagnostics.values()) < -1e-12:
            raise RuntimeError("conservation or nonnegative state check failed")
        residual = residuals_early(data, predicted)
        loss = float(residual @ residual)
        if loss < best["L"]:
            best.update(L=loss, x=x.tolist())
        if count % 10 == 0:
            print(f"ND L start {args.start_index}: calls {count}/200; "
                  f"best L={best['L']:.6g}", flush=True)
        return residual

    status, message, optimizer = "", "", {}
    try:
        result = least_squares(objective, x0, bounds=(lower, upper), method="trf",
                               jac="2-point", x_scale=1.0, ftol=1e-8, xtol=1e-8,
                               gtol=1e-8, max_nfev=200)
        optimizer = {"scipy_success": bool(result.success),
                     "scipy_status": int(result.status),
                     "scipy_nfev": int(result.nfev), "message": str(result.message)}
        status = "solver_success" if result.success else "max_eval_or_solver_stop"
        message = str(result.message)
    except BudgetReached:
        status = "max_eval"
        message = "200 actual residual calls reached, including finite differences"
    except Exception as exc:
        status = "FAIL"
        message = repr(exc)
    finally:
        pool.shutdown(wait=True)
    record = {"model": "ND", "task": "L_early_le_240_min",
              "start_index": args.start_index, "start_transformed": x0.tolist(),
              "bounds_transformed": {"lower": lower.tolist(), "upper": upper.tolist()},
              "fixed_beta_tau0_b_from_full_B": list(material),
              "fixed_full_B_status": full["status"],
              "settings": {"method": "trf", "jac": "2-point", "x_scale": 1.0,
                           "ftol": 1e-8, "xtol": 1e-8, "gtol": 1e-8,
                           "actual_objective_call_limit": 200, "parallel_workers": 6},
              "status": status, "message": message, "actual_objective_calls": count,
              "best_transformed": best["x"],
              "best_L": best["L"] if np.isfinite(best["L"]) else None,
              "elapsed_s": time.perf_counter() - started, "optimizer": optimizer}
    with out.open("a", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(record, ensure_ascii=False) + "\n")
    print(f"ND L start {args.start_index}: {status}; calls={count}; "
          f"best L={record['best_L']}; elapsed_s={record['elapsed_s']:.1f}")
    if status == "FAIL":
        raise RuntimeError(message)


if __name__ == "__main__":
    main()
