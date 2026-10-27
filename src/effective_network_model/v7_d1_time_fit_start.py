"""One independent frozen B time-holdout fit start (t <= 60 min)."""

from __future__ import annotations

import argparse
import json
import math
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.optimize import least_squares

from v7_d1_objective import H_VALUES, OUT, curves, init_forward_worker, predict_parallel, starts


class BudgetReached(Exception):
    pass


def time_residuals(curves_by_h, predictions, train: bool):
    residual = []
    for h in H_VALUES:
        points = curves_by_h[h]
        baseline = points[0][1]
        chosen = [(i, k - baseline) for i, (t, k, _) in enumerate(points)
                  if i > 0 and ((t <= 60.0) if train else (t > 60.0))]
        assert len(chosen) == 3
        residual.extend((math.asinh(predictions[h][i]) - math.asinh(y)) /
                        math.sqrt(len(chosen) * len(H_VALUES)) for i, y in chosen)
    result = np.array(residual)
    assert len(result) == 18 and np.all(np.isfinite(result))
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("model", choices=("N0", "NS"))
    parser.add_argument("start_index", type=int)
    args = parser.parse_args()
    if not 0 <= args.start_index < 4:
        raise ValueError("time-holdout start_index must be 0..3")
    out = OUT / "v7_d1_fit_time_runs.jsonl"
    if out.exists():
        for line in out.read_text(encoding="utf-8").splitlines():
            saved = json.loads(line)
            if (saved["model"], saved["start_index"]) == (args.model, args.start_index):
                raise RuntimeError("This frozen time start already has a recorded result")
    data = curves()
    x0 = starts(args.model)[args.start_index]
    lower = np.array([0.0, 0.0] + ([0.0] if args.model == "NS" else []))
    upper = np.array([5.0, 6.0] + ([6.0] if args.model == "NS" else []))
    count = 0
    best = {"L": float("inf"), "x": None}
    started = time.perf_counter()
    pool = ProcessPoolExecutor(max_workers=6, initializer=init_forward_worker)

    def objective(x):
        nonlocal count
        if count >= 200:
            raise BudgetReached()
        count += 1
        alpha, q = 10 ** x[0], 10 ** x[1]
        beta = x[2] if args.model == "NS" else 0.0
        predicted, diagnostics = predict_parallel(pool, alpha, beta, q)
        if max(d["max_mass_relative_error"] for d in diagnostics.values()) > 1e-8 or \
           min(d["minimum_u"] for d in diagnostics.values()) < -1e-12:
            raise RuntimeError("conservation or nonnegative state check failed")
        residual = time_residuals(data, predicted, train=True)
        loss = float(residual @ residual)
        if loss < best["L"]:
            best.update(L=loss, x=x.tolist())
        if count % 10 == 0:
            print(f"{args.model} time start {args.start_index}: calls {count}/200; "
                  f"best L={best['L']:.6g}", flush=True)
        return residual

    status = ""
    message = ""
    optimizer = {}
    try:
        result = least_squares(objective, x0, bounds=(lower, upper), method="trf",
                               jac="2-point", x_scale=1.0, ftol=1e-8, xtol=1e-8,
                               gtol=1e-8, max_nfev=200)
        optimizer = {"scipy_success": bool(result.success), "scipy_status": int(result.status),
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
    record = {"model": args.model, "task": "time_B_train_le_60_min",
              "start_index": args.start_index, "start_transformed": x0.tolist(),
              "bounds_transformed": {"lower": lower.tolist(), "upper": upper.tolist()},
              "settings": {"method": "trf", "jac": "2-point", "x_scale": 1.0,
                           "ftol": 1e-8, "xtol": 1e-8, "gtol": 1e-8,
                           "actual_objective_call_limit": 200, "parallel_workers": 6},
              "status": status, "message": message, "actual_objective_calls": count,
              "best_transformed": best["x"],
              "best_L": best["L"] if np.isfinite(best["L"]) else None,
              "elapsed_s": time.perf_counter() - started, "optimizer": optimizer}
    with out.open("a", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(record, ensure_ascii=False) + "\n")
    print(f"{args.model} time start {args.start_index}: {status}; calls={count}; "
          f"best L={record['best_L']}; elapsed_s={record['elapsed_s']:.1f}")
    if status == "FAIL":
        raise RuntimeError(message)


if __name__ == "__main__":
    main()
