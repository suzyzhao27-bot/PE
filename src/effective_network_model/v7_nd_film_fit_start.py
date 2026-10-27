"""One frozen ND whole-membrane B holdout fit: five films train."""

from __future__ import annotations

import argparse
import json
import math
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.optimize import least_squares

from v7_d1_objective import H_VALUES, OUT, curves
from v7_nd_objective import BOUNDS, init_worker, predict_nd, starts_nd


class BudgetReached(Exception):
    pass


def film_residuals(curves_by_h, predictions, hold_h):
    assert hold_h in H_VALUES
    residual = []
    for h in H_VALUES:
        if h == hold_h:
            continue
        points = curves_by_h[h]
        baseline = points[0][1]
        for i, (_, k, _) in enumerate(points[1:], start=1):
            residual.append((math.asinh(predictions[h][i]) - math.asinh(k - baseline)) /
                            math.sqrt(6 * (len(H_VALUES) - 1)))
    result = np.array(residual)
    assert len(result) == 30 and np.all(np.isfinite(result))
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("hold_h_nm", type=int, choices=H_VALUES)
    parser.add_argument("start_index", type=int)
    args = parser.parse_args()
    if not 0 <= args.start_index < 4:
        raise ValueError("whole-film start_index must be 0..3")
    decision = json.loads((OUT / "v7_D2_decision.json").read_text(encoding="utf-8"))
    if not decision["ND_enabled_by_frozen_D2"]:
        raise RuntimeError("Frozen D2 did not activate ND")
    out = OUT / "v7_d1_fit_ND_film_runs.jsonl"
    if out.exists():
        for line in out.read_text(encoding="utf-8").splitlines():
            saved = json.loads(line)
            if (saved["model"], saved["hold_h_nm"], saved["start_index"]) == (
                    "ND", args.hold_h_nm, args.start_index):
                raise RuntimeError("This frozen film-holdout start already has a recorded result")
    data = curves()
    x0 = starts_nd()[args.start_index]
    lower = np.array([lo for lo, _ in BOUNDS])
    upper = np.array([hi for _, hi in BOUNDS])
    count = 0
    best = {"L": float("inf"), "x": None}
    started = time.perf_counter()
    pool = ProcessPoolExecutor(max_workers=6, initializer=init_worker)

    def objective(x):
        nonlocal count
        if count >= 200:
            raise BudgetReached()
        count += 1
        alpha, q, beta, tau0, b = 10 ** x[0], 10 ** x[1], x[2], 10 ** x[3], x[4]
        predicted, diagnostics = predict_nd(pool, alpha, q, beta, tau0, b)
        if max(d["max_mass_relative_error"] for d in diagnostics.values()) > 1e-8 or \
           min(d["minimum_u"] for d in diagnostics.values()) < -1e-12:
            raise RuntimeError("conservation or nonnegative state check failed")
        residual = film_residuals(data, predicted, args.hold_h_nm)
        loss = float(residual @ residual)
        if loss < best["L"]:
            best.update(L=loss, x=x.tolist())
        if count % 10 == 0:
            print(f"ND hold {args.hold_h_nm} start {args.start_index}: "
                  f"calls {count}/200; best L={best['L']:.6g}", flush=True)
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
    record = {"model": "ND", "task": "B_whole_film_holdout",
              "hold_h_nm": args.hold_h_nm, "start_index": args.start_index,
              "start_transformed": x0.tolist(), "bounds_transformed": {
                  "lower": lower.tolist(), "upper": upper.tolist()},
              "settings": {"method": "trf", "jac": "2-point", "x_scale": 1.0,
                           "ftol": 1e-8, "xtol": 1e-8, "gtol": 1e-8,
                           "actual_objective_call_limit": 200, "parallel_workers": 6},
              "status": status, "message": message, "actual_objective_calls": count,
              "best_transformed": best["x"],
              "best_L": best["L"] if np.isfinite(best["L"]) else None,
              "elapsed_s": time.perf_counter() - started, "optimizer": optimizer}
    with out.open("a", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(record, ensure_ascii=False) + "\n")
    print(f"ND hold {args.hold_h_nm} start {args.start_index}: {status}; "
          f"calls={count}; best L={record['best_L']}; elapsed_s={record['elapsed_s']:.1f}")
    if status == "FAIL":
        raise RuntimeError(message)


if __name__ == "__main__":
    main()
