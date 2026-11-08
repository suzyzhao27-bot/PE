"""One permitted continuation of a frozen whole-film start, to 400 total calls."""

from __future__ import annotations

import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.optimize import least_squares

from v7_d1_film_fit_start import film_residuals
from v7_d1_objective import (H_VALUES, OUT, curves, init_forward_worker,
                             predict_parallel)


class BudgetReached(Exception):
    pass


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("model", choices=("N0", "NS"))
    parser.add_argument("hold_h_nm", type=int, choices=H_VALUES)
    parser.add_argument("start_index", type=int)
    args = parser.parse_args()
    raw = [json.loads(line) for line in
           (OUT / "v7_d1_fit_film_runs.jsonl").read_text(encoding="utf-8").splitlines()]
    fold = [r for r in raw if r["model"] == args.model and
            r["hold_h_nm"] == args.hold_h_nm]
    if len(fold) != 4 or sorted(r["start_index"] for r in fold) != list(range(4)):
        raise RuntimeError("All four original frozen starts must finish before extension")
    unfinished = sorted((r for r in fold if
                         (r["status"] == "max_eval" or
                          (r["status"] == "max_eval_or_solver_stop" and
                           r.get("optimizer", {}).get("scipy_status") == 0)) and
                         r["best_L"] is not None and
                         np.isfinite(r["best_L"]) and r["best_transformed"]),
                        key=lambda r: (r["best_L"], r["start_index"]))[:2]
    if args.start_index not in [r["start_index"] for r in unfinished]:
        raise RuntimeError("Start is not among the two lowest-loss unfinished starts")
    original = next(r for r in unfinished if r["start_index"] == args.start_index)
    out = OUT / "v7_d1_fit_film_extensions.jsonl"
    if out.exists():
        for line in out.read_text(encoding="utf-8").splitlines():
            saved = json.loads(line)
            if (saved["model"], saved["hold_h_nm"], saved["start_index"]) == (
                    args.model, args.hold_h_nm, args.start_index):
                raise RuntimeError("This start has already used its one extension")
    initial_calls = original["actual_objective_calls"]
    remaining = 400 - initial_calls
    if not 0 < remaining <= 400:
        raise RuntimeError("Invalid cumulative 400-call budget")
    x0 = np.array(original["best_transformed"], dtype=float)
    lower = np.array(original["bounds_transformed"]["lower"], dtype=float)
    upper = np.array(original["bounds_transformed"]["upper"], dtype=float)
    data = curves()
    count = 0
    best = {"L": float(original["best_L"]), "x": x0.tolist()}
    started = time.perf_counter()
    pool = ProcessPoolExecutor(max_workers=6, initializer=init_forward_worker)

    def objective(x):
        nonlocal count
        if count >= remaining:
            raise BudgetReached()
        count += 1
        alpha, q = 10 ** x[0], 10 ** x[1]
        beta = x[2] if args.model == "NS" else 0.0
        predicted, diagnostics = predict_parallel(pool, alpha, beta, q)
        if max(d["max_mass_relative_error"] for d in diagnostics.values()) > 1e-8 or \
           min(d["minimum_u"] for d in diagnostics.values()) < -1e-12:
            raise RuntimeError("conservation or nonnegative state check failed")
        residual = film_residuals(data, predicted, args.hold_h_nm)
        loss = float(residual @ residual)
        if loss < best["L"]:
            best.update(L=loss, x=x.tolist())
        if count % 10 == 0:
            print(f"{args.model} hold {args.hold_h_nm} start {args.start_index} "
                  f"extension: calls {initial_calls + count}/400; "
                  f"best L={best['L']:.6g}", flush=True)
        return residual

    status, message, optimizer = "", "", {}
    try:
        result = least_squares(objective, x0, bounds=(lower, upper), method="trf",
                               jac="2-point", x_scale=1.0, ftol=1e-8, xtol=1e-8,
                               gtol=1e-8, max_nfev=remaining)
        optimizer = {"scipy_success": bool(result.success),
                     "scipy_status": int(result.status),
                     "scipy_nfev": int(result.nfev), "message": str(result.message)}
        status = "solver_success" if result.success else "max_eval_or_solver_stop"
        message = str(result.message)
    except BudgetReached:
        status = "max_eval"
        message = "400 cumulative actual residual calls reached, including finite differences"
    except Exception as exc:
        status = "FAIL"
        message = repr(exc)
    finally:
        pool.shutdown(wait=True)
    record = {"model": args.model, "task": "B_whole_film_holdout_extension",
              "hold_h_nm": args.hold_h_nm, "start_index": args.start_index,
              "initial_status": original["status"],
              "original_start_transformed": original["start_transformed"],
              "resume_from_best_transformed": x0.tolist(),
              "initial_actual_objective_calls": initial_calls,
              "additional_actual_objective_calls": count,
              "cumulative_actual_objective_calls": initial_calls + count,
              "cumulative_call_limit": 400,
              "settings": {**original["settings"],
                           "extension_actual_call_limit": remaining,
                           "cumulative_actual_call_limit": 400},
              "status": status, "message": message,
              "best_transformed": best["x"], "best_L": best["L"],
              "elapsed_s": time.perf_counter() - started, "optimizer": optimizer}
    with out.open("a", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(record, ensure_ascii=False) + "\n")
    print(f"{args.model} hold {args.hold_h_nm} start {args.start_index} extension: "
          f"{status}; cumulative calls={initial_calls + count}; best L={best['L']}")
    if status == "FAIL":
        raise RuntimeError(message)


if __name__ == "__main__":
    main()
