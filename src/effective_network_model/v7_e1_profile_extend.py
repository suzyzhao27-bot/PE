"""One rule-limited E1 continuation from an unfinished start, to 400 calls total."""

from __future__ import annotations

import argparse
import json
import math
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.optimize import least_squares

from v7_d1_objective import ROOT, curves, init_forward_worker, predict_parallel, residuals


class BudgetReached(Exception):
    pass


def unfinished(record: dict) -> bool:
    return record["status"] == "max_eval" or (record["status"] == "max_eval_or_solver_stop" and
            record.get("optimizer", {}).get("scipy_status") == 0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("parameter_index", type=int, choices=range(3))
    parser.add_argument("position_index", type=int, choices=range(9))
    parser.add_argument("start_index", type=int, choices=range(3))
    args = parser.parse_args()
    path = ROOT / "WP5/v7_E1_profile_runs.jsonl"
    raw = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]
    point = [r for r in raw if (r["parameter_index"], r["position_index"]) ==
             (args.parameter_index, args.position_index)]
    if len(point) != 3 or sorted(r["start_index"] for r in point) != [0, 1, 2]:
        raise RuntimeError("All three original starts must finish first")
    candidates = sorted((r for r in point if unfinished(r) and r["best_transformed"] is not None
                         and r["best_L"] is not None and math.isfinite(r["best_L"])),
                        key=lambda r: (r["best_L"], r["start_index"]))[:2]
    if args.start_index not in [r["start_index"] for r in candidates]:
        raise RuntimeError("Start is not among the two lowest-loss unfinished starts")
    good = [r for r in point if r["status"] == "solver_success" and r["best_L"] is not None]
    best_good = min(good, key=lambda r: r["best_L"]) if good else None
    stable = bool(best_good and any(r is not best_good and
                  abs(r["best_L"] - best_good["best_L"]) / best_good["best_L"] <= .01 for r in good))
    if stable and min(r["best_L"] for r in candidates) >= best_good["best_L"]:
        raise RuntimeError("Selected successful point is stable; extension not needed")
    out = ROOT / "WP5/v7_E1_profile_extensions.jsonl"
    existing = ([json.loads(line) for line in out.read_text(encoding="utf-8").splitlines()]
                if out.exists() else [])
    existing = [r for r in existing if (r["parameter_index"], r["position_index"]) ==
                (args.parameter_index, args.position_index)]
    if any(r["start_index"] == args.start_index for r in existing):
        raise RuntimeError("This start already used its one extension")
    effective = [dict(r) for r in point]
    for extension in existing:
        target = next(r for r in effective if r["start_index"] == extension["start_index"])
        target.update(status=extension["status"], best_L=extension["best_L"])
    effective_good = [r for r in effective if r["status"] == "solver_success" and
                      r["best_L"] is not None]
    effective_best = min(effective_good, key=lambda r: r["best_L"]) if effective_good else None
    effective_stable = bool(effective_best and any(r is not effective_best and
                            abs(r["best_L"] - effective_best["best_L"]) /
                            effective_best["best_L"] <= .01 for r in effective_good))
    if effective_stable and not any(unfinished(r) and r["best_L"] is not None and
                                    r["best_L"] < effective_best["best_L"] for r in effective):
        raise RuntimeError("Point is already stable after earlier extension")
    original = next(r for r in candidates if r["start_index"] == args.start_index)
    initial_calls = original["actual_objective_calls"]
    remaining = 400 - initial_calls
    if not 0 < remaining <= 400:
        raise RuntimeError("Invalid cumulative 400-call budget")
    fixed = original["fixed_value_transformed"]
    lower, upper = np.array([0., 0., 0.]), np.array([5., 6., 6.])
    free = [i for i in range(3) if i != args.parameter_index]
    x0_full = np.array(original["best_transformed"], dtype=float)
    if not math.isclose(x0_full[args.parameter_index], fixed, abs_tol=1e-12):
        raise RuntimeError("Original best does not satisfy fixed coordinate")
    x0 = x0_full[free]
    data = curves()
    best = {"L": float(original["best_L"]), "x": x0_full.tolist()}
    count = 0
    started = time.perf_counter()
    pool = ProcessPoolExecutor(max_workers=6, initializer=init_forward_worker)

    def objective(z):
        nonlocal count
        if count >= remaining:
            raise BudgetReached()
        count += 1
        x = np.empty(3)
        x[args.parameter_index] = fixed
        x[free] = z
        predicted, diagnostics = predict_parallel(pool, 10 ** x[0], x[2], 10 ** x[1])
        if (max(d["max_mass_relative_error"] for d in diagnostics.values()) > 1e-8 or
                min(d["minimum_u"] for d in diagnostics.values()) < -1e-12):
            raise RuntimeError("conservation or nonnegative state check failed")
        residual = residuals(data, predicted)
        if not np.all(np.isfinite(residual)):
            raise RuntimeError("nonfinite profile residual")
        loss = float(residual @ residual)
        if loss < best["L"]:
            best.update(L=loss, x=x.tolist())
        if count % 20 == 0:
            print(f"NS E1 extension parameter {args.parameter_index} position {args.position_index} "
                  f"start {args.start_index}: cumulative calls {initial_calls + count}/400; "
                  f"best L={best['L']:.7g}", flush=True)
        return residual

    status, message, optimizer = "", "", {}
    try:
        result = least_squares(objective, x0, bounds=(lower[free], upper[free]),
                               method="trf", jac="2-point", x_scale=1.0,
                               ftol=1e-8, xtol=1e-8, gtol=1e-8, max_nfev=remaining)
        optimizer = {"scipy_success": bool(result.success), "scipy_status": int(result.status),
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
    record = {
        "model": "NS", "task": "E1_fixed_parameter_profile_extension",
        "parameter_index": args.parameter_index, "position_index": args.position_index,
        "start_index": args.start_index, "fixed_value_transformed": fixed,
        "initial_status": original["status"], "original_start_transformed": original["start_transformed"],
        "resume_from_best_transformed": x0_full.tolist(),
        "initial_actual_objective_calls": initial_calls,
        "additional_actual_objective_calls": count,
        "cumulative_actual_objective_calls": initial_calls + count,
        "cumulative_call_limit": 400,
        "settings": {**original["settings"], "extension_actual_call_limit": remaining,
                     "cumulative_actual_call_limit": 400},
        "status": status, "message": message,
        "best_transformed": best["x"], "best_L": best["L"],
        "elapsed_s": time.perf_counter() - started, "optimizer": optimizer,
    }
    with out.open("a", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(record, ensure_ascii=False) + "\n")
    print(f"NS E1 extension parameter {args.parameter_index} position {args.position_index} "
          f"start {args.start_index}: {status}; cumulative calls={initial_calls + count}; L={best['L']}")
    if status == "FAIL":
        raise RuntimeError(message)


if __name__ == "__main__":
    main()
