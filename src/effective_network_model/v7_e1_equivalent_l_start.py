"""One frozen four-start L early calibration for an E1 B-equivalent beta."""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.optimize import least_squares

from v7_d1_objective import ROOT, starts
from v7_d1_l_fit_start import (BudgetReached, curves_l, init_worker,
                               predict_early, residuals_early)


OUT = ROOT / "WP5"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("equivalent_rank", type=int, choices=(1, 2, 3))
    parser.add_argument("start_index", type=int, choices=range(4))
    args = parser.parse_args()
    with (OUT / "v7_E1_equivalent_selection.csv").open(
            encoding="utf-8-sig", newline="") as stream:
        selection = list(csv.DictReader(stream))
    assert len(selection) == 4
    label = f"equivalent_{args.equivalent_rank}"
    chosen = next(r for r in selection if r["label"] == label)
    beta = float(chosen["beta"])
    record_path = OUT / "v7_E1_equivalent_L_fit_runs.jsonl"
    if record_path.exists():
        for line in record_path.read_text(encoding="utf-8").splitlines():
            r = json.loads(line)
            if (r["label"], r["start_index"]) == (label, args.start_index):
                raise RuntimeError("E1 L early fit start already recorded")
    x0 = starts("NS")[args.start_index, :2].copy()
    lower = np.array([0.0, 0.0])
    upper = np.array([5.0, 6.0])
    data = curves_l()
    count = 0
    best = {"L": float("inf"), "x": None}
    started = time.perf_counter()
    pool = ProcessPoolExecutor(max_workers=6, initializer=init_worker,
                               initargs=(beta,))

    def objective(x):
        nonlocal count
        if count >= 200:
            raise BudgetReached()
        count += 1
        predicted, diagnostics = predict_early(pool, 10 ** x[0], 10 ** x[1])
        if (max(v[0] for v in diagnostics.values()) > 1e-8
                or min(v[1] for v in diagnostics.values()) < -1e-12):
            raise RuntimeError("conservation or nonnegative state check failed")
        r = residuals_early(data, predicted)
        loss = float(r @ r)
        if loss < best["L"]:
            best.update(L=loss, x=x.tolist())
        if count % 10 == 0:
            print(f"E1 {label} L start {args.start_index}: "
                  f"calls {count}/200; best L={best['L']:.6g}", flush=True)
        return r

    status, message, optimizer = "", "", {}
    try:
        result = least_squares(objective, x0, bounds=(lower, upper),
                               method="trf", jac="2-point", x_scale=1.0,
                               ftol=1e-8, xtol=1e-8, gtol=1e-8, max_nfev=200)
        optimizer = {"scipy_success": bool(result.success),
                     "scipy_status": int(result.status),
                     "scipy_nfev": int(result.nfev),
                     "message": str(result.message)}
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
    record = {
        "model": "NS", "task": "E1_equivalent_L_early_le_240_min",
        "label": label, "start_index": args.start_index,
        "start_transformed": x0.tolist(),
        "fixed_beta_from_selected_B": beta,
        "bounds_transformed": {"lower": lower.tolist(), "upper": upper.tolist()},
        "settings": {"method": "trf", "jac": "2-point", "x_scale": 1.0,
                     "ftol": 1e-8, "xtol": 1e-8, "gtol": 1e-8,
                     "actual_objective_call_limit": 200,
                     "parallel_workers": 6},
        "status": status, "message": message,
        "actual_objective_calls": count,
        "best_transformed": best["x"],
        "best_L": best["L"] if math.isfinite(best["L"]) else None,
        "elapsed_s": time.perf_counter() - started,
        "optimizer": optimizer,
    }
    with record_path.open("a", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(record, ensure_ascii=False) + "\n")
    print(f"E1 {label} L start {args.start_index}: {status}; "
          f"calls={count}; best L={record['best_L']}")
    if status == "FAIL":
        raise RuntimeError(message)


if __name__ == "__main__":
    main()
