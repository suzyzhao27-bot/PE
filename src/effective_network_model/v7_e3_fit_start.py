"""One frozen E3 NS refit start on perturbed B raw readings."""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.optimize import least_squares

from v7_d1_objective import H_VALUES, starts
from v7_e2_inverse_start import OUT, init_worker, predict, residuals


class BudgetReached(Exception):
    pass


def data_for_task(eta: float, seed: int):
    with (OUT / "v7_E3_perturbed_B_readings.csv").open(
            encoding="utf-8-sig", newline="") as stream:
        rows = [r for r in csv.DictReader(stream)
                if float(r["eta"]) == eta and int(r["seed"]) == seed]
    assert len(rows) == 42
    data = {}
    for h in H_VALUES:
        group = sorted((r for r in rows if int(r["h_nm"]) == h),
                       key=lambda r: float(r["time_min"]))
        assert len(group) == 7
        times = np.array([float(r["time_min"]) for r in group])
        observed = np.array(
            [float(r["increment_perturbed_uS_cm"]) for r in group])
        assert times[0] == 0 and observed[0] == 0
        data[h] = {"times": times, "observed": observed}
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("eta", type=float, choices=(0.01, 0.05))
    parser.add_argument("seed", type=int, choices=range(8101, 8107))
    parser.add_argument("start_index", type=int, choices=range(4))
    args = parser.parse_args()
    key = (args.eta, args.seed, args.start_index)
    path = OUT / "v7_E3_fit_runs.jsonl"
    if path.exists():
        for line in path.read_text(encoding="utf-8").splitlines():
            old = json.loads(line)
            if (old["eta"], old["seed"], old["start_index"]) == key:
                raise RuntimeError("This E3 frozen start is already recorded")
    data = data_for_task(args.eta, args.seed)
    x0 = starts("NS")[args.start_index]
    lower, upper = np.zeros(3), np.array([5.0, 6.0, 6.0])
    count = 0
    best = {"L": float("inf"), "x": None}
    started = time.perf_counter()
    pool = ProcessPoolExecutor(max_workers=6, initializer=init_worker,
                               initargs=(data,))

    def objective(x):
        nonlocal count
        if count >= 200:
            raise BudgetReached()
        count += 1
        predicted, diagnostics = predict(pool, 10 ** x[0], x[2], 10 ** x[1])
        if (max(v[0] for v in diagnostics.values()) > 1e-8
                or min(v[1] for v in diagnostics.values()) < -1e-12):
            raise RuntimeError("conservation or nonnegative state check failed")
        r = residuals(data, predicted)
        loss = float(r @ r)
        if loss < best["L"]:
            best.update(L=loss, x=x.tolist())
        if count % 10 == 0:
            print(f"E3 {key}: calls {count}/200; best L={best['L']:.6g}",
                  flush=True)
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
        "model": "NS", "task": "E3_raw_reading_perturbation",
        "eta": args.eta, "seed": args.seed, "start_index": args.start_index,
        "start_transformed": x0.tolist(),
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
    with path.open("a", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(record, ensure_ascii=False) + "\n")
    print(f"E3 {key}: {status}; calls={count}; best L={record['best_L']}")
    if status == "FAIL":
        raise RuntimeError(message)


if __name__ == "__main__":
    main()
