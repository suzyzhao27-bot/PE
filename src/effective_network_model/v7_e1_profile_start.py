"""One frozen NS E1 profile start at one of nine transformed positions."""

from __future__ import annotations

import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.optimize import least_squares

from v7_d1_objective import OUT, ROOT, curves, init_forward_worker, predict_parallel, residuals, starts


class BudgetReached(Exception):
    pass


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("parameter_index", type=int, choices=range(3))
    parser.add_argument("position_index", type=int, choices=range(9))
    parser.add_argument("start_index", type=int, choices=range(3))
    args = parser.parse_args()
    spec_path = ROOT / "原始资料/赵珈绮/14_论文补强学生任务书_v7.0.0/configs/frozen_spec.json"
    spec = json.loads(spec_path.read_text(encoding="utf-8-sig"))
    assert spec["optimization"]["starts_profile_sensitivity"] == 3
    assert spec["optimization"]["max_nfev_per_start"] == 200
    assert spec["matrices"]["E1"]["profile_positions"] == 9
    for model in ("N0", "NS"):
        full = json.loads((OUT / f"v7_{model}_fullB_summary.json").read_text(encoding="utf-8"))
        film = json.loads((OUT / f"v7_{model}_film_summary.json").read_text(encoding="utf-8"))
        assert full["status"] == "multi_start_stable"
        assert film["status"] == "all_folds_multi_start_stable"
    nd_film = json.loads((OUT / "v7_ND_film_summary.json").read_text(encoding="utf-8"))
    assert nd_film["status"] == "some_folds_unstable"
    time_n0 = json.loads((OUT / "v7_N0_time_summary.json").read_text(encoding="utf-8"))["E_time"]
    time_ns = json.loads((OUT / "v7_NS_time_summary.json").read_text(encoding="utf-8"))["E_time"]
    assert time_ns < time_n0
    out = ROOT / "WP5/v7_E1_profile_runs.jsonl"
    if out.exists():
        for line in out.read_text(encoding="utf-8").splitlines():
            old = json.loads(line)
            if (old["parameter_index"], old["position_index"], old["start_index"]) == (
                    args.parameter_index, args.position_index, args.start_index):
                raise RuntimeError("This frozen E1 start already has a result")
    names = ("log10_alpha", "log10_Q", "beta")
    lower = np.array([0.0, 0.0, 0.0])
    upper = np.array([5.0, 6.0, 6.0])
    fixed = lower[args.parameter_index] + (upper[args.parameter_index] - lower[args.parameter_index]) * args.position_index / 8
    free = [i for i in range(3) if i != args.parameter_index]
    x0_full = starts("NS")[args.start_index].copy()
    x0_full[args.parameter_index] = fixed
    x0 = x0_full[free]
    data = curves()
    best = {"L": float("inf"), "x": None}
    count = 0
    started = time.perf_counter()
    pool = ProcessPoolExecutor(max_workers=6, initializer=init_forward_worker)

    def objective(z):
        nonlocal count
        if count >= 200:
            raise BudgetReached()
        count += 1
        x = np.empty(3)
        x[args.parameter_index] = fixed
        x[free] = z
        alpha, q, beta = 10 ** x[0], 10 ** x[1], x[2]
        predicted, diagnostics = predict_parallel(pool, alpha, beta, q)
        if (max(d["max_mass_relative_error"] for d in diagnostics.values()) > 1e-8 or
                min(d["minimum_u"] for d in diagnostics.values()) < -1e-12):
            raise RuntimeError("conservation or nonnegative state check failed")
        r = residuals(data, predicted)
        if not np.all(np.isfinite(r)):
            raise RuntimeError("nonfinite profile residual")
        loss = float(r @ r)
        if loss < best["L"]:
            best.update(L=loss, x=x.tolist())
        if count % 10 == 0:
            print(f"NS E1 {names[args.parameter_index]} position {args.position_index} "
                  f"start {args.start_index}: calls {count}/200; best L={best['L']:.7g}", flush=True)
        return r

    status, message, optimizer = "", "", {}
    try:
        result = least_squares(objective, x0, bounds=(lower[free], upper[free]),
                               method="trf", jac="2-point", x_scale=1.0,
                               ftol=1e-8, xtol=1e-8, gtol=1e-8, max_nfev=200)
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
    record = {
        "model": "NS", "task": "E1_fixed_parameter_profile",
        "parameter_index": args.parameter_index, "parameter": names[args.parameter_index],
        "position_index": args.position_index, "position_count": 9,
        "fixed_value_transformed": fixed, "start_index": args.start_index,
        "start_transformed": x0_full.tolist(),
        "settings": {"method": "trf", "jac": "2-point", "x_scale": 1.0,
                     "ftol": 1e-8, "xtol": 1e-8, "gtol": 1e-8,
                     "actual_objective_call_limit": 200, "parallel_workers": 6,
                     "three_starts": "first_three_frozen_fullB_starts_projected_to_free_coordinates",
                     "positions": "nine_inclusive_equal_positions_in_transformed_bounds"},
        "status": status, "message": message, "actual_objective_calls": count,
        "best_transformed": best["x"],
        "best_L": best["L"] if np.isfinite(best["L"]) else None,
        "elapsed_s": time.perf_counter() - started, "optimizer": optimizer,
    }
    with out.open("a", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(record, ensure_ascii=False) + "\n")
    print(f"NS E1 {record['parameter']} position {args.position_index} start {args.start_index}: "
          f"{status}; calls={count}; best L={record['best_L']}")
    if status == "FAIL":
        raise RuntimeError(message)


if __name__ == "__main__":
    main()
