"""One frozen E2 inverse start on actual sparse or paired dense synthetic data."""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.optimize import least_squares

from v7_d1_objective import H_VALUES, ROOT, graphs
from v7_network_transport import SPEC_PATH, closed_baths, system


OUT = ROOT / "WP5"
_DATA = None
_GRAPHS = None
_SPEC = None


class BudgetReached(Exception):
    pass


def init_worker(data):
    global _DATA, _GRAPHS, _SPEC
    _DATA = data
    _GRAPHS = {(h, g["seed"]): g for h, group in graphs().items()
               for g in group}
    _SPEC = json.loads(SPEC_PATH.read_text(encoding="utf-8-sig"))


def forward_one(task):
    h, seed, alpha, beta = task
    times = _DATA[h]["times"] * 60.0
    capacities, edges, _ = system(_GRAPHS[h, seed], beta, alpha, _SPEC)
    run = closed_baths(capacities, edges, times)
    return h, seed, run["u_R"], run["mass_relative_error"], run["minimum_u"]


def predict(pool, alpha, beta, q):
    group = {h: [] for h in H_VALUES}
    tasks = [(h, seed, alpha, beta) for h in H_VALUES
             for seed in range(7001, 7009)]
    for h, seed, right, mass, minimum in pool.map(forward_one, tasks):
        group[h].append((seed, right, mass, minimum))
    predictions, diagnostics = {}, {}
    for h in H_VALUES:
        rows = sorted(group[h], key=lambda r: r[0])
        assert [r[0] for r in rows] == list(range(7001, 7009))
        predictions[h] = q * np.mean(np.vstack([r[1] for r in rows]), axis=0)
        diagnostics[h] = (max(r[2] for r in rows),
                          min(r[3] for r in rows))
    return predictions, diagnostics


def data_for_task(truth_label, noise_eta, sampling):
    with (OUT / "v7_E2_synthetic_observations.csv").open(
            encoding="utf-8-sig", newline="") as stream:
        rows = [r for r in csv.DictReader(stream)
                if (r["truth_label"], float(r["noise_eta"]), r["sampling"]) ==
                (truth_label, noise_eta, sampling)]
    expected = 42 if sampling == "sparse_actual" else 726
    assert len(rows) == expected
    data = {}
    for h in H_VALUES:
        group = sorted((r for r in rows if int(r["h_nm"]) == h),
                       key=lambda r: float(r["time_min"]))
        times = np.array([float(r["time_min"]) for r in group])
        observed = np.array([float(r["increment_uS_cm"]) for r in group])
        assert times[0] == 0 and observed[0] == 0
        assert len(group) == (7 if sampling == "sparse_actual" else 121)
        data[h] = {"times": times, "observed": observed}
    return data


def residuals(data, predicted):
    result = []
    for h in H_VALUES:
        obs = data[h]["observed"][1:]
        pred = predicted[h][1:]
        result.extend((np.arcsinh(pred) - np.arcsinh(obs)) /
                      math.sqrt(len(obs) * len(H_VALUES)))
    result = np.asarray(result)
    assert np.all(np.isfinite(result))
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("truth_label", choices=("fullB_best", "equivalent_1"))
    parser.add_argument("noise_eta", type=float, choices=(0.0, 0.01))
    parser.add_argument("sampling", choices=("sparse_actual", "dense_1min"))
    parser.add_argument("startset_index", type=int, choices=(0, 1))
    parser.add_argument("start_index", type=int, choices=range(4))
    args = parser.parse_args()
    key = (args.truth_label, args.noise_eta, args.sampling,
           args.startset_index, args.start_index)
    path = OUT / "v7_E2_inverse_runs.jsonl"
    if path.exists():
        for line in path.read_text(encoding="utf-8").splitlines():
            old = json.loads(line)
            if (old["truth_label"], old["noise_eta"], old["sampling"],
                    old["startset_index"], old["start_index"]) == key:
                raise RuntimeError("This E2 frozen inverse start is already recorded")
    with (OUT / "v7_E2_start_sets.csv").open(
            encoding="utf-8-sig", newline="") as stream:
        starts = [r for r in csv.DictReader(stream)
                  if (int(r["startset_index"]), int(r["start_index"])) ==
                  (args.startset_index, args.start_index)]
    assert len(starts) == 1
    x0 = np.array([float(starts[0][k])
                   for k in ("log10_alpha", "log10_Q", "beta")])
    data = data_for_task(args.truth_label, args.noise_eta, args.sampling)
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
            print(f"E2 {key}: calls {count}/200; best L={best['L']:.6g}",
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
        "model": "NS", "task": "E2_synthetic_inverse",
        "truth_label": args.truth_label, "noise_eta": args.noise_eta,
        "sampling": args.sampling,
        "startset_index": args.startset_index, "start_index": args.start_index,
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
    print(f"E2 {key}: {status}; calls={count}; best L={record['best_L']}")
    if status == "FAIL":
        raise RuntimeError(message)


if __name__ == "__main__":
    main()
