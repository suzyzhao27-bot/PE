"""Refit the six frozen endpoints at preselected Robin boundary values.

This is a sensitivity check, not a fit of r or evidence for a real PE boundary.
The reduced resistance uses the existing core field and independent adaptive
quadrature. No raw source or result file is changed.
"""

import csv
import math
import warnings
from pathlib import Path

import numpy as np
from scipy.integrate import IntegrationWarning, quad
from scipy.optimize import differential_evolution, minimize

from core.theory_model import field, log_resistance_nm


SOURCE = Path(__file__).resolve().parent / "core" / "inputs" / "Book1__Sheet1.csv"
THICKNESS = (65.0, 120.0, 200.0, 350.0, 700.0, 2600.0)
COLUMNS = ("col_2_B", "col_4_D", "col_5_E", "col_6_F", "col_7_G", "col_8_H")
BOUNDARIES = (-0.5, 0.0, 0.25, 0.5, 0.75, 1.0)
BOUNDS = ((math.log(0.01), math.log(100.0)),
          (math.log(0.1), math.log(10000.0)))
ACCESS_NM = math.pi * 25.0 / 2.0


def observations():
    with SOURCE.open(encoding="utf-8-sig", newline="") as source:
        rows = list(csv.DictReader(source))
    if len(rows) != 7 or any(float(rows[index][key]) != time
            for index, time in ((0, 0), (6, 120))
            for key in ("col_1_A", "col_3_C")):
        raise ValueError("Expected the frozen six-curve 0/120 min input.")
    increments = np.array([float(rows[-1][key]) - float(rows[0][key])
                           for key in COLUMNS])
    if not np.all(np.isfinite(increments)) or np.any(increments <= 0):
        raise ValueError("Endpoint increments must be finite and positive.")
    return np.log(increments)


def log_resistance(h, ws, ell, boundary, require_adaptive=False):
    shift = field(0.0, h, ws, ell, boundary)
    with warnings.catch_warnings():
        warnings.simplefilter("error", IntegrationWarning)
        try:
            scaled_integral = quad(
                lambda z: math.exp(field(z, h, ws, ell, boundary) - shift),
                0.0, h, epsabs=1e-11, epsrel=1e-10, limit=200,
                points=[h / 2.0],
            )[0]
        except IntegrationWarning:
            if require_adaptive:
                raise RuntimeError("Adaptive quadrature failed at the selected fit.")
            # Use the existing core quadrature only for difficult trial points.
            return log_resistance_nm(h, ws, ell, 25.0, boundary, n=8192)
    internal = shift + math.log(scaled_integral)
    high = max(internal, math.log(ACCESS_NM))
    return high + math.log(math.exp(internal - high)
                           + math.exp(math.log(ACCESS_NM) - high))


def fit_at_boundary(boundary, observed_log, train_indices=None):
    train = np.arange(len(THICKNESS)) if train_indices is None else np.array(train_indices)

    def evaluate(log_parameters):
        ws, ell = np.exp(log_parameters)
        log_resistances = np.array([
            log_resistance(h, ws, ell, boundary) for h in THICKNESS
        ])
        log_amplitude = float(np.mean((observed_log + log_resistances)[train]))
        residual = log_amplitude - log_resistances - observed_log
        return float(np.dot(residual[train], residual[train]))

    candidates = [
        (evaluate(np.log([ws, ell])), np.log([ws, ell]))
        for ws in (1.0, 3.0, 7.0, 20.0, 80.0)
        for ell in (30.0, 100.0, 300.0, 900.0, 3000.0, 9000.0)
    ]
    starts = [item[1] for item in sorted(candidates, key=lambda item: item[0])[:8]]
    best = None
    for start in starts:
        result = minimize(evaluate, start, method="L-BFGS-B", bounds=BOUNDS,
                          options={"ftol": 1e-13, "gtol": 1e-9, "maxiter": 400})
        if best is None or result.fun < best.fun:
            best = result
    if not best.success:
        polished = minimize(evaluate, best.x, method="Powell", bounds=BOUNDS,
                            options={"xtol": 1e-8, "ftol": 1e-12, "maxiter": 1000})
        if polished.success and polished.fun <= best.fun + 1e-8:
            best = polished
        else:
            global_result = differential_evolution(
                evaluate, BOUNDS, seed=0, tol=1e-9, maxiter=500, polish=True
            )
            if global_result.success and global_result.fun <= best.fun + 1e-8:
                best = global_result
            else:
                raise RuntimeError(
                    f"Optimizer failed at r={boundary}, train={train}: {best.message}"
                )
    ws, ell = np.exp(best.x)
    log_resistances = np.array([log_resistance(h, ws, ell, boundary,
                                              require_adaptive=True)
                                for h in THICKNESS])
    reference = np.array([log_resistance_nm(h, ws, ell, 25.0, boundary)
                          for h in THICKNESS])
    quadrature_difference = float(np.max(np.abs(log_resistances - reference)))
    if quadrature_difference > 1e-5:
        raise RuntimeError(f"Quadratures disagree at r={boundary}: {quadrature_difference}")
    result = dict(r=boundary, ws=float(ws), ell_nm=float(ell),
                train_log_rmse=math.sqrt(float(best.fun) / len(train)),
                ell_at_upper_bound=ell >= 9999.99,
                ws_at_upper_bound=ws >= 99.999,
                quadrature_max_log_difference=quadrature_difference)
    if len(train) == len(THICKNESS) - 1:
        held = next(index for index in range(len(THICKNESS)) if index not in train)
        log_amplitude = float(np.mean((observed_log + log_resistances)[train]))
        result["held_h_nm"] = THICKNESS[held]
        result["held_prediction_over_observed"] = math.exp(
            log_amplitude - log_resistances[held] - observed_log[held]
        )
    return result


def main():
    observed_log = observations()
    results = [fit_at_boundary(boundary, observed_log) for boundary in BOUNDARIES]
    baseline = results[1]
    if abs(baseline["ws"] - 7.0397905) > 0.01 or abs(baseline["ell_nm"] - 893.9539) > 1.0:
        raise RuntimeError("r=0 refit disagrees with the delivered v5 reference.")
    for result in results:
        print(" ".join(f"{key}={value}" for key, value in result.items()), flush=True)
    for boundary in (0.0, 0.75):
        held = [fit_at_boundary(boundary, observed_log,
                                [index for index in range(len(THICKNESS)) if index != omitted])
                for omitted in range(len(THICKNESS))]
        held_log_errors = [math.log(row["held_prediction_over_observed"])
                           for row in held]
        print(f"leave_one_out_r={boundary} held_log_rmse="
              f"{math.sqrt(sum(value * value for value in held_log_errors) / len(held))}",
              flush=True)
        for row in held:
            print(f"held_r={boundary} h={row['held_h_nm']} "
                  f"prediction_over_observed={row['held_prediction_over_observed']} "
                  f"ws={row['ws']} ell_nm={row['ell_nm']}", flush=True)
    nested = []
    for omitted in range(len(THICKNESS)):
        train = [index for index in range(len(THICKNESS)) if index != omitted]
        options = [fit_at_boundary(boundary, observed_log, train)
                   for boundary in BOUNDARIES]
        selected = min(options, key=lambda row: row["train_log_rmse"])
        nested.append(selected)
        print(f"nested_h={THICKNESS[omitted]} selected_r={selected['r']} "
              f"train_log_rmse={selected['train_log_rmse']} "
              f"prediction_over_observed={selected['held_prediction_over_observed']}",
              flush=True)
    errors = [math.log(row["held_prediction_over_observed"]) for row in nested]
    print(f"nested_held_log_rmse={math.sqrt(sum(error * error for error in errors) / len(errors))}",
          flush=True)


if __name__ == "__main__":
    main()
