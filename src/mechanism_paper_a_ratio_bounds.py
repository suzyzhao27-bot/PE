"""Fixed 80 paper A ratio bounds from archived length envelopes and LP."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

from mechanism_paper_a_thickness_bound import length_range


ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "原始资料/赵珈绮/25_论文贡献审计与执行裁剪_v1.0.0"
SOURCE = ROOT / "原始资料/赵珈绮/25_数据约束机制研究_学生任务书_v1.0.1/inputs/20/interval_transport.csv"
OUT = ROOT / "data/processed/salt_transport/25_论文裁剪A五目标界_20261031.csv"
SAMPLES = ("B65", "B700", "B2600", "R65", "R700", "R2600")
TARGETS = (("B2600_over_B65", 0, 2), ("R2600_over_R65", 3, 5),
           ("R65_over_B65", 0, 3), ("R700_over_B700", 1, 4),
           ("R2600_over_B2600", 2, 5))
FIELDS = ("case_id", "window", "target", "family", "scenario", "delta_rho",
          "external", "Gamma", "mobility_sensitivity_fraction",
          "lower_b_ratio", "upper_b_ratio", "lower_status", "upper_status",
          "lower_scaled_primal_residual", "upper_scaled_primal_residual",
          "LP_calls_for_row", "numeric_scope")


def solve_ratio(rho: list[float], delta: float, gamma: int, i: int, j: int,
                maximize: bool) -> tuple[str, float | None, float | None]:
    scale = float(np.median(rho))
    lower = [(1 - delta) * x / scale for x in rho]
    upper = [(1 + delta) * x / scale for x in rho]
    # Charnes-Cooper: V_k=v_k*t, R_g=r_g*t, X_k=xi_k*t, t=1/v_j.
    # Variable order V[0:6], R_B, R_R, X[0:6], t.
    n = 15
    matrix, rhs = [], []
    for k in range(6):
        group = 6 if k < 3 else 7
        row = np.zeros(n)
        row[8 + k], row[group], row[k] = lower[k], -1, -1
        matrix.append(row)
        rhs.append(0.)
        row = np.zeros(n)
        row[k], row[group], row[8 + k] = 1, 1, -upper[k]
        matrix.append(row)
        rhs.append(0.)
        row = np.zeros(n)
        row[8 + k], row[14] = 1, -gamma
        matrix.append(row)
        rhs.append(0.)
        row = np.zeros(n)
        row[14], row[8 + k] = 1, -1
        matrix.append(row)
        rhs.append(0.)
    a_ub = np.asarray(matrix)
    b_ub = np.asarray(rhs)
    a_eq = np.zeros((1, n))
    a_eq[0, j] = 1
    objective = np.zeros(n)
    objective[i] = -1 if maximize else 1
    bounds = [(0, None)] * n
    bounds[14] = (1 / (upper[j] * gamma), None)
    result = linprog(objective, A_ub=a_ub, b_ub=b_ub,
                     A_eq=a_eq, b_eq=[1.], bounds=bounds,
                     method="highs", options={"primal_feasibility_tolerance": 1e-9,
                                              "dual_feasibility_tolerance": 1e-9})
    if result.status == 3 and maximize:
        return "UNBOUNDED", None, None
    if result.status != 0 or result.x is None:
        return f"NUMERIC_UNRESOLVED_{result.status}", None, None
    residual = max(0., float(np.max(a_ub @ result.x - b_ub)),
                   float(np.max(np.abs(a_eq @ result.x - 1))))
    if residual > 1e-8 or result.x[14] <= 0:
        return "NUMERIC_UNRESOLVED_RESIDUAL", None, residual
    return ("FINITE" if result.x[i] > 1e-10 else "ZERO_INFIMUM",
            max(0., float(result.x[i])), residual)


def main() -> None:
    with (PACKAGE / "原矩阵论文用途子集.json").open(encoding="utf-8-sig") as stream:
        cases = json.load(stream)
    with (PACKAGE / "论文贡献可行性核查.json").open(encoding="utf-8-sig") as stream:
        extrema = json.load(stream)["operator_extrema"]
    with SOURCE.open(encoding="utf-8-sig", newline="") as stream:
        source = list(csv.DictReader(stream))
    lookup = {(r["scenario"], r["window_id"], r["sample_id"]): r for r in source}
    if len(cases) != 14 or len(lookup) != len(source):
        raise ValueError("Input identity mismatch")
    output = []
    for case in cases:
        family = case["family"]
        delta = case["delta_rho"]
        gamma = case["Gamma"]
        sensitivity = case["mobility_sensitivity_fraction"]
        lengths = {h: length_range(extrema, family, h, sensitivity)
                   for h in (65, 700, 2600)}
        for window in case["paper_windows"]:
            frozen = [lookup[(case["scenario"], window, sample)] for sample in SAMPLES]
            if any(r["status"] != "VALID" or float(r["lambda_s_inv"]) <= 0 for r in frozen):
                raise ValueError("Invalid frozen rate")
            rho = [1 / float(r["lambda_s_inv"]) for r in frozen]
            for target, i, j in TARGETS:
                low_status, low_value, low_resid = solve_ratio(rho, delta, gamma, i, j, False)
                high_status, high_value, high_resid = solve_ratio(rho, delta, gamma, i, j, True)
                if low_status.startswith("NUMERIC") or high_status.startswith("NUMERIC"):
                    raise ArithmeticError(f"LP unresolved: {case['case_id']}, {window}, {target}")
                h_i, h_j = int(SAMPLES[i][1:]), int(SAMPLES[j][1:])
                low = low_value * lengths[h_j][0] / lengths[h_i][1]
                high = (high_value * lengths[h_j][1] / lengths[h_i][0]
                        if high_value is not None else None)
                output.append({
                    "case_id": case["case_id"], "window": window, "target": target,
                    "family": family, "scenario": case["scenario"],
                    "delta_rho": delta, "external": case["external"],
                    "Gamma": gamma, "mobility_sensitivity_fraction": sensitivity,
                    "lower_b_ratio": low, "upper_b_ratio": high if high is not None else "UNBOUNDED",
                    "lower_status": low_status, "upper_status": high_status,
                    "lower_scaled_primal_residual": low_resid,
                    "upper_scaled_primal_residual": high_resid if high_resid is not None else "NA",
                    "LP_calls_for_row": 2,
                    "numeric_scope": "ARCHIVED_OPERATOR_T12_UNRESOLVED_OUTER_LENGTH_ENVELOPE",
                })
    if len(output) != 80:
        raise ValueError("Expected 80 fixed A target rows")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    print(f"rows={len(output)} LP_calls={2*len(output)} "
          f"unbounded_upper={sum(r['upper_status'] == 'UNBOUNDED' for r in output)}")


if __name__ == "__main__":
    main()
