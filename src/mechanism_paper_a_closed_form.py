"""Independent endpoint formulas for all 80 paper A outer-envelope bounds."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

from mechanism_paper_a_thickness_bound import length_range


ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "原始资料/赵珈绮/25_论文贡献审计与执行裁剪_v1.0.0"
SOURCE = ROOT / "原始资料/赵珈绮/25_数据约束机制研究_学生任务书_v1.0.1/inputs/20/interval_transport.csv"
LP = ROOT / "data/processed/salt_transport/25_论文裁剪A五目标界_20261031.csv"
OUT = ROOT / "data/processed/salt_transport/25_论文裁剪A五目标解析复核_20261031.csv"
GROUPS = {"B": ("B65", "B700", "B2600"), "R": ("R65", "R700", "R2600")}
TARGETS = {"B2600_over_B65": ("B65", "B2600"),
           "R2600_over_R65": ("R65", "R2600"),
           "R65_over_B65": ("B65", "R65"),
           "R700_over_B700": ("B700", "R700"),
           "R2600_over_B2600": ("B2600", "R2600")}
FIELDS = ("case_id", "window", "target", "lower_closed", "upper_closed",
          "lower_branch", "upper_branch", "lower_relative_error_vs_LP",
          "upper_relative_error_vs_LP", "status_agrees", "scope")


def ratio_extrema(a: dict[str, float], b: dict[str, float],
                  numerator: str, denominator: str) -> tuple[float, float, str, str]:
    same_group = numerator[0] == denominator[0]
    if not same_group:
        cap_i = min(b[k] for k in GROUPS[numerator[0]] if k != numerator)
        cap_j = min(b[k] for k in GROUPS[denominator[0]] if k != denominator)
        lower = max(0., a[numerator] - cap_i) / b[denominator]
        upper = (math.inf if a[denominator] <= cap_j
                 else b[numerator] / (a[denominator] - cap_j))
        return (lower, upper, "CROSS_GROUP_NUMERATOR_CAP",
                "CROSS_GROUP_DENOMINATOR_ZERO" if math.isinf(upper)
                else "CROSS_GROUP_DENOMINATOR_CAP")

    other = next(k for k in GROUPS[numerator[0]] if k not in (numerator, denominator))
    ai, bj = a[numerator], b[denominator]
    if ai >= bj:
        lower, lower_branch = ai / bj, "SAME_GROUP_R_ZERO"
    elif a[numerator] <= b[other]:
        lower, lower_branch = 0., "SAME_GROUP_NUMERATOR_ZERO"
    else:
        lower = (ai - b[other]) / (bj - b[other])
        lower_branch = "SAME_GROUP_OTHER_CURVE_CAP"

    bi, aj = b[numerator], a[denominator]
    if bi <= aj:
        upper, upper_branch = bi / aj, "SAME_GROUP_R_ZERO"
    elif aj <= b[other]:
        upper, upper_branch = math.inf, "SAME_GROUP_DENOMINATOR_ZERO"
    else:
        upper = (bi - b[other]) / (aj - b[other])
        upper_branch = "SAME_GROUP_OTHER_CURVE_CAP"
    return lower, upper, lower_branch, upper_branch


def main() -> None:
    with (PACKAGE / "原矩阵论文用途子集.json").open(encoding="utf-8-sig") as stream:
        cases = json.load(stream)
    with (PACKAGE / "论文贡献可行性核查.json").open(encoding="utf-8-sig") as stream:
        extrema = json.load(stream)["operator_extrema"]
    with SOURCE.open(encoding="utf-8-sig", newline="") as stream:
        source = list(csv.DictReader(stream))
    with LP.open(encoding="utf-8-sig", newline="") as stream:
        lp_rows = list(csv.DictReader(stream))
    lookup = {(r["scenario"], r["window_id"], r["sample_id"]): r for r in source}
    lp = {(r["case_id"], r["window"], r["target"]): r for r in lp_rows}
    if len(cases) != 14 or len(lookup) != len(source) or len(lp) != 80:
        raise ValueError("Source or LP identity mismatch")
    output = []
    for case in cases:
        lengths = {h: length_range(extrema, case["family"], h,
                                   case["mobility_sensitivity_fraction"])
                   for h in (65, 700, 2600)}
        for window in case["paper_windows"]:
            a, b = {}, {}
            for group in GROUPS.values():
                for sample in group:
                    frozen = lookup[(case["scenario"], window, sample)]
                    if frozen["status"] != "VALID" or float(frozen["lambda_s_inv"]) <= 0:
                        raise ValueError("Invalid frozen rate")
                    rho = 1 / float(frozen["lambda_s_inv"])
                    a[sample] = (1 - case["delta_rho"]) * rho
                    b[sample] = (1 + case["delta_rho"]) * rho * case["Gamma"]
            for target, (i, j) in TARGETS.items():
                vlo, vhi, branch_lo, branch_hi = ratio_extrema(a, b, i, j)
                li, lj = lengths[int(i[1:])], lengths[int(j[1:])]
                lower = vlo * lj[0] / li[1]
                upper = vhi * lj[1] / li[0]
                saved = lp[(case["case_id"], window, target)]
                saved_lo = float(saved["lower_b_ratio"])
                err_lo = abs(lower - saved_lo) / max(1., abs(lower), abs(saved_lo))
                if saved["upper_status"] == "UNBOUNDED":
                    err_hi = 0. if math.isinf(upper) else math.inf
                else:
                    saved_hi = float(saved["upper_b_ratio"])
                    err_hi = abs(upper - saved_hi) / max(1., abs(upper), abs(saved_hi))
                agree = err_lo <= 1e-8 and err_hi <= 1e-8
                output.append({
                    "case_id": case["case_id"], "window": window, "target": target,
                    "lower_closed": lower,
                    "upper_closed": "UNBOUNDED" if math.isinf(upper) else upper,
                    "lower_branch": branch_lo, "upper_branch": branch_hi,
                    "lower_relative_error_vs_LP": err_lo,
                    "upper_relative_error_vs_LP": err_hi,
                    "status_agrees": agree,
                    "scope": "ARCHIVED_OPERATOR_OUTER_LENGTH_ENVELOPE_ONLY",
                })
    if len(output) != 80 or any(not r["status_agrees"] for r in output):
        bad = sum(not r["status_agrees"] for r in output)
        raise ArithmeticError(f"Closed-form and LP mismatch: {bad}")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    print(f"rows={len(output)} max_relative_error="
          f"{max(max(r['lower_relative_error_vs_LP'], r['upper_relative_error_vs_LP']) for r in output):.3g}")


if __name__ == "__main__":
    main()
