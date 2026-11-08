"""Official 12-row paper-scope B rate-range necessity check, no PB or ODE."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "原始资料/赵珈绮/25_数据约束机制研究_学生任务书_v1.0.1/inputs/20/interval_transport.csv"
PILOT = ROOT / "原始资料/赵珈绮/25_论文贡献审计与执行裁剪_v1.0.0/论文贡献可行性核查.json"
SCOPE = ROOT / "原始资料/赵珈绮/25_论文贡献审计与执行裁剪_v1.0.0/原矩阵论文用途子集.json"
OUT = ROOT / "data/processed/salt_transport/25_论文裁剪B速率范围诊断_20261031.csv"
FIELDS = ("case_id", "sample_id", "h_nm", "scenario", "delta_rho",
          "max_lambda_lower_s_inv", "max_lower_window",
          "min_lambda_upper_s_inv", "min_upper_window",
          "required_rate_ratio", "G1_Lmin_nm", "G1_Lmax_nm",
          "G1_allowed_rate_ratio", "decision", "numeric_scope")


def main() -> None:
    with PILOT.open(encoding="utf-8-sig") as stream:
        pilot = json.load(stream)
    with SCOPE.open(encoding="utf-8-sig") as stream:
        selected = [x for x in json.load(stream) if x["paper_rate_range_diagnostic"]]
    if len(selected) != 2 or {x["scenario"] for x in selected} != {"C1_T85", "C2_T85"}:
        raise ValueError("Unexpected paper B scope")
    with SOURCE.open(encoding="utf-8-sig", newline="") as stream:
        source = list(csv.DictReader(stream))
    lookup = {(r["scenario"], r["window_id"], r["sample_id"]): r for r in source}
    if len(lookup) != len(source):
        raise ValueError("Duplicate source identity")

    # Rebuild length ranges from archived operator extrema and fixed geometry.
    limits = {}
    for h in (65, 700, 2600):
        candidates = []
        for e in pilot["operator_extrema"]:
            a = e["a_nm"]
            if a not in (10, 25, 50):
                continue
            for tau in (1, 2, 4):
                if e["sigma_C_m2"] != 0 and tau * h / a < 5:
                    continue
                entry = math.pi * a / (2 * tau)
                candidates.append((h / e["Dbar_max_over_Ds"] + entry,
                                   h / e["Dbar_min_over_Ds"] + entry))
        low = min(x[0] for x in candidates)
        high = max(x[1] for x in candidates)
        saved = pilot["length_bounds"]["G1"][str(h)]
        if not (math.isclose(low, saved["min_nm"], rel_tol=1e-12) and
                math.isclose(high, saved["max_nm"], rel_tol=1e-12)):
            raise ValueError(f"G1 length bound mismatch: {h}")
        limits[h] = (low, high)

    output = []
    for case in selected:
        delta = case["delta_rho"]
        if case["family"] != "G1" or delta != .2 or case["external"] != "E1":
            raise ValueError("Unexpected B condition")
        for sample in ("B65", "B700", "B2600", "R65", "R700", "R2600"):
            h = int(sample[1:])
            rates = []
            for window in ("I1", "I2", "I3", "I4", "I5"):
                row = lookup[(case["scenario"], window, sample)]
                rate = float(row["lambda_s_inv"])
                if row["status"] != "VALID" or rate <= 0:
                    raise ValueError(f"Invalid rate: {case['scenario']}, {sample}, {window}")
                rates.append((rate, window))
            fastest = max(rates)
            slowest = min(rates)
            max_lower = fastest[0] / (1 + delta)
            min_upper = slowest[0] / (1 - delta)
            required = max_lower / min_upper
            lmin, lmax = limits[h]
            allowed = lmax / lmin
            output.append({
                "case_id": case["case_id"], "sample_id": sample, "h_nm": h,
                "scenario": case["scenario"], "delta_rho": delta,
                "max_lambda_lower_s_inv": max_lower, "max_lower_window": fastest[1],
                "min_lambda_upper_s_inv": min_upper, "min_upper_window": slowest[1],
                "required_rate_ratio": required, "G1_Lmin_nm": lmin,
                "G1_Lmax_nm": lmax, "G1_allowed_rate_ratio": allowed,
                "decision": ("INCOMPATIBLE_NECESSARY_RATE_RANGE" if required > allowed * (1 + 1e-12)
                             else "NOT_EXCLUDED_BY_RATE_RANGE"),
                "numeric_scope": "ARCHIVED_PB_OPERATOR_T12_UNRESOLVED",
            })
    if len(output) != 12:
        raise ValueError("Expected twelve official B rows")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    print(f"rows={len(output)} incompatible="
          f"{sum(r['decision'] == 'INCOMPATIBLE_NECESSARY_RATE_RANGE' for r in output)}")


if __name__ == "__main__":
    main()
