"""Prescribed ±5% Dbar sensitivity of the fixed paper D projection."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

from mechanism_paper_d_neutral_projection import neutral_ratio_bounds


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data/processed/salt_transport"
PACKAGE = ROOT / "原始资料/赵珈绮/25_论文贡献审计与执行裁剪_v1.0.0"
RATIO_OUT = BASE / "25_论文裁剪D迁移系数敏感性_速率比_20261031.csv"
CONC_OUT = BASE / "25_论文裁剪D迁移系数敏感性_尿素浓度_20261031.csv"
GAIN_OUT = BASE / "25_论文裁剪D迁移系数敏感性_处理响应_20261031.csv"


def read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def write(path: Path, rows: list[dict], fields: tuple[str, ...]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    with (PACKAGE / "论文贡献可行性核查.json").open(encoding="utf-8-sig") as stream:
        extrema = json.load(stream)["operator_extrema"]
    with (PACKAGE / "原矩阵论文用途子集.json").open(encoding="utf-8-sig") as stream:
        core = [x for x in json.load(stream) if x["paper_neutral_projection"]]
    if len(core) != 2 or {x["scenario"] for x in core} != {"C1_T85", "C2_T85"}:
        raise ValueError("Unexpected D scope")
    if {x["delta_rho"] for x in core} != {.2}:
        raise ValueError("Unexpected D rho tolerance")
    delta = .2
    ratio_rows = []
    ratios = {}
    for h in (65, 700, 2600):
        base = neutral_ratio_bounds(extrema, h)
        wide = neutral_ratio_bounds(extrema, h, .05)
        if not (0 < wide[0] <= base[0] <= base[1] <= wide[1]):
            raise ArithmeticError("Mobility sensitivity is not an envelope")
        ratios[h] = (base, wide)
        ratio_rows.append(dict(h_nm=h, baseline_min=base[0], baseline_max=base[1],
                               mobility_5pct_min=wide[0], mobility_5pct_max=wide[1],
                               scope="ARCHIVED_G1_OPERATOR_SENSITIVITY_NOT_ERROR_LIMIT"))

    source_conc = [r for r in read(BASE / "25_论文裁剪D尿素浓度条件范围_20261031.csv")
                   if r["external_relation"] == "DIFFUSION_SCALED_SHARED"]
    if len(source_conc) != 24:
        raise ValueError("Expected 24 baseline D concentration rows")
    conc_rows = []
    for r in source_conc:
        h = int(r["h_nm"])
        base, wide = ratios[h]
        if not (math.isclose(float(r["lambda_u_over_s_min"]), base[0], rel_tol=1e-12)
                and math.isclose(float(r["lambda_u_over_s_max"]), base[1], rel_tol=1e-12)):
            raise ArithmeticError("Baseline ratio identity mismatch")
        lo = 5 * -math.expm1(-wide[0] * float(r["salt_exposure_min"]))
        hi = 5 * -math.expm1(-wide[1] * float(r["salt_exposure_max"]))
        if not (0 <= lo <= float(r["urea_receiver_min_mM"])
                <= float(r["urea_receiver_max_mM"]) <= hi <= 5):
            raise ArithmeticError("Concentration envelope mismatch")
        conc_rows.append(dict(scenario=r["scenario"], sample_id=r["sample_id"],
                              h_nm=h, time_min=r["time_min"],
                              baseline_min_mM=r["urea_receiver_min_mM"],
                              baseline_max_mM=r["urea_receiver_max_mM"],
                              mobility_5pct_min_mM=lo, mobility_5pct_max_mM=hi,
                              scope="OLD_DEVICE_CONDITIONAL_SENSITIVITY_ONLY"))

    source_gain = [r for r in read(BASE / "25_论文裁剪D尿素处理响应范围_20261031.csv")
                   if r["external_relation"] == "DIFFUSION_SCALED_SHARED"]
    if len(source_gain) != 12:
        raise ValueError("Expected 12 baseline D gain rows")
    gain_rows = []
    for r in source_gain:
        h = int(r["h_nm"])
        base, wide = ratios[h]
        salt_gain = float(r["salt_gain_R_over_B"])
        expected_lo = salt_gain * base[0] / base[1] * (1 - delta) / (1 + delta)
        expected_hi = salt_gain * base[1] / base[0] * (1 + delta) / (1 - delta)
        if not (math.isclose(expected_lo, float(r["urea_gain_min"]), rel_tol=1e-12)
                and math.isclose(expected_hi, float(r["urea_gain_max"]), rel_tol=1e-12)):
            raise ArithmeticError("Baseline gain identity mismatch")
        lo = salt_gain * wide[0] / wide[1] * (1 - delta) / (1 + delta)
        hi = salt_gain * wide[1] / wide[0] * (1 + delta) / (1 - delta)
        if not (0 < lo <= expected_lo <= expected_hi <= hi):
            raise ArithmeticError("Gain envelope mismatch")
        gain_rows.append(dict(scenario=r["scenario"], h_nm=h, window=r["window"],
                              salt_gain_R_over_B=salt_gain,
                              baseline_min=expected_lo, baseline_max=expected_hi,
                              mobility_5pct_min=lo, mobility_5pct_max=hi,
                              excludes_unity_under_5pct=lo > 1,
                              scope="OLD_DEVICE_GROUP_RATE_SENSITIVITY_ONLY"))
    if len({(r["scenario"], r["sample_id"], r["time_min"]) for r in conc_rows}) != 24:
        raise ValueError("Duplicate concentration identity")
    if len({(r["scenario"], r["h_nm"], r["window"]) for r in gain_rows}) != 12:
        raise ValueError("Duplicate gain identity")
    write(RATIO_OUT, ratio_rows, ("h_nm", "baseline_min", "baseline_max",
                                  "mobility_5pct_min", "mobility_5pct_max", "scope"))
    write(CONC_OUT, conc_rows, ("scenario", "sample_id", "h_nm", "time_min",
                                "baseline_min_mM", "baseline_max_mM",
                                "mobility_5pct_min_mM", "mobility_5pct_max_mM", "scope"))
    write(GAIN_OUT, gain_rows, ("scenario", "h_nm", "window", "salt_gain_R_over_B",
                                "baseline_min", "baseline_max", "mobility_5pct_min",
                                "mobility_5pct_max", "excludes_unity_under_5pct", "scope"))
    print(f"ratio={len(ratio_rows)} concentration={len(conc_rows)} gain={len(gain_rows)} "
          f"gain_excludes_unity={sum(r['excludes_unity_under_5pct'] for r in gain_rows)}")


if __name__ == "__main__":
    main()
