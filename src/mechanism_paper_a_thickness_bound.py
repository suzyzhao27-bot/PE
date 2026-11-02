"""Conservative B2600/B65 bound for the 16 fixed paper A case-windows."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "原始资料/赵珈绮/25_论文贡献审计与执行裁剪_v1.0.0"
SOURCE = ROOT / "原始资料/赵珈绮/25_数据约束机制研究_学生任务书_v1.0.1/inputs/20/interval_transport.csv"
OUT = ROOT / "data/processed/salt_transport/25_论文裁剪A厚薄保守界_20261031.csv"
FIELDS = ("case_id", "window", "family", "scenario", "delta_rho", "external",
          "Gamma", "mobility_sensitivity_fraction", "lambda_B65_s_inv",
          "lambda_B2600_s_inv", "L_B65_max_nm", "L_B2600_min_nm",
          "lower_b_B2600_over_B65", "bound_kind", "numeric_scope")


def length_range(extrema: list[dict], family: str, h: int, sensitivity: float) -> tuple[float, float]:
    if family == "G0":
        entries = [e for e in extrema if e["a_nm"] == 25 and e["sigma_C_m2"] == 0]
        taus = (2,)
    else:
        radii = (10, 25, 50) if family == "G1" else (2, 10, 25, 50)
        entries = [e for e in extrema if e["a_nm"] in radii]
        taus = (1, 2, 4)
    candidates = []
    for e in entries:
        for tau in taus:
            if family == "G1" and e["sigma_C_m2"] != 0 and tau * h / e["a_nm"] < 5:
                continue
            entry = math.pi * e["a_nm"] / (2 * tau)
            d_low = e["Dbar_min_over_Ds"] * (1 - sensitivity)
            d_high = e["Dbar_max_over_Ds"] * (1 + sensitivity)
            if d_low <= 0:
                raise ValueError("Nonpositive mobility bound")
            candidates.append((h / d_high + entry, h / d_low + entry))
    return min(x[0] for x in candidates), max(x[1] for x in candidates)


def main() -> None:
    with (PACKAGE / "原矩阵论文用途子集.json").open(encoding="utf-8-sig") as stream:
        cases = json.load(stream)
    with (PACKAGE / "论文贡献可行性核查.json").open(encoding="utf-8-sig") as stream:
        pilot = json.load(stream)
    with SOURCE.open(encoding="utf-8-sig", newline="") as stream:
        source = list(csv.DictReader(stream))
    lookup = {(r["scenario"], r["window_id"], r["sample_id"]): r for r in source}
    if len(lookup) != len(source) or len(cases) != 14:
        raise ValueError("Source or scope identity mismatch")

    output = []
    for case in cases:
        if not case["spatial_bounds"] or case["external"] not in ("E1", "E2"):
            raise ValueError("Unexpected paper A case")
        family = case["family"]
        sensitivity = case["mobility_sensitivity_fraction"]
        low_thick, _ = length_range(pilot["operator_extrema"], family, 2600, sensitivity)
        _, high_thin = length_range(pilot["operator_extrema"], family, 65, sensitivity)
        if sensitivity == 0:
            for h, calculated in ((65, length_range(pilot["operator_extrema"], family, 65, 0)),
                                  (2600, length_range(pilot["operator_extrema"], family, 2600, 0))):
                saved = pilot["length_bounds"][family][str(h)]
                if not (math.isclose(calculated[0], saved["min_nm"], rel_tol=1e-12) and
                        math.isclose(calculated[1], saved["max_nm"], rel_tol=1e-12)):
                    raise ValueError(f"Length mismatch: {family}, {h}")
        delta = case["delta_rho"]
        gamma = case["Gamma"]
        for window in case["paper_windows"]:
            thin = lookup[(case["scenario"], window, "B65")]
            thick = lookup[(case["scenario"], window, "B2600")]
            if thin["status"] != "VALID" or thick["status"] != "VALID":
                raise ValueError("Invalid frozen interval")
            lam_thin, lam_thick = float(thin["lambda_s_inv"]), float(thick["lambda_s_inv"])
            if min(lam_thin, lam_thick) <= 0:
                raise ValueError("Nonpositive rate")
            rate_factor = (lam_thick / lam_thin) * (1 - delta) / (1 + delta) / gamma
            if rate_factor <= 1:
                raise ValueError("Shared external monotonicity not proved")
            output.append({
                "case_id": case["case_id"], "window": window, "family": family,
                "scenario": case["scenario"], "delta_rho": delta,
                "external": case["external"], "Gamma": gamma,
                "mobility_sensitivity_fraction": sensitivity,
                "lambda_B65_s_inv": lam_thin, "lambda_B2600_s_inv": lam_thick,
                "L_B65_max_nm": high_thin, "L_B2600_min_nm": low_thick,
                "lower_b_B2600_over_B65": rate_factor * low_thick / high_thin,
                "bound_kind": "CONSERVATIVE_NECESSARY_LOWER_ONLY",
                "numeric_scope": "ARCHIVED_PB_OPERATOR_T12_UNRESOLVED",
            })
    if len(output) != 16 or len({(r["case_id"], r["window"]) for r in output}) != 16:
        raise ValueError("Expected sixteen distinct A case-windows")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    print(f"case_windows={len(output)} min_bound={min(r['lower_b_B2600_over_B65'] for r in output):.4f}")


if __name__ == "__main__":
    main()
