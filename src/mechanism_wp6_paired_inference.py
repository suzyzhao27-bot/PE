"""Strict paired-log WP6 inference for three quality-passed parent films."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


T95_DF2 = math.sqrt(2 * .95**2 / (1 - .95**2))
N = 3
PAIR_MAP = (Path(__file__).resolve().parents[1] /
            "data/processed/salt_transport/22_WP6母膜B_E配对运行入口_20261029.csv")


def mean(values: list[float]) -> float:
    return sum(values) / len(values)


def covariance(x: list[float], y: list[float]) -> float:
    mx, my = mean(x), mean(y)
    return sum((a - mx) * (b - my) for a, b in zip(x, y, strict=True)) / (N - 1)


def log_interval(values: list[float]) -> dict[str, float]:
    center = mean(values)
    radius = T95_DF2 * math.sqrt(covariance(values, values) / N)
    return {"mean_log": center, "sample_sd_log": math.sqrt(covariance(values, values)),
            "lower_log": center - radius, "upper_log": center + radius,
            "geometric_mean_ratio": math.exp(center),
            "lower_ratio": math.exp(center - radius), "upper_ratio": math.exp(center + radius)}


def fieller_ratio(numerator: list[float], denominator: list[float]) -> dict[str, float | str]:
    u, s = mean(numerator), mean(denominator)
    q = T95_DF2**2 / N
    a = s * s - q * covariance(denominator, denominator)
    b = -2 * (s * u - q * covariance(denominator, numerator))
    c = u * u - q * covariance(numerator, numerator)
    if a <= 0:
        return {"status": "UNBOUNDED_OR_DISJOINT_DENOMINATOR_NOT_SEPARATED_FROM_ZERO"}
    disc = b * b - 4 * a * c
    if disc < -1e-12 * max(1, b * b, abs(4 * a * c)):
        raise ArithmeticError("Invalid Fieller discriminant")
    width = math.sqrt(max(0.0, disc))
    return {"status": "BOUNDED_95PCT_FIELLER", "estimate": u / s,
            "lower": (-b - width) / (2 * a), "upper": (-b + width) / (2 * a)}


def validate_fixed_pairs(rows: list[dict[str, str]]) -> tuple[str, str, str]:
    with PAIR_MAP.open(encoding="utf-8-sig", newline="") as stream:
        fixed = list(csv.DictReader(stream))
    if len(fixed) != 15:
        raise ArithmeticError("Fixed WP6 pair map changed")
    by_runs = {(r["B_run_id"], r["E_run_id"]): r for r in fixed}
    matched = []
    for row in rows:
        pair = by_runs.get((row.get("B_run_id"), row.get("E_run_id")))
        if pair is None or row.get("parent_film_id") != pair["parent_film_id"]:
            raise ValueError("B/E run IDs do not match the fixed mother-film pair map")
        matched.append(pair)
    conditions = {(r["thickness_level"], r["salt_level"]) for r in matched}
    if len(conditions) != 1:
        raise ValueError("Three mother films must belong to one fixed condition")
    level, salt = next(iter(conditions))
    expected = {r["parent_film_id"] for r in fixed
                if (r["thickness_level"], r["salt_level"]) == (level, salt)}
    if {r["parent_film_id"] for r in matched} != expected:
        raise ValueError("P1/P2/P3 coverage of the fixed condition is incomplete")
    if salt == "low":
        try:
            choices = {float(r["global_low_salt_choice_mM"]) for r in rows}
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError("One global measured-LOQ low-salt choice is required") from exc
        if choices not in ({1.0}, {10.0}):
            raise ValueError("One global measured-LOQ low-salt choice is required")
        choice = str(int(next(iter(choices))))
    else:
        if any(r.get("global_low_salt_choice_mM") not in (None, "") for r in rows):
            raise ValueError("High-salt condition must not declare a low-salt choice")
        choice = "NA_HIGH_600_PLAN"
    return level, salt, choice


def analyze(rows: list[dict[str, str]]) -> dict[str, object]:
    if len(rows) != N or len({r.get("parent_film_id") for r in rows}) != N:
        raise ValueError("Exactly three distinct parent films are required")
    if any(r.get("quality_status") != "PASS" or r.get("metric_scope") != "MEMBRANE_P"
           or not r.get("source_record") for r in rows):
        raise ValueError("All rows need quality PASS, MEMBRANE_P scope and source record")
    level, salt_level, choice = validate_fixed_pairs(rows)
    fields = ("B_salt_P_m_s", "E_salt_P_m_s", "B_urea_P_m_s", "E_urea_P_m_s")
    values = []
    for row in rows:
        try:
            v = [float(row[field]) for field in fields]
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError("Missing or censored membrane permeability") from exc
        if any(not math.isfinite(x) or x <= 0 for x in v):
            raise ValueError("All membrane permeabilities must be finite and positive")
        values.append(v)
    salt = [math.log(v[1] / v[0]) for v in values]
    urea = [math.log(v[3] / v[2]) for v in values]
    pi = [s - u for s, u in zip(salt, urea, strict=True)]
    salt_interval = log_interval(salt)
    # The fixed plan defines F_common only when E_s is clearly above one.
    # A positive 95% lower log bound is a conservative interpretation gate.
    f_common = (fieller_ratio(urea, salt) if salt_interval["lower_log"] > 0 else
                {"status": "NOT_INTERPRETABLE_ES_95PCT_LOWER_NOT_ABOVE_ONE"})
    result = {
        "status": "METHOD_RESULT_REQUIRES_SOURCE_REVIEW",
        "n_independent_parent_films": N,
        "parent_film_ids": [r["parent_film_id"] for r in rows],
        "fixed_condition": {"thickness_level": level, "salt_level": salt_level,
                            "global_low_salt_choice_mM": choice},
        "salt_E_over_B": salt_interval,
        "urea_E_over_B": log_interval(urea),
        "Pi_paired": log_interval(pi),
        "F_common": f_common,
        "note": "No LOQ censoring, interval-valued P, failed quality gate, or external-resistance ambiguity accepted.",
    }
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.input.open(encoding="utf-8-sig", newline="") as stream:
        result = analyze(list(csv.DictReader(stream)))
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(result["status"])


if __name__ == "__main__":
    main()
