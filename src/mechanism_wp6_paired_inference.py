"""Strict paired-log WP6 inference for three quality-passed parent films."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


T95_DF2 = math.sqrt(2 * .95**2 / (1 - .95**2))
N = 3


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


def analyze(rows: list[dict[str, str]]) -> dict[str, object]:
    if len(rows) != N or len({r.get("parent_film_id") for r in rows}) != N:
        raise ValueError("Exactly three distinct parent films are required")
    if any(r.get("quality_status") != "PASS" or r.get("metric_scope") != "MEMBRANE_P"
           or not r.get("source_record") for r in rows):
        raise ValueError("All rows need quality PASS, MEMBRANE_P scope and source record")
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
    result = {
        "status": "METHOD_RESULT_REQUIRES_SOURCE_REVIEW",
        "n_independent_parent_films": N,
        "parent_film_ids": [r["parent_film_id"] for r in rows],
        "salt_E_over_B": log_interval(salt),
        "urea_E_over_B": log_interval(urea),
        "Pi_paired": log_interval(pi),
        "F_common": fieller_ratio(urea, salt),
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
