from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from validate_liquid_compatibility import (
    TEMPLATE_PATH as RECORD_TEMPLATE_PATH,
    load_registry as load_record_registry,
    validate_registry as validate_record_registry,
)
from validate_liquid_compatibility_thresholds import (
    TEMPLATE_PATH as THRESHOLD_TEMPLATE_PATH,
    load_registry as load_threshold_registry,
    validate_registry as validate_threshold_registry,
)


CONDITION_FIELDS = [
    "protocol_liquid_identity",
    "protocol_liquid_batch",
    "exposure_duration",
    "duration_unit",
    "temperature",
    "temperature_unit",
    "agitation",
    "measurement_protocol_id",
]


def absolute_percent_change(pre_value: str, post_value: str) -> float:
    pre = float(pre_value)
    post = float(post_value)
    if not math.isfinite(pre) or not math.isfinite(post) or pre <= 0 or post < 0:
        raise ValueError("Pre value must be positive and post value nonnegative")
    return abs(post - pre) / pre * 100.0


def evaluate_records(
    record_rows: list[dict[str, str]],
    threshold_rows: list[dict[str, str]],
) -> dict[str, object]:
    record_report = validate_record_registry(record_rows)
    threshold_report = validate_threshold_registry(threshold_rows)
    errors = list(record_report["errors"]) + list(threshold_report["errors"])
    ready_thresholds = {
        row["compatibility_threshold_id"].strip(): row
        for row in threshold_rows
        if row["status"].strip() == "ready"
    }
    results: list[dict[str, object]] = []
    evaluated_ready_compatible_ids: list[str] = []

    for row in record_rows:
        if row["status"].strip() != "ready":
            continue
        record_id = row["compatibility_record_id"].strip()
        threshold_id = row["compatibility_threshold_id"].strip()
        threshold = ready_thresholds.get(threshold_id)
        if threshold is None:
            errors.append(f"{record_id}: referenced compatibility threshold is not ready")
            continue
        mismatched = [
            field
            for field in CONDITION_FIELDS
            if row[field].strip() != threshold[field].strip()
        ]
        if mismatched:
            errors.append(f"{record_id}: threshold condition mismatch {mismatched}")
            continue

        try:
            metrics = {
                "abs_mass_change_pct": absolute_percent_change(
                    row["pre_mass"], row["post_mass"]
                ),
                "abs_length_change_pct": absolute_percent_change(
                    row["pre_length"], row["post_length"]
                ),
                "abs_width_change_pct": absolute_percent_change(
                    row["pre_width"], row["post_width"]
                ),
            }
            checks = {
                "mass_change": metrics["abs_mass_change_pct"]
                <= float(threshold["max_abs_mass_change_pct"]),
                "length_change": metrics["abs_length_change_pct"]
                <= float(threshold["max_abs_length_change_pct"]),
                "width_change": metrics["abs_width_change_pct"]
                <= float(threshold["max_abs_width_change_pct"]),
                "visual_damage": row["visual_damage"].strip() == "no"
                or threshold["visual_damage_allowed"].strip() == "yes",
            }
        except ValueError as exc:
            errors.append(f"{record_id}: cannot evaluate compatibility metrics: {exc}")
            continue

        passed = all(checks.values())
        derived_decision = (
            "compatible" if passed else threshold["failure_decision"].strip()
        )
        decision_matches = row["decision"].strip() == derived_decision
        results.append(
            {
                "compatibility_record_id": record_id,
                "compatibility_threshold_id": threshold_id,
                **metrics,
                "checks": checks,
                "passed": passed,
                "derived_decision": derived_decision,
                "recorded_decision": row["decision"].strip(),
                "decision_matches": decision_matches,
                "experiment_authorized": False,
            }
        )
        if passed and decision_matches:
            evaluated_ready_compatible_ids.append(record_id)
        if not decision_matches:
            errors.append(
                f"{record_id}: recorded decision does not match predeclared thresholds"
            )

    return {
        "valid": not errors,
        "record_count": len(record_rows),
        "evaluated_count": len(results),
        "evaluated_ready_compatible_ids": sorted(evaluated_ready_compatible_ids),
        "results": results,
        "errors": errors,
        "compatibility_ready": bool(evaluated_ready_compatible_ids) and not errors,
        "experiment_authorized": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Evaluate liquid-compatibility records against predeclared thresholds."
    )
    parser.add_argument("records", nargs="?", type=Path, default=RECORD_TEMPLATE_PATH)
    parser.add_argument(
        "thresholds", nargs="?", type=Path, default=THRESHOLD_TEMPLATE_PATH
    )
    args = parser.parse_args()
    report = evaluate_records(
        load_record_registry(args.records), load_threshold_registry(args.thresholds)
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))
    if not report["valid"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
