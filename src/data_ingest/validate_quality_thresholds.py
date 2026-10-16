from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path


TEMPLATE_PATH = Path(__file__).with_name("quality_threshold_template.csv")

FIELDS = [
    "threshold_id",
    "metric_role",
    "metric_name",
    "comparison",
    "lower_bound",
    "upper_bound",
    "unit",
    "basis_record_ids",
    "action",
    "input_level",
    "status",
    "source_record_reference",
    "approved_by",
    "approval_time",
    "notes",
]

METRIC_ROLES = {
    "calibration_drift",
    "mass_closure_error",
    "no_membrane_blank",
    "impermeable_blank",
    "temperature_deviation",
}
REQUIRED_ROLES = set(METRIC_ROLES)
COMPARISONS = {"abs_le", "le", "ge", "between"}
ACTIONS = {"review", "exclude"}
INPUT_LEVELS = {"measured", "reported", "provisional"}
STATUSES = {"planned", "ready", "retired"}
REQUIRED_FOR_READY = FIELDS[:14]


def load_registry(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != FIELDS:
            raise ValueError("Threshold columns do not match the frozen interface")
        return list(reader)


def parse_bound(value: str) -> float | None:
    normalized = value.strip().lower()
    if normalized in {"unknown", "not_applicable"}:
        return None
    result = float(value)
    if not math.isfinite(result):
        raise ValueError
    return result


def evaluate_threshold(
    value: float,
    comparison: str,
    lower_bound: float | None,
    upper_bound: float | None,
) -> bool:
    if not math.isfinite(value):
        raise ValueError("Metric value must be finite")
    if comparison == "abs_le":
        if upper_bound is None:
            raise ValueError("abs_le requires upper_bound")
        return abs(value) <= upper_bound
    if comparison == "le":
        if upper_bound is None:
            raise ValueError("le requires upper_bound")
        return value <= upper_bound
    if comparison == "ge":
        if lower_bound is None:
            raise ValueError("ge requires lower_bound")
        return value >= lower_bound
    if comparison == "between":
        if lower_bound is None or upper_bound is None:
            raise ValueError("between requires both bounds")
        if lower_bound > upper_bound:
            raise ValueError("lower_bound exceeds upper_bound")
        return lower_bound <= value <= upper_bound
    raise ValueError(f"Unknown comparison {comparison!r}")


def evaluate_metrics(
    metrics: dict[str, float],
    threshold_rows: list[dict[str, str]],
) -> list[dict[str, object]]:
    report = validate_registry(threshold_rows)
    if not report["threshold_set_ready"]:
        raise ValueError(f"Threshold registry is not ready: {report['errors']}")
    results: list[dict[str, object]] = []
    for row in threshold_rows:
        if row["status"] != "ready":
            continue
        role = row["metric_role"]
        if role not in metrics:
            raise ValueError(f"Missing metric value for {role}")
        passed = evaluate_threshold(
            metrics[role],
            row["comparison"],
            parse_bound(row["lower_bound"]),
            parse_bound(row["upper_bound"]),
        )
        results.append(
            {
                "threshold_id": row["threshold_id"],
                "metric_role": role,
                "value": metrics[role],
                "unit": row["unit"],
                "passed": passed,
                "decision": "pass" if passed else row["action"],
                "experiment_authorized": False,
            }
        )
    return results


def validate_registry(rows: list[dict[str, str]]) -> dict[str, object]:
    errors: list[str] = []
    warnings: list[str] = []
    threshold_ids = [row["threshold_id"].strip() for row in rows]
    if len(threshold_ids) != len(set(threshold_ids)):
        errors.append("threshold_id values must be unique")

    ready_roles: set[str] = set()
    status_counts: Counter[str] = Counter()
    for index, row in enumerate(rows, start=2):
        threshold_id = row["threshold_id"].strip() or f"row {index}"
        role = row["metric_role"].strip()
        comparison = row["comparison"].strip()
        status = row["status"].strip()

        if role not in METRIC_ROLES:
            errors.append(f"{threshold_id}: invalid metric_role {role!r}")
        if comparison not in COMPARISONS:
            errors.append(f"{threshold_id}: invalid comparison {comparison!r}")
        if row["action"].strip() not in ACTIONS:
            errors.append(f"{threshold_id}: invalid action")
        if row["input_level"].strip() not in INPUT_LEVELS:
            errors.append(f"{threshold_id}: invalid input_level")
        if status not in STATUSES:
            errors.append(f"{threshold_id}: invalid status {status!r}")
        else:
            status_counts[status] += 1

        blank_fields = [field for field in FIELDS[:-1] if not row[field].strip()]
        if blank_fields:
            errors.append(f"{threshold_id}: blank fields {blank_fields}; use 'unknown'")
        if status == "ready":
            unresolved = [
                field
                for field in REQUIRED_FOR_READY
                if row[field].strip().lower() == "unknown"
            ]
            if unresolved:
                errors.append(f"{threshold_id}: ready row still has unknown required fields")

        try:
            lower = parse_bound(row["lower_bound"])
            upper = parse_bound(row["upper_bound"])
        except ValueError:
            errors.append(f"{threshold_id}: bounds must be finite numbers or declared missing")
            lower = None
            upper = None

        if status == "ready" and comparison in COMPARISONS:
            if comparison in {"abs_le", "le"} and upper is None:
                errors.append(f"{threshold_id}: {comparison} requires upper_bound")
            if comparison == "ge" and lower is None:
                errors.append(f"{threshold_id}: ge requires lower_bound")
            if comparison == "between":
                if lower is None or upper is None:
                    errors.append(f"{threshold_id}: between requires both bounds")
                elif lower > upper:
                    errors.append(f"{threshold_id}: lower_bound exceeds upper_bound")
        if status == "ready" and role in METRIC_ROLES:
            ready_roles.add(role)

    missing_roles = sorted(REQUIRED_ROLES - ready_roles)
    if rows and missing_roles:
        warnings.append(f"missing ready quality-threshold roles {missing_roles}")

    return {
        "valid": not errors,
        "row_count": len(rows),
        "status_counts": dict(status_counts),
        "ready_roles": sorted(ready_roles),
        "warnings": warnings,
        "errors": errors,
        "threshold_set_ready": bool(rows) and not errors and not missing_roles,
        "experiment_authorized": False,
        "authorization_note": (
            "Threshold registration does not authorize experiments. Thresholds must be "
            "frozen before viewing treatment differences and applied with recorded blanks."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate predeclared quality thresholds.")
    parser.add_argument("path", nargs="?", type=Path, default=TEMPLATE_PATH)
    args = parser.parse_args()
    report = validate_registry(load_registry(args.path))
    print(json.dumps(report, ensure_ascii=False, indent=2))
    if not report["valid"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
