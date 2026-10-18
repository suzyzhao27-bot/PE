from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path


TEMPLATE_PATH = Path(__file__).with_name("liquid_compatibility_threshold_template.csv")

FIELDS = [
    "compatibility_threshold_id",
    "protocol_liquid_identity",
    "protocol_liquid_batch",
    "exposure_duration",
    "duration_unit",
    "temperature",
    "temperature_unit",
    "agitation",
    "measurement_protocol_id",
    "max_abs_mass_change_pct",
    "max_abs_length_change_pct",
    "max_abs_width_change_pct",
    "visual_damage_allowed",
    "failure_decision",
    "input_level",
    "status",
    "basis_record_ids",
    "approved_by",
    "approval_time",
    "notes",
]

NUMERIC_FIELDS = [
    "exposure_duration",
    "temperature",
    "max_abs_mass_change_pct",
    "max_abs_length_change_pct",
    "max_abs_width_change_pct",
]
INPUT_LEVELS = {"measured", "reported", "provisional"}
STATUSES = {"planned", "ready", "retired"}
FAILURE_DECISIONS = {"review", "incompatible"}
BOOLEAN_VALUES = {"yes", "no"}
REQUIRED_FOR_READY = FIELDS[:-1]


def load_registry(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != FIELDS:
            raise ValueError("Compatibility-threshold columns do not match the frozen interface")
        return list(reader)


def validate_registry(rows: list[dict[str, str]]) -> dict[str, object]:
    errors: list[str] = []
    threshold_ids = [row["compatibility_threshold_id"].strip() for row in rows]
    if len(threshold_ids) != len(set(threshold_ids)):
        errors.append("compatibility_threshold_id values must be unique")

    status_counts: Counter[str] = Counter()
    ready_threshold_ids: list[str] = []
    for index, row in enumerate(rows, start=2):
        threshold_id = row["compatibility_threshold_id"].strip() or f"row {index}"
        status = row["status"].strip()
        if row["input_level"].strip() not in INPUT_LEVELS:
            errors.append(f"{threshold_id}: invalid input_level")
        if status not in STATUSES:
            errors.append(f"{threshold_id}: invalid status {status!r}")
        else:
            status_counts[status] += 1
        if row["failure_decision"].strip() not in FAILURE_DECISIONS:
            errors.append(f"{threshold_id}: invalid failure_decision")
        if row["visual_damage_allowed"].strip() not in BOOLEAN_VALUES:
            errors.append(f"{threshold_id}: invalid visual_damage_allowed")

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

        for field in NUMERIC_FIELDS:
            raw = row[field].strip().lower()
            if raw == "unknown":
                continue
            try:
                value = float(raw)
                if not math.isfinite(value):
                    raise ValueError
                if field != "temperature" and value < 0:
                    raise ValueError
            except ValueError:
                errors.append(f"{threshold_id}: {field} must be finite and nonnegative")

        if status == "ready":
            ready_threshold_ids.append(threshold_id)

    return {
        "valid": not errors,
        "row_count": len(rows),
        "status_counts": dict(status_counts),
        "ready_threshold_ids": sorted(ready_threshold_ids),
        "errors": errors,
        "thresholds_ready": bool(ready_threshold_ids) and not errors,
        "experiment_authorized": False,
        "authorization_note": (
            "Compatibility thresholds are condition-specific and must be frozen before "
            "post-exposure measurements are reviewed."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate predeclared protocol-liquid compatibility thresholds."
    )
    parser.add_argument("path", nargs="?", type=Path, default=TEMPLATE_PATH)
    args = parser.parse_args()
    report = validate_registry(load_registry(args.path))
    print(json.dumps(report, ensure_ascii=False, indent=2))
    if not report["valid"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
