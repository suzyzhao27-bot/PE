from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path


TEMPLATE_PATH = Path(__file__).with_name("blank_registry_template.csv")

FIELDS = [
    "blank_id",
    "blank_type",
    "applies_to",
    "fixture_id",
    "calibration_id",
    "left_solution",
    "right_solution",
    "left_volume",
    "right_volume",
    "volume_unit",
    "temperature",
    "temperature_unit",
    "agitation_left",
    "agitation_right",
    "duration",
    "duration_unit",
    "sampling_plan",
    "raw_file",
    "input_level",
    "status",
    "exclusion_reason",
    "notes",
]

BLANK_TYPES = {"no_membrane", "impermeable", "liquid_container", "electrode_bridge"}
APPLICATIONS = {"salt_transport", "potential", "both"}
INPUT_LEVELS = {"measured", "reported", "provisional", "unknown"}
STATUSES = {"planned", "ready", "excluded"}
REQUIRED_FOR_READY = FIELDS[:20]


def load_registry(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != FIELDS:
            raise ValueError("Blank-registry columns do not match the frozen interface")
        return list(reader)


def validate_registry(rows: list[dict[str, str]]) -> dict[str, object]:
    errors: list[str] = []
    warnings: list[str] = []
    blank_ids = [row["blank_id"].strip() for row in rows]
    if len(blank_ids) != len(set(blank_ids)):
        errors.append("blank_id values must be unique")

    ready_types_by_application: dict[str, set[str]] = {
        "salt_transport": set(),
        "potential": set(),
    }
    status_counts: Counter[str] = Counter()
    type_counts: Counter[str] = Counter()
    for index, row in enumerate(rows, start=2):
        blank_id = row["blank_id"].strip() or f"row {index}"
        blank_type = row["blank_type"].strip()
        application = row["applies_to"].strip()
        status = row["status"].strip()
        input_level = row["input_level"].strip()

        if blank_type not in BLANK_TYPES:
            errors.append(f"{blank_id}: invalid blank_type {blank_type!r}")
        else:
            type_counts[blank_type] += 1
        if application not in APPLICATIONS:
            errors.append(f"{blank_id}: invalid applies_to {application!r}")
        if status not in STATUSES:
            errors.append(f"{blank_id}: invalid status {status!r}")
        else:
            status_counts[status] += 1
        if input_level not in INPUT_LEVELS:
            errors.append(f"{blank_id}: invalid input_level {input_level!r}")

        blank_fields = [field for field in FIELDS[:-1] if not row[field].strip()]
        if blank_fields:
            errors.append(f"{blank_id}: blank fields {blank_fields}; use 'unknown'")
        if status == "ready":
            unresolved = [
                field
                for field in REQUIRED_FOR_READY
                if row[field].strip().lower() == "unknown"
            ]
            if unresolved:
                errors.append(f"{blank_id}: ready row still has unknown required fields")
        if status == "excluded" and row["exclusion_reason"].strip().lower() == "unknown":
            errors.append(f"{blank_id}: excluded row requires an exclusion_reason")

        for field in ("left_volume", "right_volume", "duration"):
            if row[field].strip().lower() == "unknown":
                continue
            try:
                value = float(row[field])
                if not math.isfinite(value) or value < 0:
                    raise ValueError
            except ValueError:
                errors.append(f"{blank_id}: {field} must be finite and non-negative")
        if row["temperature"].strip().lower() != "unknown":
            try:
                if not math.isfinite(float(row["temperature"])):
                    raise ValueError
            except ValueError:
                errors.append(f"{blank_id}: temperature must be numeric")

        if status == "ready" and blank_type in BLANK_TYPES:
            if application in {"salt_transport", "both"}:
                ready_types_by_application["salt_transport"].add(blank_type)
            if application in {"potential", "both"}:
                ready_types_by_application["potential"].add(blank_type)

    salt_required = {"no_membrane", "impermeable"}
    salt_missing = sorted(salt_required - ready_types_by_application["salt_transport"])
    potential_missing = sorted(
        {"electrode_bridge"} - ready_types_by_application["potential"]
    )
    if rows and salt_missing:
        warnings.append(f"missing ready salt-transport blanks {salt_missing}")
    if rows and potential_missing:
        warnings.append(f"missing ready potential blanks {potential_missing}")
    if rows and "liquid_container" not in ready_types_by_application["salt_transport"]:
        warnings.append("no ready liquid/container blank for salt transport")

    return {
        "valid": not errors,
        "row_count": len(rows),
        "status_counts": dict(status_counts),
        "blank_type_counts": dict(type_counts),
        "warnings": warnings,
        "errors": errors,
        "salt_blank_set_ready": bool(rows) and not errors and not salt_missing,
        "potential_blank_set_ready": bool(rows) and not errors and not potential_missing,
        "experiment_authorized": False,
        "authorization_note": (
            "Blank registration does not establish numerical acceptance thresholds or "
            "authorize experiments. Blank results must still be reviewed."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate experimental blank registry.")
    parser.add_argument("path", nargs="?", type=Path, default=TEMPLATE_PATH)
    args = parser.parse_args()
    report = validate_registry(load_registry(args.path))
    print(json.dumps(report, ensure_ascii=False, indent=2))
    if not report["valid"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
