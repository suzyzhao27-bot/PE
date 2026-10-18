from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path


TEMPLATE_PATH = Path(__file__).with_name("osw_sample_registry_template.csv")

FIELDS = [
    "batch_id",
    "parent_membrane_id",
    "coupon_id",
    "position_x",
    "position_y",
    "machine_direction",
    "intended_use",
    "treatment_group",
    "allocation_seed",
    "allocation_time",
    "protocol_liquid_identity",
    "protocol_liquid_batch",
    "compatibility_record_id",
    "coupon_length",
    "coupon_width",
    "dimension_unit",
    "instrument_or_fixture_id",
    "initial_liquid_volume",
    "membrane_exposed_area",
    "area_volume_ratio",
    "start_time",
    "end_time",
    "temperature",
    "agitation",
    "change_count",
    "change_times",
    "final_measurement_liquid",
    "equilibration_start",
    "equilibration_end",
    "operator",
    "input_level",
    "source_record_reference",
    "status",
    "exclusion_reason",
    "notes",
]

TREATMENTS = {"O", "S", "W"}
INPUT_LEVELS = {"measured", "reported", "provisional", "unknown"}
STATUSES = {"planned", "ready", "excluded"}
REQUIRED_FOR_READY = FIELDS[:-3]


def load_registry(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != FIELDS:
            raise ValueError("Registry columns do not match the frozen O/S/W interface")
        return list(reader)


def validate_registry(rows: list[dict[str, str]]) -> dict[str, object]:
    errors: list[str] = []
    warnings: list[str] = []
    coupon_ids = [row["coupon_id"].strip() for row in rows]

    if len(coupon_ids) != len(set(coupon_ids)):
        errors.append("coupon_id values must be unique")

    groups_by_parent: dict[str, set[str]] = defaultdict(set)
    status_counts: Counter[str] = Counter()
    level_counts: Counter[str] = Counter()
    for index, row in enumerate(rows, start=2):
        coupon_id = row["coupon_id"].strip() or f"row {index}"
        treatment = row["treatment_group"].strip()
        status = row["status"].strip()
        input_level = row["input_level"].strip()

        if treatment not in TREATMENTS:
            errors.append(f"{coupon_id}: invalid treatment_group {treatment!r}")
        elif status != "excluded":
            groups_by_parent[row["parent_membrane_id"].strip()].add(treatment)
        if status not in STATUSES:
            errors.append(f"{coupon_id}: invalid status {status!r}")
        else:
            status_counts[status] += 1
        if input_level not in INPUT_LEVELS:
            errors.append(f"{coupon_id}: invalid input_level {input_level!r}")
        else:
            level_counts[input_level] += 1

        blank_fields = [field for field in FIELDS[:-1] if not row[field].strip()]
        if blank_fields:
            errors.append(f"{coupon_id}: blank fields {blank_fields}; use 'unknown'")

        unresolved = [
            field
            for field in REQUIRED_FOR_READY
            if row[field].strip().lower() == "unknown"
        ]
        if status == "ready" and unresolved:
            errors.append(f"{coupon_id}: ready row still has unknown required fields")
        if status == "excluded" and row["exclusion_reason"].strip().lower() == "unknown":
            errors.append(f"{coupon_id}: excluded row requires an exclusion_reason")

    active_rows = [row for row in rows if row["status"].strip() != "excluded"]
    active_parents = {
        row["parent_membrane_id"].strip()
        for row in active_rows
        if row["parent_membrane_id"].strip().lower() != "unknown"
    }
    if rows and len(active_parents) < 2:
        warnings.append("fewer than two identified active parent membranes")
    for parent in sorted(active_parents):
        missing = sorted(TREATMENTS - groups_by_parent[parent])
        if missing:
            warnings.append(f"{parent}: missing active treatment groups {missing}")

    protocol_ready = (
        bool(rows)
        and not errors
        and not warnings
        and status_counts["ready"] == len(rows)
    )
    return {
        "valid": not errors,
        "row_count": len(rows),
        "parent_membrane_count": len(active_parents),
        "status_counts": dict(status_counts),
        "input_level_counts": dict(level_counts),
        "warnings": warnings,
        "errors": errors,
        "protocol_ready": protocol_ready,
        "experiment_authorized": False,
        "authorization_note": (
            "Registry completeness does not replace instrument, blank, compatibility, "
            "and signed protocol review."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate O/S/W sample registry structure.")
    parser.add_argument("path", nargs="?", type=Path, default=TEMPLATE_PATH)
    args = parser.parse_args()
    report = validate_registry(load_registry(args.path))
    print(json.dumps(report, ensure_ascii=False, indent=2))
    if not report["valid"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
