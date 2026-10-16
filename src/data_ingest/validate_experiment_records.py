from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from normalize_salt_transport import CURVES, OUTPUT_ROOT


TRACKER_PATH = OUTPUT_ROOT / "salt_transport_record_tracker.csv"

FIELDS = [
    "curve_id",
    "status",
    "source_record_reference",
    "parent_membrane_id",
    "coupon_id",
    "time_unit_confirmation",
    "conductivity_unit_confirmation",
    "conductivity_calibration_id",
    "temperature_record_id",
    "effective_area_record_id",
    "left_volume_record_id",
    "right_volume_record_id",
    "sampling_replenishment_record_id",
    "wetting_record_id",
    "mixing_record_id",
    "blank_record_id",
    "independent_replicate_id",
    "exclusion_reason",
    "notes",
]

REQUIRED_FOR_COMPLETE = FIELDS[2:17]
ALLOWED_STATUS = {"pending", "complete", "excluded"}
ALLOWED_CONFIRMATION = {"unknown", "confirmed", "rejected"}


def expected_curve_ids() -> set[str]:
    return {
        f"{curve.project}:{curve.workbook}:{curve.worksheet}:{curve.sample_id}"
        for curve in CURVES
    }


def load_tracker(path: Path = TRACKER_PATH) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != FIELDS:
            raise ValueError("Tracker columns do not match the frozen record interface")
        return list(reader)


def validate_tracker(rows: list[dict[str, str]]) -> dict[str, object]:
    errors: list[str] = []
    curve_ids = [row["curve_id"].strip() for row in rows]
    expected = expected_curve_ids()

    if len(curve_ids) != len(set(curve_ids)):
        errors.append("curve_id values must be unique")
    if set(curve_ids) != expected:
        missing = sorted(expected - set(curve_ids))
        extra = sorted(set(curve_ids) - expected)
        errors.append(f"curve_id coverage mismatch; missing={missing}; extra={extra}")

    status_counts = {status: 0 for status in sorted(ALLOWED_STATUS)}
    unresolved_by_curve: dict[str, list[str]] = {}
    for row in rows:
        curve_id = row["curve_id"].strip()
        status = row["status"].strip()
        if status not in ALLOWED_STATUS:
            errors.append(f"{curve_id}: invalid status {status!r}")
            continue
        status_counts[status] += 1

        blank_fields = [field for field in FIELDS[:-1] if not row[field].strip()]
        if blank_fields:
            errors.append(f"{curve_id}: blank fields {blank_fields}; use 'unknown'")

        for field in ("time_unit_confirmation", "conductivity_unit_confirmation"):
            value = row[field].strip()
            if value not in ALLOWED_CONFIRMATION:
                errors.append(f"{curve_id}: invalid {field} {value!r}")

        unresolved = [
            field
            for field in REQUIRED_FOR_COMPLETE
            if row[field].strip().lower() == "unknown"
        ]
        if unresolved:
            unresolved_by_curve[curve_id] = unresolved

        if status == "complete":
            if unresolved:
                errors.append(f"{curve_id}: complete row still has unknown required fields")
            for field in ("time_unit_confirmation", "conductivity_unit_confirmation"):
                if row[field].strip() != "confirmed":
                    errors.append(f"{curve_id}: complete row requires {field}=confirmed")
        elif status == "excluded" and row["exclusion_reason"].strip().lower() == "unknown":
            errors.append(f"{curve_id}: excluded row requires an exclusion_reason")

    return {
        "valid": not errors,
        "curve_count": len(rows),
        "status_counts": status_counts,
        "unresolved_curve_count": len(unresolved_by_curve),
        "unresolved_by_curve": unresolved_by_curve,
        "errors": errors,
        "mechanistic_fit_authorized": False,
        "authorization_note": (
            "Metadata completeness is necessary but not sufficient. Scientific review, "
            "calibration validity, blanks, and stage-gate approval remain required."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate curve-level experimental record references."
    )
    parser.add_argument(
        "--details",
        action="store_true",
        help="include unresolved field names for every curve",
    )
    args = parser.parse_args()
    report = validate_tracker(load_tracker())
    output = report if args.details else {
        key: value
        for key, value in report.items()
        if key != "unresolved_by_curve"
    }
    print(json.dumps(output, ensure_ascii=False, indent=2))
    if not report["valid"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
