from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path


TEMPLATE_PATH = Path(__file__).with_name("liquid_compatibility_template.csv")

FIELDS = [
    "compatibility_record_id",
    "protocol_liquid_identity",
    "protocol_liquid_batch",
    "coupon_id",
    "parent_membrane_id",
    "exposure_duration",
    "duration_unit",
    "temperature",
    "temperature_unit",
    "agitation",
    "measurement_protocol_id",
    "pre_mass",
    "post_mass",
    "mass_unit",
    "pre_length",
    "post_length",
    "pre_width",
    "post_width",
    "dimension_unit",
    "visual_damage",
    "wetting_observation",
    "raw_file",
    "input_level",
    "status",
    "compatibility_threshold_id",
    "decision",
    "decision_basis",
    "approved_by",
    "approval_time",
    "exclusion_reason",
    "notes",
]

INPUT_LEVELS = {"measured", "reported", "provisional", "unknown"}
STATUSES = {"planned", "ready", "excluded"}
DECISIONS = {"compatible", "incompatible", "review", "unknown"}
VISUAL_DAMAGE = {"no", "yes", "unknown"}
NUMERIC_FIELDS = [
    "exposure_duration",
    "temperature",
    "pre_mass",
    "post_mass",
    "pre_length",
    "post_length",
    "pre_width",
    "post_width",
]
REQUIRED_FOR_READY = FIELDS[:-2]


def load_registry(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != FIELDS:
            raise ValueError("Compatibility columns do not match the frozen interface")
        return list(reader)


def validate_registry(rows: list[dict[str, str]]) -> dict[str, object]:
    errors: list[str] = []
    warnings: list[str] = []
    record_ids = [row["compatibility_record_id"].strip() for row in rows]
    if len(record_ids) != len(set(record_ids)):
        errors.append("compatibility_record_id values must be unique")

    status_counts: Counter[str] = Counter()
    ready_compatible_ids: list[str] = []
    for index, row in enumerate(rows, start=2):
        record_id = row["compatibility_record_id"].strip() or f"row {index}"
        status = row["status"].strip()
        decision = row["decision"].strip()
        if row["input_level"].strip() not in INPUT_LEVELS:
            errors.append(f"{record_id}: invalid input_level")
        if status not in STATUSES:
            errors.append(f"{record_id}: invalid status {status!r}")
        else:
            status_counts[status] += 1
        if decision not in DECISIONS:
            errors.append(f"{record_id}: invalid decision {decision!r}")
        if row["visual_damage"].strip() not in VISUAL_DAMAGE:
            errors.append(f"{record_id}: invalid visual_damage")

        blank_fields = [field for field in FIELDS[:-1] if not row[field].strip()]
        if blank_fields:
            errors.append(f"{record_id}: blank fields {blank_fields}; use 'unknown'")
        if status == "ready":
            unresolved = [
                field
                for field in REQUIRED_FOR_READY
                if row[field].strip().lower() == "unknown"
            ]
            if unresolved:
                errors.append(f"{record_id}: ready row still has unknown required fields")
            if decision == "unknown":
                errors.append(f"{record_id}: ready row requires a decision")
        if status == "excluded" and row["exclusion_reason"].strip().lower() == "unknown":
            errors.append(f"{record_id}: excluded row requires an exclusion_reason")

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
                errors.append(f"{record_id}: {field} must be finite")

        if status == "ready" and decision == "compatible":
            ready_compatible_ids.append(record_id)
        if status == "ready" and row["visual_damage"].strip() == "yes":
            warnings.append(f"{record_id}: visual_damage=yes requires review")

    return {
        "valid": not errors,
        "row_count": len(rows),
        "status_counts": dict(status_counts),
        "ready_compatible_ids": sorted(ready_compatible_ids),
        "warnings": warnings,
        "errors": errors,
        "compatibility_ready": bool(ready_compatible_ids) and not errors,
        "experiment_authorized": False,
        "authorization_note": (
            "A compatible decision applies only to the recorded liquid, batch, exposure, "
            "coupon, and acceptance basis. It does not authorize O/S/W experiments."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate protocol-liquid compatibility records.")
    parser.add_argument("path", nargs="?", type=Path, default=TEMPLATE_PATH)
    args = parser.parse_args()
    report = validate_registry(load_registry(args.path))
    print(json.dumps(report, ensure_ascii=False, indent=2))
    if not report["valid"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
