from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path


TEMPLATE_PATH = Path(__file__).with_name("independent_readout_template.csv")

FIELDS = [
    "method_id",
    "readout_role",
    "analyte_or_tracer",
    "charge_state_basis",
    "instrument",
    "calibration_id",
    "matrix",
    "sample_volume",
    "dilution",
    "detection_limit",
    "quantification_limit",
    "blank_id",
    "recovery_check",
    "adsorption_check",
    "stability_check",
    "raw_file",
    "operator",
    "input_level",
    "status",
    "exclusion_reason",
    "notes",
]

ROLES = {
    "independent_na",
    "independent_cl",
    "neutral_tracer",
    "symmetric_conductance",
}
REQUIRED_ROLES = {"independent_na", "independent_cl", "neutral_tracer"}
INPUT_LEVELS = {"measured", "reported", "provisional", "unknown"}
STATUSES = {"planned", "ready", "excluded"}
REQUIRED_FOR_READY = FIELDS[:18]


def load_registry(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != FIELDS:
            raise ValueError("Registry columns do not match the frozen readout interface")
        return list(reader)


def validate_registry(rows: list[dict[str, str]]) -> dict[str, object]:
    errors: list[str] = []
    warnings: list[str] = []
    method_ids = [row["method_id"].strip() for row in rows]
    if len(method_ids) != len(set(method_ids)):
        errors.append("method_id values must be unique")

    status_counts: Counter[str] = Counter()
    role_counts: Counter[str] = Counter()
    active_ready_roles: set[str] = set()
    for index, row in enumerate(rows, start=2):
        method_id = row["method_id"].strip() or f"row {index}"
        role = row["readout_role"].strip()
        status = row["status"].strip()
        input_level = row["input_level"].strip()

        if role not in ROLES:
            errors.append(f"{method_id}: invalid readout_role {role!r}")
        else:
            role_counts[role] += 1
        if status not in STATUSES:
            errors.append(f"{method_id}: invalid status {status!r}")
        else:
            status_counts[status] += 1
        if input_level not in INPUT_LEVELS:
            errors.append(f"{method_id}: invalid input_level {input_level!r}")

        blank_fields = [field for field in FIELDS[:-1] if not row[field].strip()]
        if blank_fields:
            errors.append(f"{method_id}: blank fields {blank_fields}; use 'unknown'")

        unresolved = [
            field
            for field in REQUIRED_FOR_READY
            if row[field].strip().lower() == "unknown"
        ]
        if status == "ready" and unresolved:
            errors.append(f"{method_id}: ready row still has unknown required fields")
        if status == "excluded" and row["exclusion_reason"].strip().lower() == "unknown":
            errors.append(f"{method_id}: excluded row requires an exclusion_reason")

        if role == "neutral_tracer" and status == "ready":
            tracer_checks = [
                "charge_state_basis",
                "recovery_check",
                "adsorption_check",
                "stability_check",
            ]
            if any(row[field].strip().lower() == "unknown" for field in tracer_checks):
                errors.append(f"{method_id}: ready neutral tracer lacks required checks")

        if role in ROLES and status == "ready":
            active_ready_roles.add(role)

    missing_roles = sorted(REQUIRED_ROLES - active_ready_roles)
    if rows and missing_roles:
        warnings.append(f"missing ready minimum readout roles {missing_roles}")

    selection_ready = bool(rows) and not errors and not warnings
    return {
        "valid": not errors,
        "row_count": len(rows),
        "status_counts": dict(status_counts),
        "role_counts": dict(role_counts),
        "warnings": warnings,
        "errors": errors,
        "selection_ready": selection_ready,
        "experiment_authorized": False,
        "authorization_note": (
            "Readout selection still requires instrument availability, blank performance, "
            "sample-volume compatibility, and signed protocol review."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate independent ion and wet-path readout registry."
    )
    parser.add_argument("path", nargs="?", type=Path, default=TEMPLATE_PATH)
    args = parser.parse_args()
    report = validate_registry(load_registry(args.path))
    print(json.dumps(report, ensure_ascii=False, indent=2))
    if not report["valid"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
