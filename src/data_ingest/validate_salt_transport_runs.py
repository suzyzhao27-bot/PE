from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path


TEMPLATE_PATH = Path(__file__).with_name("salt_transport_run_template.csv")

FIELDS = [
    "run_id",
    "coupon_id",
    "parent_membrane_id",
    "treatment_group",
    "receiver_side",
    "time",
    "time_unit",
    "record_phase",
    "left_conductivity",
    "right_conductivity",
    "conductivity_unit",
    "calibration_id",
    "left_volume",
    "right_volume",
    "volume_unit",
    "left_sampled_volume",
    "right_sampled_volume",
    "left_replenished_volume",
    "right_replenished_volume",
    "left_replenishment_concentration",
    "right_replenishment_concentration",
    "concentration_unit",
    "temperature",
    "temperature_unit",
    "effective_area",
    "area_unit",
    "blank_id",
    "fixture_id",
    "leakage_flag",
    "input_level",
    "source_record_reference",
    "notes",
]

TREATMENTS = {"O", "S", "W", "control"}
INPUT_LEVELS = {"measured", "reported", "provisional", "unknown"}
LEAKAGE_FLAGS = {"no", "yes", "unknown"}
NUMERIC_NONNEGATIVE = [
    "time",
    "left_conductivity",
    "right_conductivity",
    "left_volume",
    "right_volume",
    "left_sampled_volume",
    "right_sampled_volume",
    "left_replenished_volume",
    "right_replenished_volume",
    "left_replenishment_concentration",
    "right_replenishment_concentration",
]
CONSTANT_WITHIN_RUN = [
    "coupon_id",
    "parent_membrane_id",
    "treatment_group",
    "receiver_side",
    "time_unit",
    "record_phase",
    "conductivity_unit",
    "calibration_id",
    "volume_unit",
    "concentration_unit",
    "temperature_unit",
    "effective_area",
    "area_unit",
    "blank_id",
    "fixture_id",
]


def load_registry(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != FIELDS:
            raise ValueError("Run-log columns do not match the frozen interface")
        return list(reader)


def validate_registry(rows: list[dict[str, str]]) -> dict[str, object]:
    errors: list[str] = []
    warnings: list[str] = []
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)

    for index, row in enumerate(rows, start=2):
        label = f"row {index}"
        blank_fields = [field for field in FIELDS[:-1] if not row[field].strip()]
        if blank_fields:
            errors.append(f"{label}: blank fields {blank_fields}; use 'unknown'")
        if row["treatment_group"].strip() not in TREATMENTS:
            errors.append(f"{label}: invalid treatment_group")
        if row["receiver_side"].strip() not in {"left", "right"}:
            errors.append(f"{label}: receiver_side must be 'left' or 'right'")
        if row["input_level"].strip() not in INPUT_LEVELS:
            errors.append(f"{label}: invalid input_level")
        if row["leakage_flag"].strip() not in LEAKAGE_FLAGS:
            errors.append(f"{label}: invalid leakage_flag")
        if row["record_phase"].strip() != "pre_event":
            errors.append(f"{label}: record_phase must be 'pre_event'")
        for field in NUMERIC_NONNEGATIVE:
            try:
                value = float(row[field])
                if not math.isfinite(value) or value < 0:
                    raise ValueError
            except ValueError:
                errors.append(f"{label}: {field} must be finite and non-negative")
        try:
            temperature = float(row["temperature"])
            effective_area = float(row["effective_area"])
            if not math.isfinite(temperature):
                raise ValueError
            if not math.isfinite(effective_area) or effective_area <= 0:
                errors.append(f"{label}: effective_area must be finite and positive")
        except ValueError:
            errors.append(f"{label}: temperature and effective_area must be numeric")
        grouped[row["run_id"].strip()].append(row)

    run_status: dict[str, dict[str, object]] = {}
    for run_id, group in sorted(grouped.items()):
        run_errors: list[str] = []
        if not run_id or run_id.lower() == "unknown":
            run_errors.append("run_id must be identified")
        for field in CONSTANT_WITHIN_RUN:
            values = {row[field].strip() for row in group}
            if len(values) != 1:
                run_errors.append(f"mixed {field} values")
        try:
            times = [float(row["time"]) for row in group]
            if len(times) != len(set(times)):
                run_errors.append("duplicate time values")
            if times != sorted(times):
                run_errors.append("rows are not in increasing time order")
        except ValueError:
            pass
        if any(row["leakage_flag"].strip() == "yes" for row in group):
            warnings.append(f"{run_id}: leakage_flag=yes requires review")
        if run_errors:
            errors.extend(f"{run_id}: {error}" for error in run_errors)
        run_status[run_id] = {
            "point_count": len(group),
            "valid": not run_errors,
        }

    return {
        "valid": not errors,
        "row_count": len(rows),
        "run_count": len(grouped),
        "runs": run_status,
        "warnings": warnings,
        "errors": errors,
        "run_log_ready": bool(rows) and not errors,
        "mechanistic_fit_authorized": False,
        "authorization_note": (
            "Run-log validity does not prove calibration, blank, mass closure, or sample "
            "independence. leakage_flag=yes always requires review."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate salt-transport run logs.")
    parser.add_argument("path", nargs="?", type=Path, default=TEMPLATE_PATH)
    args = parser.parse_args()
    report = validate_registry(load_registry(args.path))
    print(json.dumps(report, ensure_ascii=False, indent=2))
    if not report["valid"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
