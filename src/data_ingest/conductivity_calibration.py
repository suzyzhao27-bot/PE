from __future__ import annotations

import argparse
import csv
import json
import math
from bisect import bisect_left
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


TEMPLATE_PATH = Path(__file__).with_name("conductivity_calibration_template.csv")

FIELDS = [
    "calibration_id",
    "standard_id",
    "concentration",
    "concentration_unit",
    "conductivity",
    "conductivity_unit",
    "temperature",
    "temperature_unit",
    "temperature_compensation",
    "input_level",
    "source_record_reference",
    "notes",
]

INPUT_LEVELS = {"measured", "reported", "provisional", "unknown"}


@dataclass(frozen=True)
class CalibrationPoint:
    concentration: float
    conductivity: float


def validate_points(points: list[CalibrationPoint]) -> list[CalibrationPoint]:
    if len(points) < 2:
        raise ValueError("Calibration requires at least two points")
    if any(
        not math.isfinite(point.concentration)
        or not math.isfinite(point.conductivity)
        or point.concentration < 0
        or point.conductivity < 0
        for point in points
    ):
        raise ValueError("Calibration values must be finite and non-negative")

    ordered = sorted(points, key=lambda point: point.conductivity)
    for previous, current in zip(ordered, ordered[1:]):
        if current.conductivity <= previous.conductivity:
            raise ValueError("Conductivity values must be unique and increasing")
        if current.concentration <= previous.concentration:
            raise ValueError("Concentration must increase with conductivity")
    return ordered


def interpolate_concentration(
    conductivity: float,
    points: list[CalibrationPoint],
) -> float:
    ordered = validate_points(points)
    if not math.isfinite(conductivity) or conductivity < 0:
        raise ValueError("Conductivity must be finite and non-negative")
    conductivities = [point.conductivity for point in ordered]
    if conductivity < conductivities[0] or conductivity > conductivities[-1]:
        raise ValueError("Conductivity is outside the calibration range")

    index = bisect_left(conductivities, conductivity)
    if index < len(ordered) and conductivities[index] == conductivity:
        return ordered[index].concentration
    lower = ordered[index - 1]
    upper = ordered[index]
    fraction = (conductivity - lower.conductivity) / (
        upper.conductivity - lower.conductivity
    )
    return lower.concentration + fraction * (
        upper.concentration - lower.concentration
    )


def load_registry(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != FIELDS:
            raise ValueError("Calibration columns do not match the frozen interface")
        return list(reader)


def validate_registry(rows: list[dict[str, str]]) -> dict[str, object]:
    errors: list[str] = []
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for index, row in enumerate(rows, start=2):
        row_name = row["standard_id"].strip() or f"row {index}"
        blank_fields = [field for field in FIELDS[:-1] if not row[field].strip()]
        if blank_fields:
            errors.append(f"{row_name}: blank fields {blank_fields}; use 'unknown'")
        if row["input_level"].strip() not in INPUT_LEVELS:
            errors.append(f"{row_name}: invalid input_level {row['input_level']!r}")
        grouped[row["calibration_id"].strip()].append(row)

    calibration_status: dict[str, dict[str, object]] = {}
    for calibration_id, group in sorted(grouped.items()):
        group_errors: list[str] = []
        for field in (
            "concentration_unit",
            "conductivity_unit",
            "temperature",
            "temperature_unit",
            "temperature_compensation",
        ):
            values = {row[field].strip() for row in group}
            if len(values) != 1:
                group_errors.append(f"mixed {field} values")
        try:
            points = [
                CalibrationPoint(
                    concentration=float(row["concentration"]),
                    conductivity=float(row["conductivity"]),
                )
                for row in group
            ]
            validate_points(points)
        except ValueError as error:
            group_errors.append(str(error))
        if group_errors:
            errors.extend(f"{calibration_id}: {error}" for error in group_errors)
        calibration_status[calibration_id] = {
            "point_count": len(group),
            "valid": not group_errors,
        }

    return {
        "valid": not errors,
        "row_count": len(rows),
        "calibration_count": len(grouped),
        "calibrations": calibration_status,
        "errors": errors,
        "conversion_ready": bool(rows) and not errors,
        "mechanistic_fit_authorized": False,
        "authorization_note": (
            "Calibration validity does not close area, volume, sampling, wetting, blank, "
            "or replicate records. Interpolation never extrapolates."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate conductivity-to-concentration calibration records."
    )
    parser.add_argument("path", nargs="?", type=Path, default=TEMPLATE_PATH)
    args = parser.parse_args()
    report = validate_registry(load_registry(args.path))
    print(json.dumps(report, ensure_ascii=False, indent=2))
    if not report["valid"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
