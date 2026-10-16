from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

from conductivity_calibration import (
    CalibrationPoint,
    interpolate_concentration,
    validate_registry as validate_calibration_registry,
)
from mass_balance import interval_flux, relative_closure_error
from validate_salt_transport_runs import validate_registry as validate_run_registry


CONCENTRATION_TO_MOL_M3 = {"mol/L": 1000.0, "mmol/L": 1.0}
VOLUME_TO_M3 = {"L": 1e-3, "mL": 1e-6}
AREA_TO_M2 = {"m2": 1.0, "cm2": 1e-4}
TIME_TO_S = {"s": 1.0, "min": 60.0, "h": 3600.0}


def group_calibrations(rows: list[dict[str, str]]) -> dict[str, dict[str, object]]:
    report = validate_calibration_registry(rows)
    if not report["conversion_ready"]:
        raise ValueError(f"Calibration registry is not ready: {report['errors']}")
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["calibration_id"]].append(row)

    result: dict[str, dict[str, object]] = {}
    for calibration_id, group in grouped.items():
        result[calibration_id] = {
            "points": [
                CalibrationPoint(float(row["concentration"]), float(row["conductivity"]))
                for row in group
            ],
            "concentration_unit": group[0]["concentration_unit"],
            "conductivity_unit": group[0]["conductivity_unit"],
            "temperature": float(group[0]["temperature"]),
            "temperature_unit": group[0]["temperature_unit"],
        }
    return result


def process_runs(
    run_rows: list[dict[str, str]],
    calibration_rows: list[dict[str, str]],
) -> list[dict[str, object]]:
    run_report = validate_run_registry(run_rows)
    if not run_report["run_log_ready"]:
        raise ValueError(f"Run registry is not ready: {run_report['errors']}")
    calibrations = group_calibrations(calibration_rows)
    grouped_runs: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in run_rows:
        grouped_runs[row["run_id"]].append(row)

    outputs: list[dict[str, object]] = []
    for run_id, rows in sorted(grouped_runs.items()):
        rows.sort(key=lambda row: float(row["time"]))
        calibration_id = rows[0]["calibration_id"]
        if calibration_id not in calibrations:
            raise ValueError(f"{run_id}: missing calibration {calibration_id}")
        calibration = calibrations[calibration_id]
        concentration_unit = rows[0]["concentration_unit"]
        conductivity_unit = rows[0]["conductivity_unit"]
        volume_unit = rows[0]["volume_unit"]
        area_unit = rows[0]["area_unit"]
        time_unit = rows[0]["time_unit"]
        if concentration_unit != calibration["concentration_unit"]:
            raise ValueError(f"{run_id}: concentration unit differs from calibration")
        if conductivity_unit != calibration["conductivity_unit"]:
            raise ValueError(f"{run_id}: conductivity unit differs from calibration")
        for value, mapping, name in (
            (concentration_unit, CONCENTRATION_TO_MOL_M3, "concentration"),
            (volume_unit, VOLUME_TO_M3, "volume"),
            (area_unit, AREA_TO_M2, "area"),
            (time_unit, TIME_TO_S, "time"),
        ):
            if value not in mapping:
                raise ValueError(f"{run_id}: unsupported {name} unit {value!r}")

        concentration_factor = CONCENTRATION_TO_MOL_M3[concentration_unit]
        volume_factor = VOLUME_TO_M3[volume_unit]
        area_m2 = float(rows[0]["effective_area"]) * AREA_TO_M2[area_unit]
        time_factor = TIME_TO_S[time_unit]
        receiver = rows[0]["receiver_side"]
        removed = {"left": 0.0, "right": 0.0}
        added = {"left": 0.0, "right": 0.0}
        initial_total: float | None = None
        initial_receiver: float | None = None
        initial_time_s: float | None = None

        for row in rows:
            if row["temperature_unit"] != calibration["temperature_unit"]:
                raise ValueError(f"{run_id}: temperature unit differs from calibration")
            if float(row["temperature"]) != calibration["temperature"]:
                raise ValueError(f"{run_id}: temperature differs from calibration")

            concentration = {
                "left": interpolate_concentration(
                    float(row["left_conductivity"]), calibration["points"]
                ),
                "right": interpolate_concentration(
                    float(row["right_conductivity"]), calibration["points"]
                ),
            }
            amount = {
                side: concentration[side]
                * concentration_factor
                * float(row[f"{side}_volume"])
                * volume_factor
                for side in ("left", "right")
            }
            corrected_total = sum(amount.values()) + sum(removed.values()) - sum(added.values())
            corrected_receiver = amount[receiver] + removed[receiver] - added[receiver]
            time_s = float(row["time"]) * time_factor
            if initial_total is None:
                initial_total = corrected_total
                initial_receiver = corrected_receiver
                initial_time_s = time_s
            elapsed_s = time_s - initial_time_s
            transfer = corrected_receiver - initial_receiver
            outputs.append(
                {
                    "run_id": run_id,
                    "time": float(row["time"]),
                    "time_unit": time_unit,
                    "left_concentration": concentration["left"],
                    "right_concentration": concentration["right"],
                    "concentration_unit": concentration_unit,
                    "corrected_total_amount_mol": corrected_total,
                    "relative_closure_error": relative_closure_error(
                        corrected_total, initial_total
                    ),
                    "corrected_receiver_transfer_mol": transfer,
                    "cumulative_average_flux_mol_m2_s": (
                        None if elapsed_s == 0 else interval_flux(transfer, area_m2, elapsed_s)
                    ),
                    "input_level": row["input_level"],
                    "mechanistic_fit_authorized": False,
                }
            )

            for side in ("left", "right"):
                sampled_volume_m3 = float(row[f"{side}_sampled_volume"]) * volume_factor
                replenished_volume_m3 = (
                    float(row[f"{side}_replenished_volume"]) * volume_factor
                )
                replenishment_concentration = (
                    float(row[f"{side}_replenishment_concentration"])
                    * concentration_factor
                )
                removed[side] += concentration[side] * concentration_factor * sampled_volume_m3
                added[side] += replenishment_concentration * replenished_volume_m3
    return outputs


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert validated salt-transport logs and calculate corrected mass balance."
    )
    parser.add_argument("run_log", type=Path)
    parser.add_argument("calibration", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output = process_runs(load_csv(args.run_log), load_csv(args.calibration))
    text = json.dumps(output, ensure_ascii=False, indent=2)
    if args.output:
        args.output.write_text(text, encoding="utf-8")
    else:
        print(text)


if __name__ == "__main__":
    main()
