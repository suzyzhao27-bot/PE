from __future__ import annotations

import argparse
import json
from pathlib import Path

from conductivity_calibration import (
    TEMPLATE_PATH as CALIBRATION_TEMPLATE_PATH,
    load_registry as load_calibration_registry,
    validate_registry as validate_calibration_registry,
)
from evaluate_liquid_compatibility import evaluate_records as evaluate_compatibility_records
from validate_blank_registry import (
    TEMPLATE_PATH as BLANK_TEMPLATE_PATH,
    load_registry as load_blank_registry,
    validate_registry as validate_blank_registry,
)
from validate_experiment_records import TRACKER_PATH, load_tracker, validate_tracker
from validate_independent_readouts import (
    TEMPLATE_PATH as READOUT_TEMPLATE_PATH,
    load_registry as load_readout_registry,
    validate_registry as validate_readout_registry,
)
from validate_liquid_compatibility import (
    TEMPLATE_PATH as COMPATIBILITY_TEMPLATE_PATH,
    load_registry as load_compatibility_registry,
    validate_registry as validate_compatibility_registry,
)
from validate_liquid_compatibility_thresholds import (
    TEMPLATE_PATH as COMPATIBILITY_THRESHOLD_TEMPLATE_PATH,
    load_registry as load_compatibility_threshold_registry,
    validate_registry as validate_compatibility_threshold_registry,
)
from validate_osw_protocol import (
    TEMPLATE_PATH as OSW_TEMPLATE_PATH,
    load_registry as load_osw_registry,
    validate_registry as validate_osw_registry,
)
from validate_salt_transport_runs import (
    TEMPLATE_PATH as RUN_TEMPLATE_PATH,
    load_registry as load_run_registry,
    validate_registry as validate_run_registry,
)
from validate_quality_thresholds import (
    TEMPLATE_PATH as THRESHOLD_TEMPLATE_PATH,
    load_registry as load_threshold_registry,
    validate_registry as validate_threshold_registry,
)


def compatibility_links_ready(
    compatibility_rows: list[dict[str, str]],
    ready_compatible_ids: list[str],
    osw_rows: list[dict[str, str]],
) -> bool:
    ready_ids = set(ready_compatible_ids)
    compatible_liquids = {
        row["compatibility_record_id"].strip(): (
            row["protocol_liquid_identity"].strip(),
            row["protocol_liquid_batch"].strip(),
        )
        for row in compatibility_rows
        if row["compatibility_record_id"].strip() in ready_ids
    }
    return bool(osw_rows) and all(
        row["status"].strip() == "excluded"
        or compatible_liquids.get(row["compatibility_record_id"].strip())
        == (
            row["protocol_liquid_identity"].strip(),
            row["protocol_liquid_batch"].strip(),
        )
        for row in osw_rows
    )


def assess(
    *,
    salt_tracker_path: Path = TRACKER_PATH,
    compatibility_registry_path: Path = COMPATIBILITY_TEMPLATE_PATH,
    compatibility_threshold_path: Path = COMPATIBILITY_THRESHOLD_TEMPLATE_PATH,
    osw_registry_path: Path = OSW_TEMPLATE_PATH,
    calibration_path: Path = CALIBRATION_TEMPLATE_PATH,
    run_log_path: Path = RUN_TEMPLATE_PATH,
    blank_registry_path: Path = BLANK_TEMPLATE_PATH,
    threshold_registry_path: Path = THRESHOLD_TEMPLATE_PATH,
    readout_registry_path: Path = READOUT_TEMPLATE_PATH,
) -> dict[str, object]:
    salt = validate_tracker(load_tracker(salt_tracker_path))
    compatibility_rows = load_compatibility_registry(compatibility_registry_path)
    compatibility = validate_compatibility_registry(compatibility_rows)
    compatibility_threshold_rows = load_compatibility_threshold_registry(
        compatibility_threshold_path
    )
    compatibility_thresholds = validate_compatibility_threshold_registry(
        compatibility_threshold_rows
    )
    compatibility_evaluation = evaluate_compatibility_records(
        compatibility_rows, compatibility_threshold_rows
    )
    osw_rows = load_osw_registry(osw_registry_path)
    osw = validate_osw_registry(osw_rows)
    readouts = validate_readout_registry(
        load_readout_registry(readout_registry_path)
    )
    calibration = validate_calibration_registry(
        load_calibration_registry(calibration_path)
    )
    runs = validate_run_registry(load_run_registry(run_log_path))
    blanks = validate_blank_registry(load_blank_registry(blank_registry_path))
    thresholds = validate_threshold_registry(
        load_threshold_registry(threshold_registry_path)
    )

    salt_metadata_complete = (
        salt["valid"]
        and salt["unresolved_curve_count"] == 0
        and salt["status_counts"].get("complete", 0) == salt["curve_count"]
    )
    osw_compatibility_links_ready = compatibility_links_ready(
        compatibility_rows,
        compatibility_evaluation["evaluated_ready_compatible_ids"],
        osw_rows,
    )
    osw_protocol_ready = osw["protocol_ready"] and osw_compatibility_links_ready
    next_actions: list[str] = []
    if not salt_metadata_complete:
        next_actions.append("complete curve-level salt transport records")
    if not calibration["conversion_ready"]:
        next_actions.append("provide a valid conductivity calibration table")
    if not runs["run_log_ready"]:
        next_actions.append("populate and review salt-transport run logs")
    if not blanks["salt_blank_set_ready"]:
        next_actions.append("qualify no-membrane and impermeable salt-transport blanks")
    if not thresholds["threshold_set_ready"]:
        next_actions.append("predeclare and approve quality thresholds")
    if not compatibility_thresholds["thresholds_ready"]:
        next_actions.append("predeclare protocol-liquid compatibility thresholds")
    if not compatibility_evaluation["compatibility_ready"]:
        next_actions.append("qualify protocol-liquid compatibility")
    if not osw_protocol_ready:
        next_actions.append("populate and review the O/S/W sample registry")
    if not readouts["selection_ready"]:
        next_actions.append("select and qualify Na, Cl, and neutral-tracer readouts")

    return {
        "salt_transport": {
            "tracker_valid": salt["valid"],
            "curve_count": salt["curve_count"],
            "unresolved_curve_count": salt["unresolved_curve_count"],
            "metadata_complete": salt_metadata_complete,
        },
        "liquid_compatibility": {
            "registry_valid": compatibility["valid"],
            "row_count": compatibility["row_count"],
            "evaluation_valid": compatibility_evaluation["valid"],
            "evaluated_ready_compatible_count": len(
                compatibility_evaluation["evaluated_ready_compatible_ids"]
            ),
            "compatibility_ready": compatibility_evaluation["compatibility_ready"],
        },
        "liquid_compatibility_thresholds": {
            "registry_valid": compatibility_thresholds["valid"],
            "row_count": compatibility_thresholds["row_count"],
            "thresholds_ready": compatibility_thresholds["thresholds_ready"],
        },
        "osw_protocol": {
            "registry_valid": osw["valid"],
            "row_count": osw["row_count"],
            "compatibility_links_ready": osw_compatibility_links_ready,
            "protocol_ready": osw_protocol_ready,
        },
        "conductivity_calibration": {
            "registry_valid": calibration["valid"],
            "row_count": calibration["row_count"],
            "conversion_ready": calibration["conversion_ready"],
        },
        "salt_transport_runs": {
            "registry_valid": runs["valid"],
            "row_count": runs["row_count"],
            "run_log_ready": runs["run_log_ready"],
        },
        "blanks": {
            "registry_valid": blanks["valid"],
            "row_count": blanks["row_count"],
            "salt_blank_set_ready": blanks["salt_blank_set_ready"],
            "potential_blank_set_ready": blanks["potential_blank_set_ready"],
        },
        "quality_thresholds": {
            "registry_valid": thresholds["valid"],
            "row_count": thresholds["row_count"],
            "threshold_set_ready": thresholds["threshold_set_ready"],
        },
        "independent_readouts": {
            "registry_valid": readouts["valid"],
            "row_count": readouts["row_count"],
            "selection_ready": readouts["selection_ready"],
        },
        "next_actions": next_actions,
        "mechanistic_fit_authorized": False,
        "experiment_authorized": False,
        "authorization_note": (
            "This report summarizes record readiness only. Formal scientific and signed "
            "experimental review remain separate gates."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Assess project input readiness.")
    parser.add_argument("--salt-tracker", type=Path, default=TRACKER_PATH)
    parser.add_argument(
        "--compatibility-registry", type=Path, default=COMPATIBILITY_TEMPLATE_PATH
    )
    parser.add_argument(
        "--compatibility-thresholds",
        type=Path,
        default=COMPATIBILITY_THRESHOLD_TEMPLATE_PATH,
    )
    parser.add_argument("--osw-registry", type=Path, default=OSW_TEMPLATE_PATH)
    parser.add_argument("--calibration", type=Path, default=CALIBRATION_TEMPLATE_PATH)
    parser.add_argument("--run-log", type=Path, default=RUN_TEMPLATE_PATH)
    parser.add_argument("--blank-registry", type=Path, default=BLANK_TEMPLATE_PATH)
    parser.add_argument("--threshold-registry", type=Path, default=THRESHOLD_TEMPLATE_PATH)
    parser.add_argument("--readout-registry", type=Path, default=READOUT_TEMPLATE_PATH)
    args = parser.parse_args()
    report = assess(
        salt_tracker_path=args.salt_tracker,
        compatibility_registry_path=args.compatibility_registry,
        compatibility_threshold_path=args.compatibility_thresholds,
        osw_registry_path=args.osw_registry,
        calibration_path=args.calibration,
        run_log_path=args.run_log,
        blank_registry_path=args.blank_registry,
        threshold_registry_path=args.threshold_registry,
        readout_registry_path=args.readout_registry,
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
