from __future__ import annotations

import json

from conductivity_calibration import (
    TEMPLATE_PATH as CALIBRATION_TEMPLATE_PATH,
    load_registry as load_calibration_registry,
    validate_registry as validate_calibration_registry,
)
from validate_blank_registry import (
    TEMPLATE_PATH as BLANK_TEMPLATE_PATH,
    load_registry as load_blank_registry,
    validate_registry as validate_blank_registry,
)
from validate_experiment_records import load_tracker, validate_tracker
from validate_independent_readouts import (
    TEMPLATE_PATH as READOUT_TEMPLATE_PATH,
    load_registry as load_readout_registry,
    validate_registry as validate_readout_registry,
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


def assess() -> dict[str, object]:
    salt = validate_tracker(load_tracker())
    osw = validate_osw_registry(load_osw_registry(OSW_TEMPLATE_PATH))
    readouts = validate_readout_registry(
        load_readout_registry(READOUT_TEMPLATE_PATH)
    )
    calibration = validate_calibration_registry(
        load_calibration_registry(CALIBRATION_TEMPLATE_PATH)
    )
    runs = validate_run_registry(load_run_registry(RUN_TEMPLATE_PATH))
    blanks = validate_blank_registry(load_blank_registry(BLANK_TEMPLATE_PATH))
    thresholds = validate_threshold_registry(
        load_threshold_registry(THRESHOLD_TEMPLATE_PATH)
    )

    salt_metadata_complete = (
        salt["valid"]
        and salt["unresolved_curve_count"] == 0
        and salt["status_counts"].get("complete", 0) == salt["curve_count"]
    )
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
    if not osw["protocol_ready"]:
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
        "osw_protocol": {
            "registry_valid": osw["valid"],
            "row_count": osw["row_count"],
            "protocol_ready": osw["protocol_ready"],
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
    print(json.dumps(assess(), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
