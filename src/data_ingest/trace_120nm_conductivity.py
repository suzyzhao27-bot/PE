"""Trace the Book1 120 nm conductivity curve without filling missing experiment records."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "data/processed/origin_exports/Polymer_salt_2h/Book1__Sheet1.csv"
METADATA = SOURCE.with_name("Book1__Sheet1.metadata.json")
OUTPUT = ROOT / "data/processed/salt_transport/book1_120nm_calculation_chain.json"


def main() -> None:
    with SOURCE.open(encoding="utf-8-sig", newline="") as stream:
        source_rows = list(csv.DictReader(stream))
    metadata = json.loads(METADATA.read_text(encoding="utf-8-sig"))
    time_column, conductivity_column = metadata["columns"][2:4]
    if ((time_column["short_name"], time_column["long_name"], time_column["units"])
            != ("C", "Time", "min")):
        raise ValueError("Book1 time-column evidence changed")
    if ((conductivity_column["short_name"], conductivity_column["long_name"],
         conductivity_column["units"], conductivity_column["comments"])
            != ("D", "Ionic conductivity", "μS/cm", "120 nm")):
        raise ValueError("Book1 conductivity-column evidence changed")

    points = [
        {"time_min": float(row["col_3_C"]),
         "conductivity_uS_per_cm": float(row["col_4_D"])}
        for row in source_rows if row["col_3_C"] and row["col_4_D"]
    ]
    if len(points) != 7 or any(right["time_min"] <= left["time_min"]
                               for left, right in zip(points, points[1:])):
        raise ValueError("Unexpected Book1 120 nm time series")
    first, last = points[0], points[-1]
    if first["time_min"] != 0 or last["time_min"] != 120:
        raise ValueError("Unexpected 120-minute endpoints")
    delta = round(last["conductivity_uS_per_cm"] - first["conductivity_uS_per_cm"], 12)
    result = {
        "curve_id": "Polymer_salt_2h:Book1:Sheet1:120 nm",
        "source_csv": SOURCE.relative_to(ROOT).as_posix(),
        "source_metadata": METADATA.relative_to(ROOT).as_posix(),
        "x_column": "C",
        "y_column": "D",
        "sample_label": "120 nm (Origin column comment; independent membrane ID unknown)",
        "points": points,
        "calculated_from_source": {
            "delta_time_min": last["time_min"] - first["time_min"],
            "delta_time_s": (last["time_min"] - first["time_min"]) * 60,
            "delta_conductivity_uS_per_cm": delta,
            "delta_conductivity_S_per_m": round(delta * 1e-4, 12),
            "endpoint_average_delta_conductivity_uS_per_cm_per_min":
                round(delta / (last["time_min"] - first["time_min"]), 12),
        },
        "stages": {
            "conductivity": {
                "status": "source_explicit_units; numeric delta available",
                "meaning": "receiver-side conductivity observation, not concentration",
            },
            "concentration": {
                "status": "unknown",
                "needed": ["curve-specific salt identity or preparation record",
                           "conductivity-to-concentration calibration with validity range",
                           "measurement temperature and compensation settings"],
                "calculation_when_available": "c_R(t) = calibrated_inverse(kappa_R(t), salt, T)",
            },
            "flux": {
                "status": "unknown",
                "needed_after_concentration": ["receiver volume history",
                                               "sampling and replenishment log",
                                               "effective exposed membrane area"],
                "calculation_when_available":
                    "J_interval = receiver NaCl amount balance / (effective area * interval seconds)",
            },
            "cross_membrane_or_charge_mechanism": {
                "status": "unknown",
                "needed_after_flux": ["wetting history", "mixing conditions",
                                      "parent membrane and coupon identity",
                                      "independent replicate mapping"],
            },
        },
        "related_plan_not_curve_confirmation": {
            "source": "原始资料/赵珈绮/01_研究框架与资料_v2/05_真实数据_v1核查/Polymer_salt_experiment_plan.txt",
            "salt_named": "NaCl",
            "recipe_written": "8.75 g made up to 250 mL = 35 g/L mass per solution volume",
            "warning": "Plan does not confirm this curve followed it; 35 g/L is not a verified 3.5 wt%.",
        },
    }
    OUTPUT.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"points={len(points)} delta_kappa={delta:.2f} μS/cm concentration=unknown flux=unknown")


if __name__ == "__main__":
    main()
