from __future__ import annotations

import csv
import json
from dataclasses import dataclass, asdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
EXPORT_ROOT = ROOT / "data" / "processed" / "origin_exports"
OUTPUT_ROOT = ROOT / "data" / "processed" / "salt_transport"


@dataclass(frozen=True)
class CurveMap:
    project: str
    workbook: str
    worksheet: str
    sample_id: str
    condition: str
    time_column: str
    conductivity_column: str
    time_unit: str
    conductivity_unit: str
    sample_id_evidence: str
    time_column_evidence: str
    unit_evidence: str


CURVES = [
    CurveMap("Polymer_salt_2h", "Book1", "Sheet1", "65 nm", "2h", "col_1_A", "col_2_B", "min", "μS/cm", "Origin comment on conductivity column", "paired explicit Time/Ionic conductivity columns", "explicit Origin units"),
    CurveMap("Polymer_salt_2h", "Book1", "Sheet1", "120 nm", "2h", "col_3_C", "col_4_D", "min", "μS/cm", "Origin comment on conductivity column", "paired explicit Time/Ionic conductivity columns", "explicit Origin units"),
    CurveMap("Polymer_salt_2h", "Book1", "Sheet1", "200 nm", "2h", "col_3_C", "col_5_E", "min", "μS/cm", "Origin comment on conductivity column", "saved Origin graph X/Y binding", "units inherited from explicit columns in same sheet"),
    CurveMap("Polymer_salt_2h", "Book1", "Sheet1", "350 nm", "2h", "col_3_C", "col_6_F", "min", "μS/cm", "Origin comment on conductivity column", "saved Origin graph X/Y binding", "units inherited from explicit columns in same sheet"),
    CurveMap("Polymer_salt_2h", "Book1", "Sheet1", "700 nm", "2h", "col_3_C", "col_7_G", "min", "μS/cm", "Origin comment on conductivity column", "saved Origin graph X/Y binding", "units inherited from explicit columns in same sheet"),
    CurveMap("Polymer_salt_2h", "Book1", "Sheet1", "2600 nm", "2h", "col_3_C", "col_8_H", "min", "μS/cm", "Origin comment on conductivity column", "saved Origin graph X/Y binding", "units inherited from explicit columns in same sheet"),
    CurveMap("Polymer_salt_2h", "Book1", "Sheet1", "PC filter, 75 μm thick, 200 nm pore", "2h", "col_3_C", "col_9_I", "min", "μS/cm", "Origin comment on conductivity column", "saved Origin graph X/Y binding", "units inherited from explicit columns in same sheet"),
    CurveMap("Polymer_salt_2h", "Book1", "Sheet1", "filter paper, 100 μm thick, 200 nm pore", "2h", "col_10_J", "col_11_K", "min", "μS/cm", "Origin comment on conductivity column", "adjacent time/conductivity pair inferred from table layout", "units inherited from explicit columns in same sheet"),
    CurveMap("Polymer_salt_2h", "Book2", "Sheet1", "R65 nm", "2h", "col_1_A", "col_2_B", "min", "μS/cm", "Origin comment on conductivity column", "paired explicit Time/Ionic conductivity columns", "explicit Origin units"),
    CurveMap("Polymer_salt_2h", "Book2", "Sheet1", "R700 nm", "2h", "col_1_A", "col_3_C", "min", "μS/cm", "Origin comment on conductivity column", "shared time column inferred from table layout", "units inherited from explicit columns in same sheet"),
    CurveMap("Polymer_salt_2h", "Book2", "Sheet1", "R2600 nm", "2h", "col_1_A", "col_4_D", "min", "μS/cm", "Origin comment on conductivity column", "shared time column inferred from table layout", "units inherited from explicit columns in same sheet"),
    CurveMap("Polymer_salt_2h", "Book2", "Sheet1", "65 nm", "2h", "col_5_E", "col_6_F", "min", "μS/cm", "Origin comment on conductivity column", "adjacent time/conductivity pair inferred from table layout", "units inherited from explicit columns in same sheet"),
    CurveMap("Polymer_salt_2h", "Book2", "Sheet1", "700 nm", "2h", "col_7_G", "col_8_H", "min", "μS/cm", "Origin comment on conductivity column", "adjacent time/conductivity pair inferred from table layout", "units inherited from explicit columns in same sheet"),
    CurveMap("Polymer_salt_2h", "Book2", "Sheet1", "2600 nm", "2h", "col_7_G", "col_9_I", "min", "μS/cm", "Origin comment on conductivity column", "shared time column inferred from table layout", "units inherited from explicit columns in same sheet"),
    CurveMap("Polymer_salt_dry_6h", "Book1", "未润湿体系数据", "700 nm", "dry_6h", "col_1_A", "col_2_B", "min", "μS/cm", "Origin comment on conductivity column", "adjacent paired columns inferred from raw-sheet layout", "units inferred from companion salt project; requires record confirmation"),
    CurveMap("Polymer_salt_dry_6h", "Book1", "未润湿体系数据", "350 nm", "dry_6h", "col_3_C", "col_4_D", "min", "μS/cm", "Origin comment on conductivity column", "adjacent paired columns inferred from raw-sheet layout", "units inferred from companion salt project; requires record confirmation"),
    CurveMap("Polymer_salt_dry_6h", "Book1", "未润湿体系数据", "200 nm", "dry_6h", "col_5_E", "col_6_F", "min", "μS/cm", "Origin comment on conductivity column", "adjacent paired columns inferred from raw-sheet layout", "units inferred from companion salt project; requires record confirmation"),
]


FIELDS = [
    "source_opju",
    "workbook",
    "worksheet",
    "curve_id",
    "sample_id",
    "condition",
    "point_index",
    "time",
    "time_unit",
    "conductivity",
    "conductivity_unit",
    "sample_id_evidence",
    "time_column_evidence",
    "unit_evidence",
    "eligible_for_mechanistic_fit",
    "exclusion_reason",
]


def load_rows(curve: CurveMap) -> list[dict[str, str]]:
    path = EXPORT_ROOT / curve.project / f"{curve.workbook}__{curve.worksheet}.csv"
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows or curve.time_column not in rows[0] or curve.conductivity_column not in rows[0]:
        raise ValueError(f"Mapping does not match exported table: {curve}")
    return rows


def normalize_curve(curve: CurveMap) -> list[dict[str, object]]:
    result: list[dict[str, object]] = []
    curve_id = f"{curve.project}:{curve.workbook}:{curve.worksheet}:{curve.sample_id}"
    for index, raw in enumerate(load_rows(curve), start=1):
        time_value = raw[curve.time_column].strip()
        conductivity_value = raw[curve.conductivity_column].strip()
        if not time_value and not conductivity_value:
            continue
        if not time_value or not conductivity_value:
            continue
        result.append(
            {
                "source_opju": f"{curve.project}.opju",
                "workbook": curve.workbook,
                "worksheet": curve.worksheet,
                "curve_id": curve_id,
                "sample_id": curve.sample_id,
                "condition": curve.condition,
                "point_index": index,
                "time": float(time_value),
                "time_unit": curve.time_unit,
                "conductivity": float(conductivity_value),
                "conductivity_unit": curve.conductivity_unit,
                "sample_id_evidence": curve.sample_id_evidence,
                "time_column_evidence": curve.time_column_evidence,
                "unit_evidence": curve.unit_evidence,
                "eligible_for_mechanistic_fit": False,
                "exclusion_reason": "conductivity calibration, temperature, effective area, bath volumes, sampling/replenishment, wetting history, and independent replicate identity are not all closed",
            }
        )
    return result


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    all_rows: list[dict[str, object]] = []
    for curve in CURVES:
        all_rows.extend(normalize_curve(curve))

    two_hour = [row for row in all_rows if row["condition"] == "2h"]
    dry = [row for row in all_rows if row["condition"] == "dry_6h"]
    write_csv(OUTPUT_ROOT / "salt_2h_tidy.csv", two_hour)
    write_csv(OUTPUT_ROOT / "salt_dry_6h_tidy.csv", dry)
    write_csv(OUTPUT_ROOT / "salt_transport_tidy.csv", all_rows)

    audit = {
        "curve_count": len(CURVES),
        "point_count": len(all_rows),
        "curve_mappings": [asdict(curve) for curve in CURVES],
        "fit_eligibility": False,
        "blocking_fields": [
            "conductivity-to-concentration calibration and temperature compensation",
            "effective membrane area and sealing geometry",
            "initial and time-varying bath volumes",
            "sampling/replenishment and mixing log",
            "curve-level wetting history",
            "independent membrane/parent-film replicate mapping",
            "confirmation of inherited units and time-column mappings not fixed by saved graphs",
        ],
        "excluded_origin_sheets": "All FitLinear* and FitPolynomial* sheets are derivative fits, not raw observations.",
    }
    (OUTPUT_ROOT / "normalization_audit.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    print(json.dumps({"curves": len(CURVES), "points": len(all_rows)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
