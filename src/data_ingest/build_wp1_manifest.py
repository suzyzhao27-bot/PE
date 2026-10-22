"""Build the v6 WP1 source and parameter inventory from the audited tidy table."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "data/processed/salt_transport/salt_transport_tidy.csv"
TARGET = ROOT / "WP1/03_数据与参数清单.csv"
RAW_DIR = "原始资料/赵珈绮/01_研究框架与资料_v2/05_真实数据_v1核查/source_snapshots"

FIELDS = [
    "record_type", "dataset_id", "source_file", "sheet", "time_column",
    "value_column", "thickness_nm", "time_min", "conductivity_raw",
    "conductivity_unit_status", "replicate_group", "use_role", "identity_note",
    "parameter_symbol", "parameter_role", "parameter_unit", "parameter_rule",
    "source_ref",
]

BOOK1 = {
    "65 nm": ("A", "B"), "120 nm": ("C", "D"),
    "200 nm": ("C", "E"), "350 nm": ("C", "F"),
    "700 nm": ("C", "G"), "2600 nm": ("C", "H"),
    "PC filter, 75 μm thick, 200 nm pore": ("C", "I"),
    "filter paper, 100 μm thick, 200 nm pore": ("J", "K"),
}
BOOK2 = {
    "R65 nm": ("A", "B"), "R700 nm": ("A", "C"),
    "R2600 nm": ("A", "D"), "65 nm": ("E", "F"),
    "700 nm": ("G", "H"), "2600 nm": ("G", "I"),
}
DRY = {"700 nm": ("A", "B"), "350 nm": ("C", "D"), "200 nm": ("E", "F")}

PARAMETERS = [
    ("h", "input_geometry_label", "m", "Book1 65/120/200/350/700/2600 nm; convert nm to m"),
    ("a", "fixed_geometry_scenario", "m", "main 25 nm; common 20/30 nm sensitivity"),
    ("r", "fixed_boundary_scenario", "1", "main 0; sensitivity 0.75"),
    ("r_V", "fixed_reservoir_scenario", "1", "main 1; sensitivity 0.5/2"),
    ("D_b", "fixed_reference_not_measured", "m2/s", "2 D_Na D_Cl/(D_Na+D_Cl); D_Na=1.33e-9; D_Cl=2.03e-9"),
    ("w_s", "shared_fit_M0_M1_M2", "kBT", "diagnostic search 0..20"),
    ("ell", "shared_fit_M0_M1_M2", "m", "diagnostic search 10..10000 nm"),
    ("D_m", "shared_fit_M0_M1_M2", "m2/s", "D_m/D_b diagnostic search 1e-12..1"),
    ("gamma=A epsilon/V_R", "shared_fit_M0_M1_M2", "1/m", "diagnostic search 1e-2..1e6"),
    ("Q", "shared_observation_fit_M0_M1_M2", "original conductivity unit", "diagnostic search 1..1e8; no per-film Q"),
    ("chi", "shared_fit_M2_only", "1", "independent zero baseline; positive diagnostic search 1e-4..1e6"),
    ("tau_x", "shared_fit_M2_only", "s", "diagnostic search 1e-3..1e6"),
    ("epsilon", "not_individually_identified", "1", "only within gamma and P_eff unless A/V_R known"),
    ("A/V_R", "not_individually_identified", "1/m", "unknown; capacity constraint A/V_R >= gamma(1+chi) for dead-end water"),
    ("C0", "reference_concentration_scale", "mol/m3", "unknown absolute calibration; solve dimensionless c/C0"),
    ("conductivity_calibration", "unknown_metadata", "unknown", "not recovered for Book1; no absolute flux inference"),
    ("temperature", "unknown_metadata", "unknown", "not recovered for Book1; D_b is reference scenario"),
    ("sampling_and_wetting_history", "unknown_metadata", "unknown", "not recovered for Book1"),
]


def classify(row: dict[str, str]) -> tuple[str, str, str]:
    book = row["workbook"]
    sample = row["sample_id"]
    source = row["source_opju"]
    if source == "Polymer_salt_dry_6h.opju":
        return "auxiliary_dry_condition", f"dry_6h_{sample}", "different condition; time/conductivity units need record confirmation"
    if book == "Book2" and sample in {"65 nm", "700 nm", "2600 nm"}:
        note = "same time and values" if sample == "65 nm" else "same conductivity sequence; one time entry differs"
        return "provenance_duplicate", f"book1_{sample}", note + "; not an independent replicate"
    if book == "Book2":
        return "auxiliary_reverse_unpaired", f"reverse_{sample}", "R label; membrane/sample pairing not established"
    if sample in {"PC filter, 75 μm thick, 200 nm pore", "filter paper, 100 μm thick, 200 nm pore"}:
        return "auxiliary_material_control", sample, "different material; do not share PE field or fit"
    return "primary_conditional_fit", f"book1_{sample}", "one thickness-labelled curve; independent membrane identity unconfirmed"


def main() -> None:
    with SOURCE.open(encoding="utf-8-sig", newline="") as stream:
        source_rows = list(csv.DictReader(stream))
    assert len(source_rows) == 167
    counts = Counter(row["curve_id"] for row in source_rows)
    assert len(counts) == 17
    output: list[dict[str, str]] = []
    for row in source_rows:
        sample = row["sample_id"]
        book = row["workbook"]
        if row["source_opju"] == "Polymer_salt_dry_6h.opju":
            time_col, value_col = DRY[sample]
        elif book == "Book1":
            time_col, value_col = BOOK1[sample]
        else:
            time_col, value_col = BOOK2[sample]
        role, group, note = classify(row)
        thickness = sample.removesuffix(" nm") if sample.endswith(" nm") and sample.removesuffix(" nm").isdigit() else ""
        output.append({
            "record_type": "observation",
            "dataset_id": row["curve_id"],
            "source_file": f"{RAW_DIR}/{row['source_opju']}",
            "sheet": f"{book}/{row['worksheet']}",
            "time_column": time_col,
            "value_column": value_col,
            "thickness_nm": thickness,
            "time_min": row["time"],
            "conductivity_raw": row["conductivity"],
            "conductivity_unit_status": row["conductivity_unit"] + "; " + row["unit_evidence"],
            "replicate_group": group,
            "use_role": role,
            "identity_note": note,
            "source_ref": "data/processed/salt_transport/salt_transport_tidy.csv; data/processed/origin_exports/Polymer_salt_2h or Polymer_salt_dry_6h metadata",
        })
    assert sum(r["use_role"] == "primary_conditional_fit" for r in output) == 42
    assert sum(r["use_role"] == "provenance_duplicate" for r in output) == 21
    for symbol, role, unit, rule in PARAMETERS:
        output.append({
            "record_type": "parameter_or_unknown", "parameter_symbol": symbol,
            "parameter_role": role, "parameter_unit": unit,
            "parameter_rule": rule,
            "source_ref": "原始资料/赵珈绮/12_论文导向研究总方案_v6.0.0/01_物理模型与推导规格.md;04_计算矩阵与验收规则.md",
        })
    TARGET.parent.mkdir(exist_ok=True)
    with TARGET.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    print(f"{TARGET}: {len(source_rows)} observations in {len(counts)} curves; {len(PARAMETERS)} parameter rows")


if __name__ == "__main__":
    main()
