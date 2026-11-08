"""Record the approved 19-package representative replay and bounded inventory check."""

import csv
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEST = ROOT / "原始资料/赵珈绮/19_学生薄化荷电成果_v9.0.0"
MANIFEST = ROOT / "WP6/v9_拟归档清单.csv"


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def main():
    records = read_csv(MANIFEST)
    missing = [r["archive_relative"] for r in records
               if not (DEST / r["archive_relative"]).is_file()]
    if missing:
        raise FileNotFoundError(f"missing {len(missing)} archived files: {missing[:3]}")
    technical = read_json(DEST / "WP3/v9_技术验收.json")
    calibration = read_json(DEST / "WP4/v9_第4步判定.json")
    step5 = read_json(DEST / "WP5/v9_第5步验收.json")
    fields = read_csv(DEST / "WP6/v9_平衡场索引.csv")
    observed = read_csv(DEST / "WP6/v9_观察量索引.csv")
    layers = read_csv(DEST / "WP6/v9_层荷电索引.csv")
    claims = read_csv(DEST / "WP6/v9_主张证据表.csv")
    runs = read_csv(DEST / "WP6/v9_运行清单.csv")
    images = ["v9_calibration_gate.svg", "v9_calibration_gate.png",
              "v9_model_only_thickness_sensitivity.svg",
              "v9_model_only_thickness_sensitivity.png"]
    if not all((DEST / "data/processed/figures" / name).is_file() for name in images):
        raise FileNotFoundError("missing representative figure after replay")
    if technical["overall"] != "TECHNICAL_PASS" or technical["passed_checks"] != 10:
        raise AssertionError("archive technical replay failed")
    if calibration["compatible_gamma_count"] != 0 or calibration["actual_objective_calls"] != 556:
        raise AssertionError("archive calibration evidence changed")
    if step5["execution_status"] != "COMPLETE_PREDECLARED_MODEL_ONLY_BRANCH":
        raise AssertionError("archive Step 5 replay failed")
    if (len(fields), len(observed), len(layers), len(claims), len(runs)) != (202, 110, 1598, 7, 9):
        raise AssertionError("archive saved-result indexes incomplete")
    record = {
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "archive_relative": "原始资料/赵珈绮/19_学生薄化荷电成果_v9.0.0",
        "manifest_file_count": len(records),
        "manifest_targets_present": len(records),
        "source_verifier_command": "python taskbook/preparation/verify_inputs.py",
        "source_verifier_observed": {"source_copies": 23, "package_files_checked": 54,
                                     "ramps": 36, "groups": 9, "status": "PASS", "exit_code": 0},
        "representative_replay": [
            {"command": "python src/surface_pb/v9_complete_technical.py",
             "exit_code": 0, "observed": "TECHNICAL_PASS 10/10"},
            {"command": "python src/surface_pb/v9_step5_close.py",
             "exit_code": 0, "observed": "COMPLETE_PREDECLARED_MODEL_ONLY_BRANCH"},
            {"command": "python src/surface_pb/v9_export_results.py",
             "exit_code": 0, "observed": "202/110/1598 saved-result rows"},
            {"command": "python src/surface_pb/v9_plot_evidence.py",
             "exit_code": 0, "observed": "2 SVG and 2 PNG regenerated"},
            {"command": "python src/surface_pb/v9_handoff_tables.py",
             "exit_code": 0, "observed": "2 figures, 7 claims, 9 run stages"},
        ],
        "calibration_repeated_in_archive": False,
        "archive_replay_status": "PASS_REPRESENTATIVE_REPLAY",
        "scientific_material_support": False,
        "mentor_review": "NOT_OCCURRED",
    }
    for path in (ROOT / "WP6/v9_包内回放记录.json", DEST / "WP6/v9_包内回放记录.json"):
        path.write_text(json.dumps(record, ensure_ascii=False, indent=2) + "\n",
                        encoding="utf-8")
    print(json.dumps({"archive_replay_status": record["archive_replay_status"],
                      "manifest_targets_present": len(records)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
