"""Build v9 figure, claim, run and inventory tables from saved evidence."""

import csv
import importlib.metadata
import json
import platform
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "WP6"
PACKAGE = ROOT / "原始资料/赵珈绮/18_薄化荷电主线任务书_v9.0.0"


def write_csv(name, rows):
    with (OUT / name).open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    package_relative = ("taskbook" if (ROOT / "taskbook").is_dir() else
                        "原始资料/赵珈绮/18_薄化荷电主线任务书_v9.0.0")
    source_manifest = f"{package_relative}/file_sha256.csv"
    figures = [
        {"figure_id": "v9_F1", "title": "冻结182校准与门槛",
         "svg": "data/processed/figures/v9_calibration_gate.svg",
         "png": "data/processed/figures/v9_calibration_gate.png",
         "source_csv": "WP4/v9_校准七档汇总.csv",
         "generation_script": "src/surface_pb/v9_plot_evidence.py",
         "role": "NO_COMPATIBLE_WITNESS; not material validation"},
        {"figure_id": "v9_F2", "title": "旧示例纯厚度与固定敏感性",
         "svg": "data/processed/figures/v9_model_only_thickness_sensitivity.svg",
         "png": "data/processed/figures/v9_model_only_thickness_sensitivity.png",
         "source_csv": "WP5/v9_MODEL_ONLY_精简厚度图.csv | WP5/v9_MODEL_ONLY_N4_N11数值误差.csv | WP5/v9_MODEL_ONLY_八项敏感性差值.csv",
         "generation_script": "src/surface_pb/v9_plot_evidence.py",
         "role": "MODEL_ONLY; not PE parameter fit"},
    ]
    for row in figures:
        if not all((ROOT / row[key]).is_file() for key in ("svg", "png")):
            raise FileNotFoundError(row["figure_id"])
    write_csv("v9_图源索引.csv", figures)

    claims = [
        ("C1", "输入为9组36Ramp且保留反向流量/低R²", "DATA_AUDIT_PASS", "WP1/v9_数据审计.json", "独立膜数、真实盐、仪器法向未认证"),
        ("C2", "四物种PB/荷电调节/观察算子在冻结案例技术通过", "TECHNICAL_PASS", "WP3/v9_技术验收.json", "技术PASS不证明真实PE化学"),
        ("C3", "182三组存在冻结相容化学", "NOT_SUPPORTED_OPTIMIZATION_UNRESOLVED", "WP4/v9_第4步判定.json", "最小E=1.0037188>1；8起点未证稳定"),
        ("C4", "共同化学预测357/604六组", "NOT_RUN_PREDECLARED_GATE", "WP4/v9_182校准与固定分支.md", "没有相容集合，不按留出集调参"),
        ("C5", "薄化增强真实PE界面荷电", "NOT_SUPPORTED", "WP5/v9_第5步验收.json", "旧示例N4/N11差低于10%和五倍误差，仅MODEL_ONLY"),
        ("C6", "有限PE介电敏感性触发预定物理风险界", "NO_TRIGGER_MODEL_ONLY", "WP5/v9_MODEL_ONLY_八项敏感性判定.json", "差值低于网格误差；不能据此测得真实ε或位点"),
        ("C7", "现有证据支持正面材料机制投稿", "NOT_SUPPORTED_MENTOR_REVIEW_PENDING", "WP6/v9_论文中文稿.md", "导师审阅未发生；理论边界稿也未获材料解释证据"),
    ]
    write_csv("v9_主张证据表.csv", [dict(claim_id=i, claim=c, status=s, evidence=e, boundary=b)
                                    for i, c, s, e, b in claims])

    def command(*paths):
        return "; ".join(".\\.venv\\Scripts\\python.exe " + path.replace("/", "\\") for path in paths)

    run_rows = [
        ("data", command("src/data_ingest/v9_data.py"), "0", "0", "WP1/v9_数据审计.json"),
        ("analytic_geometry", command("src/surface_pb/v9_analytic.py", "src/surface_pb/v9_geometry.py"), "0", "0", "WP2/v9_解析核对.json; WP3/v9_几何接口核对.json"),
        ("technical", command("src/surface_pb/v9_reference_checks.py", "src/surface_pb/v9_finite_checks.py", "src/surface_pb/v9_complete_technical.py"), "0", "technical_excluded", "WP3/v9_技术验收.json"),
        ("calibration", command("src/surface_pb/v9_calibrate.py", "src/surface_pb/v9_predict.py"), "0", "1668", "WP4/v9_校准目标调用.jsonl; WP4/v9_第4步判定.json"),
        ("model_only_maps", command("src/surface_pb/v9_model_only_map.py"), "0", "62", "WP5/v9_MODEL_ONLY_精简图判定.json"),
        ("final_numerics", command("src/surface_pb/v9_final_numerics.py", "src/surface_pb/v9_dielectric_validation.py"), "0", "25", "WP5/v9_MODEL_ONLY_十二状态数值复核.json; WP5/v9_二维介电求解器核对.json"),
        ("sensitivity", command("src/surface_pb/v9_sensitivity.py", "src/surface_pb/v9_sensitivity_analyze.py", "src/surface_pb/v9_step5_close.py"), "0", "126", "WP5/v9_第5步验收.json"),
        ("saved_result_export", command("src/surface_pb/v9_export_results.py"), "0", "0", "WP6/v9_平衡场索引.csv; WP6/v9_观察量索引.csv; WP6/v9_层荷电索引.csv"),
        ("paper", command("src/surface_pb/v9_plot_evidence.py", "src/surface_pb/v9_handoff_tables.py"), "0", "0", "WP6/v9_图源索引.csv; WP6/v9_主张证据表.csv"),
    ]
    write_csv("v9_运行清单.csv", [dict(stage=stage, actual_command=command,
                                    start_time="NOT_RECORDED", end_time="NOT_RECORDED",
                                    exit_code=exit_code, actual_PB_fields=fields,
                                    log_or_evidence=evidence, source_hash_manifest=source_manifest)
                                for stage, command, exit_code, fields, evidence in run_rows])
    manifest = {"execution_os": f"{platform.system()} {platform.release()}",
                "processor": platform.processor(),
                "python": platform.python_version(),
                "numpy": importlib.metadata.version("numpy"),
                "scipy": importlib.metadata.version("scipy"),
                "matplotlib": importlib.metadata.version("matplotlib"),
                "per_stage_start_end_times": "NOT_RECORDED_AT_EXECUTION; do not infer from file times",
                "project_branch_at_handoff": "dev-zhao",
                "frozen_input_hashes_existing_manifest": source_manifest,
                "original_package_source_manifest": f"{package_relative}/source_manifest.csv",
                "original_package_verifier": f"{package_relative}/preparation/verify_inputs.py",
                "new_file_content_hashes_computed": False,
                "geometry_hashes_computed_as_required_by_taskbook": True,
                "calibration_actual_objective_calls": 556,
                "calibration_PB_fields": 1668,
                "model_only_map_and_control_PB_fields": 62,
                "final_numerical_additional_PB_fields": 25,
                "sensitivity_PB_fields": 126,
                "technical_cases_excluded_from_research_PB_budget": True,
                "mentor_review": "NOT_OCCURRED"}
    (OUT / "v9_运行记录.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
                                           encoding="utf-8")
    print("v9 handoff tables: 2 figures, 7 claims, 9 run stages")


if __name__ == "__main__":
    main()
