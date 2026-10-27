"""Prepare the user-approved v8.0.1 formal archive; do not run before consent."""

import csv
import json
import shutil
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
TASKBOOK=ROOT/"原始资料/赵珈绮/16_物理诊断与接续任务_v8.0.1"
TARGET=ROOT/"17_学生物理约束接续成果_v8.0.1"


def copy_file(source,target):
    if not source.is_file():
        raise FileNotFoundError(source)
    target.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(source,target)


def selected_files():
    files=[]
    for folder in ("WP1","WP2","WP3","WP4","WP5","WP6"):
        files.extend((ROOT/folder).glob("v8_*"))
    files.extend([ROOT/"WP1/量纲与假设表.csv",ROOT/"WP1/七点电导往返记录.csv",
                  ROOT/"WP2/预筛全部记录.csv",ROOT/"WP2/预筛可行集.json"])
    files.extend((ROOT/"src/effective_network_model").glob("v8_*.py"))
    files.extend([ROOT/"src/effective_network_model/v7_network_geometry.py",
                  ROOT/"src/data_ingest/v8_step1_units.py"])
    return sorted(set(files))


def main():
    if TARGET.resolve().parent != ROOT.resolve() or TARGET.exists():
        raise RuntimeError("formal archive target must be new and directly under project root")
    files=selected_files()
    figures=sorted((ROOT/"data/processed/figures").glob("v8_Fig*.svg"))+sorted(
        (ROOT/"data/processed/figures").glob("v8_Fig*.png"))
    figure_sources=[ROOT/p for p in (
        "WP1/量纲与假设表.csv","WP2/预筛全部记录.csv",
        "WP4/v8_全部逐点预测.csv","WP4/v8_模型比较.csv",
        "WP5/v8_固定十情景汇总.csv","WP5/v8_域尺寸汇总.csv")]
    checks=[ROOT/p for p in ("WP3/v8_技术检查.csv","WP2/预筛可行集.json",
                             "WP5/v8_主张判定.json")]
    required=[*files,*figures,*figure_sources,*checks,
              ROOT/"WP6/v8_actual_spec.json",ROOT/"WP6/v8_source_manifest.csv",
              ROOT/"WP6/v8_run_manifest.json",ROOT/"WP6/v8_导师决策一页纸.md"]
    missing=[str(p) for p in required if not p.is_file()]
    if missing or len(figures)!=12 or not TASKBOOK.is_dir():
        raise RuntimeError(f"incomplete archive inputs: missing={missing}, figures={len(figures)}")
    TARGET.mkdir()
    shutil.copytree(TASKBOOK,TARGET/"taskbook")
    for source in files:
        copy_file(source,TARGET/source.relative_to(ROOT))
    for source in figures:
        copy_file(source,TARGET/"figures"/source.name)
    for source in figure_sources:
        copy_file(source,TARGET/"figure_data"/source.name)
    for source in checks:
        copy_file(source,TARGET/"checks"/source.name)
    for name,source in (("actual_spec.json",ROOT/"WP6/v8_actual_spec.json"),
                        ("source_manifest.csv",ROOT/"WP6/v8_source_manifest.csv"),
                        ("00_导师先读.md",ROOT/"WP6/v8_导师决策一页纸.md")):
        copy_file(source,TARGET/name)
    manifest=json.loads((ROOT/"WP6/v8_run_manifest.json").read_text(encoding="utf-8"))
    manifest["states"]["formal_17_archive"]="CREATED_PENDING_REPLAY"
    (TARGET/"run_manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    (TARGET/"environment.txt").write_text(
        f"Python {sys.version.split()[0]}\nNumPy 2.5.4\nSciPy 1.18.1\nMatplotlib 3.11.2\n",
        encoding="utf-8")
    (TARGET/"README.md").write_text(
        "# v8.0.1 formal research archive\n\n"
        "Start with `00_导师先读.md`, then `WP6/v8_manuscript.md` and "
        "`WP6/v8_supporting_information.md`.\n\n"
        "Frozen source inputs are in `taskbook/`; generated evidence is in WP1-WP6. "
        "The scripts read `taskbook/` when running from this archive. "
        "With an existing compatible interpreter, replay `python src/data_ingest/v8_step1_units.py` "
        "and `python src/effective_network_model/v8_verify.py` from this directory. "
        "The source project used its local `.venv` and did not package an environment. "
        "The archive's root `run_manifest.json` starts as CREATED_PENDING_REPLAY.\n",
        encoding="utf-8")
    print(f"created {TARGET}; task files={len(files)}; figures={len(figures)}; representative replay pending")


if __name__=="__main__":
    main()
