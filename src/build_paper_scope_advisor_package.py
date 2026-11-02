"""Assemble a compact local advisor review ZIP without changing source records."""

from __future__ import annotations

import csv
import re
from collections import Counter
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/salt_transport"
PAPER = ROOT / "paper"
SUMMARY = ROOT / "阶段总结"
TRACE = PAPER / "review-traces"
ORIGINAL = ROOT / "原始资料/赵珈绮"
PACKAGE_NAME = "25_论文贡献裁剪_导师侧计算审阅包_20261031_v1.0.0"
OUT = PAPER / f"{PACKAGE_NAME}.zip"
DRAFT = PAPER / "20261031_Stage2_B1_论文计算部分审阅稿.md"


def main() -> None:
    if OUT.exists():
        raise FileExistsError(f"Review package already exists: {OUT}")
    files: dict[Path, str] = {}
    dests: set[str] = set()

    def add(path: Path, dest: str) -> None:
        path = path.resolve()
        if not path.is_file() or dest in dests or path in files:
            raise ValueError(f"Missing or duplicate package member: {path} -> {dest}")
        files[path] = dest
        dests.add(dest)

    add(PAPER / "20261031_25号论文裁剪_导师先读.md", "00_导师先读.md")
    add(ROOT / "requirements.txt", "requirements.txt")
    for path in DATA.glob("25_论文裁剪*.csv"):
        add(path, f"figure_data/{path.name}")
    add(DATA / "23_成对效应膜阻口径勘误_20261030.csv",
        "figure_data/23_成对效应膜阻口径勘误_20261030.csv")
    for path in DATA.glob("25_论文裁剪*.png"):
        add(path, f"figures/{path.name}")
    for path in (ROOT / "src").glob("mechanism_paper_*.py"):
        add(path, f"src/{path.name}")
    for name in ("mechanism_paired_effects_entry_correction.py",
                 "mechanism_neutral_baseline.py", "render_paper_scope_figures.py",
                 "build_paper_scope_advisor_package.py"):
        add(ROOT / "src" / name, f"src/{name}")

    add(SUMMARY / "20261030_23号成对效应膜阻口径勘误.md",
        "WP1/20261030_23号成对效应膜阻口径勘误.md")
    for path in SUMMARY.glob("20261031_论文*.md"):
        name = path.name
        if "论文A_D图源预览" in name or "论文裁剪计算范围收口" in name or "论文贡献裁剪包接续核对" in name:
            folder = "checks"
        elif "论文A" in name or "论文裁剪A" in name:
            folder = "WP2"
        elif "论文裁剪B" in name:
            folder = "WP3"
        elif "论文D" in name or "论文裁剪D" in name:
            folder = "WP4"
        else:
            raise ValueError(f"Unmapped stage summary: {name}")
        add(path, f"{folder}/{name}")
    for path in TRACE.glob("20261031_*.md"):
        add(path, f"checks/{path.name}")
    add(TRACE / "20261030_23号膜阻口径勘误审查.md",
        "checks/20261030_23号膜阻口径勘误审查.md")

    inputs = {
        "20_interval_transport.csv": ORIGINAL / "25_数据约束机制研究_学生任务书_v1.0.1/inputs/20/interval_transport.csv",
        "25_原矩阵论文用途子集.json": ORIGINAL / "25_论文贡献审计与执行裁剪_v1.0.0/原矩阵论文用途子集.json",
        "25_论文贡献可行性核查.json": ORIGINAL / "25_论文贡献审计与执行裁剪_v1.0.0/论文贡献可行性核查.json",
        "22_paired_predictions.json": ORIGINAL / "22_机制预测与实验验证方案_v1.0.0/matrices/paired_predictions.json",
        "22_paired_effects_20261028.csv": DATA / "22_paired_effects_20261028.csv",
        "22_pb_main_20261028.csv": ORIGINAL / "25_数据约束机制研究_学生任务书_v1.0.1/inputs/23/22_pb_main_20261028.csv",
        "22_pb_target_recheck_20261028.csv": ORIGINAL / "25_数据约束机制研究_学生任务书_v1.0.1/inputs/23/22_pb_target_recheck_20261028.csv",
    }
    for name, path in inputs.items():
        add(path, f"figure_data/source_inputs/{name}")

    # Keep the manuscript's relative links valid within the compact ZIP layout.
    manuscript = DRAFT.read_text(encoding="utf-8")
    link_pattern = re.compile(r"(\]\()([^)]+)(\))")

    def relink(match: re.Match[str]) -> str:
        target = match.group(2)
        if "://" in target or target.startswith("#"):
            return match.group(0)
        source = (DRAFT.parent / target).resolve()
        if source not in files:
            raise ValueError(f"Manuscript link omitted from package: {target}")
        return f"]({files[source]})"

    manuscript = link_pattern.sub(relink, manuscript)
    linked = re.findall(r"\]\(([^)]+)\)", manuscript)
    if any(link not in dests for link in linked):
        raise ValueError("Rewritten manuscript link missing")
    counts = Counter(dest.split("/", 1)[0] for dest in dests)
    review_note = (
        "# 包内核对记录\n\n"
        "2026-10-31。此包仅为导师侧计算审阅，学生本人复算和WP5/6未完成。\n\n"
        f"源文件副本 {len(files)} 份；另含导师先读以外的包内导读、论文稿和本核对记录。"
        f"论文稿 {len(linked)} 个相对链接均指向包内文件。\n\n"
        "目录计数（源文件）：" + "、".join(f"{key} {value}" for key, value in sorted(counts.items())) + "。\n\n"
        "已定向核对：23号勘误54行；A五目标80行；B诊断12行；D主算子浓度48行、响应24行；"
        "D原定±5%速率比3行、浓度24行、响应12行。两张PNG仅为审阅预览。"
        "原T12十项未决、旧装置元数据、学生本人验收和实验状态均未改。"
        "无新PB、时序、LP或实验；未发送、上传或同步OneDrive。\n"
    )
    readme = (
        "# 导师审阅包导读\n\n"
        "先读`00_导师先读.md`，再看`01_论文计算部分审阅稿.md`和`figures/`。"
        "`figure_data/`存逐行结果，`figure_data/source_inputs/`存最小冻结输入副本；"
        "`WP1`为23号系数勘误，`WP2`为A能力界，`WP3`为B速率诊断，`WP4`为D中性探针；"
        "`checks/`存独立审查和完成边界，`src/`存对应源码。\n\n"
        "源输入原件位于项目`原始资料/赵珈绮/25_数据约束机制研究_学生任务书_v1.0.1`、"
        "`25_论文贡献审计与执行裁剪_v1.0.0`及22号原方案。"
        "包内源码保留项目相对路径，独立运行须置回原项目结构；本包用于导师阅读，"
        "不冒充脱离原资料即可完全复现的学生成果。\n"
    )

    # Validate current decisive table identities before creating an archive.
    expected = {
        "23_成对效应膜阻口径勘误_20261030.csv": 54,
        "25_论文裁剪A五目标界_20261031.csv": 80,
        "25_论文裁剪B速率范围诊断_20261031.csv": 12,
        "25_论文裁剪D尿素浓度条件范围_20261031.csv": 48,
        "25_论文裁剪D尿素处理响应范围_20261031.csv": 24,
        "25_论文裁剪D迁移系数敏感性_速率比_20261031.csv": 3,
        "25_论文裁剪D迁移系数敏感性_尿素浓度_20261031.csv": 24,
        "25_论文裁剪D迁移系数敏感性_处理响应_20261031.csv": 12,
    }
    for name, count in expected.items():
        with (DATA / name).open(encoding="utf-8-sig", newline="") as stream:
            if sum(1 for _ in csv.DictReader(stream)) != count:
                raise ValueError(f"Unexpected row count: {name}")
    with ZipFile(OUT, "x", ZIP_DEFLATED, compresslevel=6) as archive:
        for source, dest in sorted(files.items(), key=lambda item: item[1]):
            archive.write(source, f"{PACKAGE_NAME}/{dest}")
        for dest, content in (
            ("README.md", readme),
            ("01_论文计算部分审阅稿.md", manuscript),
            ("checks/包内核对记录.md", review_note),
        ):
            archive.writestr(f"{PACKAGE_NAME}/{dest}", content)
    with ZipFile(OUT) as archive:
        names = [name.removeprefix(f"{PACKAGE_NAME}/") for name in archive.namelist()]
        if len(names) != len(files) + 3 or len(names) != len(set(names)):
            raise ValueError("Archive member count or names mismatch")
        if any(link not in names for link in linked):
            raise ValueError("Archived manuscript link missing")
        if not archive.read(f"{PACKAGE_NAME}/01_论文计算部分审阅稿.md"):
            raise ValueError("Archived manuscript unreadable")
    print(f"package={OUT.name} members={len(files)+3} linked={len(linked)} bytes={OUT.stat().st_size}")


if __name__ == "__main__":
    main()
