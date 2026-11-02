"""Build a standalone advisor reading ZIP for the 27-review follow-up."""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/salt_transport"
PAPER = ROOT / "paper"
ORIGINAL = ROOT / "原始资料/赵珈绮"
REVIEW = ORIGINAL / "27_26号成果审阅_20261031_v1.0.0"
NAME = "28_27号审阅落实_导师侧补充审阅包_v1.0.1"
OUT = ORIGINAL / f"{NAME}.zip"
DRAFT = PAPER / "20261031_Stage2_B1_论文计算部分审阅稿.md"


def main() -> None:
    if OUT.exists():
        raise FileExistsError(OUT)
    files: dict[Path, str] = {}
    destinations: set[str] = set()

    def add(path: Path, dest: str) -> str:
        path = path.resolve()
        if path in files:
            if files[path] != dest:
                raise ValueError(f"Conflicting destinations: {path}")
            return dest
        if not path.is_file() or dest in destinations:
            raise ValueError(f"Missing or duplicate member: {path} -> {dest}")
        files[path] = dest
        destinations.add(dest)
        return dest

    def linked_dest(path: Path) -> str:
        path = path.resolve()
        rel = path.relative_to(ROOT)
        if rel.parts[:3] == ("data", "processed", "salt_transport"):
            folder = "figures" if path.suffix.lower() == ".png" else "figure_data"
            return add(path, f"{folder}/{path.name}")
        if rel.parts[0] == "阶段总结":
            return add(path, f"checks/{path.name}")
        if rel.parts[:2] == ("paper", "review-traces"):
            return add(path, f"checks/{path.name}")
        if path.name == "22_考试后实验与判别.md":
            return add(path, f"protocol/{path.name}")
        raise ValueError(f"Unmapped manuscript link: {path}")

    manuscript = DRAFT.read_text(encoding="utf-8")
    links: list[str] = []

    def relink(match: re.Match[str]) -> str:
        target = match.group(1)
        if "://" in target or target.startswith("#"):
            return match.group(0)
        dest = linked_dest(DRAFT.parent / target)
        links.append(dest)
        return f"]({dest})"

    manuscript = re.sub(r"\]\(([^)]+)\)", relink, manuscript)
    for name in (
        "25_论文裁剪A主图预览_20261031.png",
        "25_论文裁剪D主图预览_20261031.png",
    ):
        add(DATA / name, f"figures/{name}")
    for path in DATA.glob("25_论文裁剪*.csv"):
        add(path, f"figure_data/{path.name}")
    add(DATA / "23_成对效应膜阻口径勘误_20261030.csv",
        "figure_data/23_成对效应膜阻口径勘误_20261030.csv")
    for path in (ROOT / "src").glob("mechanism_paper_*.py"):
        add(path, f"src/{path.name}")
    for name in (
        "mechanism_paired_effects_entry_correction.py",
        "mechanism_neutral_baseline.py",
        "render_paper_scope_figures.py",
        "build_27_review_followup_package.py",
    ):
        add(ROOT / "src" / name, f"src/{name}")
    add(ROOT / "requirements.txt", "requirements.txt")
    for path in (ROOT / "阶段总结").glob("20261101_*.md"):
        add(path, f"checks/{path.name}")
    for name in (
        "00_26号成果审阅意见.md",
        "独立复核结果.json",
        "独立逐项数值核对.json",
        "核查源码.py.txt",
    ):
        add(REVIEW / name, f"review_source/{name}")
    source_inputs = {
        "20_interval_transport.csv": ORIGINAL / "25_数据约束机制研究_学生任务书_v1.0.1/inputs/20/interval_transport.csv",
        "25_原矩阵论文用途子集.json": ORIGINAL / "25_论文贡献审计与执行裁剪_v1.0.0/原矩阵论文用途子集.json",
        "25_论文贡献可行性核查.json": ORIGINAL / "25_论文贡献审计与执行裁剪_v1.0.0/论文贡献可行性核查.json",
        "22_paired_predictions.json": ORIGINAL / "22_机制预测与实验验证方案_v1.0.0/matrices/paired_predictions.json",
        "22_paired_effects_20261028.csv": DATA / "22_paired_effects_20261028.csv",
        "22_pb_main_20261028.csv": ORIGINAL / "25_数据约束机制研究_学生任务书_v1.0.1/inputs/23/22_pb_main_20261028.csv",
        "22_pb_target_recheck_20261028.csv": ORIGINAL / "25_数据约束机制研究_学生任务书_v1.0.1/inputs/23/22_pb_target_recheck_20261028.csv",
    }
    for name, path in source_inputs.items():
        add(path, f"figure_data/source_inputs/{name}")

    expected_rows = {
        "23_成对效应膜阻口径勘误_20261030.csv": 54,
        "25_论文裁剪A五目标界_20261031.csv": 80,
        "25_论文裁剪B速率范围诊断_20261031.csv": 12,
        "25_论文裁剪D尿素浓度条件范围_20261031.csv": 48,
        "25_论文裁剪D尿素处理响应范围_20261031.csv": 24,
        "25_论文裁剪D迁移系数敏感性_速率比_20261031.csv": 3,
        "25_论文裁剪D迁移系数敏感性_尿素浓度_20261031.csv": 24,
        "25_论文裁剪D迁移系数敏感性_处理响应_20261031.csv": 12,
    }
    for name, count in expected_rows.items():
        with (DATA / name).open(encoding="utf-8-sig", newline="") as stream:
            if sum(1 for _ in csv.DictReader(stream)) != count:
                raise ValueError(f"Unexpected row count: {name}")
    if any(not (DATA / name).is_file() for name in (
        "25_论文裁剪A主图预览_20261031.png",
        "25_论文裁剪D主图预览_20261031.png",
    )):
        raise ValueError("Missing review figure")

    cover = """# 28号导师先读：27号审阅后的补充材料

2026-11-01｜Stage 2 · B1｜导师侧补充审阅。先读本页，再看`01_论文计算部分审阅稿.md`及`figures/`。本包可独立阅读修订稿、图、关键逐行表和审阅依据；26号原审阅包与27号原件保持原状。

27号审阅后，A图从既有80行加入G0/E1/δ=.20参照，与G1同δ展示；两者几何及允许荷电范围仍不同，不能把差额单归电荷。D图保留厚膜±5%区间覆盖1，标明外阻取消条件。论文稿补图注、原WP6判据及旧R/新E身份；D浓度的等体积、无取样、无净流现明确为模型假设，旧逐运行记录未认证这些条件。

审阅者报告550项独立数值比较通过，最大绝对差约1.60×10⁻¹²；这是其表格与代数核对，不是本包重新运行PB、学生代码或实验。A的1358.48/1369.88及Γ=10下135.85/136.99为既定条件的必要能力下界，不是已测孔数。D的±5%为敏感性情景，不是已证实误差界。

**本次交付范围：** 研究者明确不要求学生运行记录和本人解释。这两项不作为28号导师侧交付门槛；本包也不声称学生已执行。A/D图已成修订审阅版并附图注草案，论文计算部分已有导师侧审阅稿。代表性本机复现入口见`03_代表性复现入口.md`。详见`02_27号意见逐项回应.md`。

T12十项目标保持数值未决，不新增PB。WP5/6仍未开展；不把旧60/100 min预测当新实验靶值。本包不是学生正式成果、机制已验证的论文或实验报告；未代替原26号包的已审阅版本。
"""
    response = """# 27号审阅意见逐项回应

| 27号指示 | 本次处理 | 交付边界 |
|---|---|---|
| 学生运行记录与解释 | 研究者决定28号导师侧交付不要求这两项；不收取、不代签 | 不声称学生已执行或解释 |
| A/D两图 | `figures/`为现有逐行表重绘的修订审阅图；A补同δ参照，D保留比值1和取消条件；正文有图注草案 | 最终版式仍供导师审阅 |
| 论文与复现 | `01_论文计算部分审阅稿.md`补理论、结果、竞争解释和原WP6接口；`03_代表性复现入口.md`列本机命令与输入输出位置 | 是导师侧计算稿；脚本仍需原项目目录结构 |

另补正旧Book2逐运行体积、取样与净流身份：D浓度公式只在等体积、无取样、无净流模型条件下成立。T12未决；按原25号裁剪不扩算，按原22号方案等待WP5/6实验。
"""
    readme = """# 包内导读

这是27号审阅后、可单独发导师阅读的28号补充ZIP。`01_论文计算部分审阅稿.md`的本地链接已改写为包内路径。`figures/`为两张修订预览；`figure_data/`为轻量结果与关键源输入；`src/`为对应脚本；`checks/`为推导、证据边界和本轮修订；`protocol/`为原WP5/6规则；`review_source/`为27号审阅原件副本。较早的审查记录按当时任务要求保留；本次交付门槛以`scope_override.json`和`02_27号意见逐项回应.md`为准。

脚本保持原项目路径结构；从此ZIP单独解压不能直接执行全部代码。代表性本机复现命令、数据根和输出位置见`03_代表性复现入口.md`。研究者决定本次导师侧交付不要求学生记录或解释；本包不声称学生已执行。原26号包未被覆盖或重打包。
"""
    reproduce = """# 代表性复现入口

本机项目根：`D:\\超薄 PE 纤维膜的界面荷电、离子选择性与盐输运`。在完整项目结构、已安装`requirements.txt`且存在项目`.venv`时，运行：

```powershell
Set-Location -LiteralPath 'D:\\超薄 PE 纤维膜的界面荷电、离子选择性与盐输运'
& '.venv\\Scripts\\python.exe' 'src\\mechanism_paper_d_mobility_sensitivity.py'
```

脚本读取项目`原始资料/赵珈绮/25_论文贡献审计与执行裁剪_v1.0.0/`的两个JSON，及`data/processed/salt_transport/`中的D浓度、响应基线CSV；输出同一数据目录的三张`25_论文裁剪D迁移系数敏感性_*.csv`。已保存的导师侧运行结果行数为速率比3、浓度24、响应12，其中8行在±5%后仍排除比值1。本命令是项目本机示例，ZIP内源码保留原路径，单独解压不能直接运行。无新PB、LP或实验。
"""
    scope = {
        "package_number": 28,
        "package_name": NAME,
        "role": "ADVISOR_SIDE_27_REVIEW_FOLLOWUP",
        "base_package": "26_数据约束机制研究_导师侧计算审阅包_v1.0.0",
        "review_package": REVIEW.name,
        "original_task": "25_数据约束机制研究_学生任务书_v1.0.1",
        "scope_override": "25_论文贡献审计与执行裁剪_v1.0.0",
        "advisor_side_figures_and_manuscript": "REVISED_FOR_REVIEW",
        "student_independent_implementation": "NOT_RUN",
        "student_signoff": "NOT_RUN",
        "student_records_and_explanation_required_for_this_delivery": False,
        "student_requirement_decision": "WAIVED_BY_RESEARCHER_FOR_28_ADVISOR_DELIVERY",
        "T12_integral_targets": "10_NUMERIC_UNRESOLVED",
        "WP5_experiment": "NOT_RUN",
        "WP6_inference": "NOT_RUN",
        "source_packages_unchanged": True,
        "sent_or_uploaded": False,
    }
    note = (
        "# 包内核对记录\n\n"
        f"源文件副本{len(files)}份，论文稿{len(links)}个相对链接均指包内文件。"
        "23号勘误54行，A五目标80行，B诊断12行，D主浓度48行、响应24行，"
        "D原定±5%速率比3行、浓度24行、响应12行。"
        "两张PNG为导师侧修订审阅图。未运行新PB、LP、时序或实验；本次交付不要求学生记录与解释，亦不声称其已执行。\n"
    )
    generated = {
        "00_导师先读.md": cover,
        "01_论文计算部分审阅稿.md": manuscript,
        "02_27号意见逐项回应.md": response,
        "03_代表性复现入口.md": reproduce,
        "README.md": readme,
        "scope_override.json": json.dumps(scope, ensure_ascii=False, indent=2) + "\n",
        "checks/包内核对记录.md": note,
    }
    if any(dest in destinations for dest in generated):
        raise ValueError("Generated member conflicts with source")
    with ZipFile(OUT, "x", ZIP_DEFLATED, compresslevel=6) as archive:
        for path, dest in sorted(files.items(), key=lambda item: item[1]):
            archive.write(path, f"{NAME}/{dest}")
        for dest, content in generated.items():
            archive.writestr(f"{NAME}/{dest}", content)
    with ZipFile(OUT) as archive:
        names = [p.removeprefix(f"{NAME}/") for p in archive.namelist()]
        if len(names) != len(files) + len(generated) or len(names) != len(set(names)):
            raise ValueError("Archive members missing or duplicated")
        if any(link not in names for link in links):
            raise ValueError("Manuscript link missing from ZIP")
        for dest in generated:
            if not archive.read(f"{NAME}/{dest}"):
                raise ValueError(f"Unreadable archive member: {dest}")
    print(f"package={OUT.name} members={len(files)+len(generated)} links={len(links)} bytes={OUT.stat().st_size}")


if __name__ == "__main__":
    main()
