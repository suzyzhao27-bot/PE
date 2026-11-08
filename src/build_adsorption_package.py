"""Assemble the next advisor-review folder without changing source records.

The builder is intentionally separate from the scientific calculation.  It can
be run only after the dry/wet simulation stage has been verified; until then it
is safe to exercise in a temporary directory through the unit tests.
"""

from __future__ import annotations

import csv
import importlib.util
import json
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PACKAGE_NAME = "32_拉伸结构与小分子吸附_导师侧计算审阅包_20261109_v1.0.1"
OUTPUT_ROOT = ROOT / "原始资料" / "赵珈绮"
DATA_ROOT = ROOT / "data" / "processed" / "surface_adsorption"
BASELINE_MODULE = ROOT / "src" / "mechanism_surface_adsorption_baseline.py"


def _resolve_work_root() -> Path:
    """Find the sibling directory that holds the supplied task materials.

    The working copy is commonly ``E:/PE/PE-main`` while a GitHub clone may
    be placed directly under ``E:/``.  Supporting both layouts keeps the
    package builder reproducible without copying source inputs into the repo.
    """
    candidates = (ROOT.parent, ROOT.parent / "PE")
    marker = "20261109_珈绮计算任务_拉伸结构与小分子吸附.md"
    for candidate in candidates:
        if (candidate / marker).is_file():
            return candidate
    return ROOT.parent


WORK_ROOT = _resolve_work_root()


def _load_baseline_module():
    spec = importlib.util.spec_from_file_location("surface_adsorption_baseline", BASELINE_MODULE)
    if spec is None or spec.loader is None:
        raise ImportError(BASELINE_MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def structure_rows() -> list[dict[str, object]]:
    return _load_baseline_module().structure_rows()


def build_baseline_rows() -> list[dict[str, object]]:
    values = _load_baseline_module().compute_baseline()
    return [
        {"comparison": "thin_over_thick", "ratio": values["thin_over_thick"], "delta_delta_g_kj_mol": values["thin_over_thick_delta_delta_g_kj_mol"]},
        {"comparison": "middle_over_thick", "ratio": values["middle_over_thick"], "delta_delta_g_kj_mol": values["middle_over_thick_delta_delta_g_kj_mol"]},
        {"comparison": "thin_over_middle", "ratio": values["thin_over_middle"], "delta_delta_g_kj_mol": values["thin_over_middle_delta_delta_g_kj_mol"]},
        {"comparison": "noncrystalline_bound", "ratio": values["noncrystalline_bound"], "delta_delta_g_kj_mol": ""},
    ]


def _copy(source: Path, destination: Path) -> None:
    if not source.is_file():
        raise FileNotFoundError(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _root_documents() -> dict[str, str]:
    return {
        "00_导师先读.txt": """32号导师先读：拉伸 PE 外表面小分子吸附

2026-11-09｜导师侧计算推进包。先看本页，再看 `01_科研推进与复算记录.md`、`02_计算边界与输入缺口.txt` 和 `figures/`。

本轮完成三组结构约束整理、复用已有纤维几何接口得到外层代理、随附九组荷电基线逐列复算、位点密度比与 -RT ln(R) 量级标尺。face_charge.py 的 9 行/20 个共同字段复现最大绝对差为 1.42×10^-14。

本轮没有生成分子吸附自由能。E:\\PE 中没有三组各三份原子表面、辛酸/辛酸根拓扑、PE 力场或 GROMACS/WHAM/MBAR，因此 W(z)、K_dry、R_dry、R_wet 均保持未产生。图2是预计算量级标尺，不是吸附曲线。

这是阶段性导师审阅材料，不把宏观结构表、ζ 条件反演或计划中的伞形采样写成辛酸吸附实测。
""",
        "01_科研推进与复算记录.md": """# 01 科研推进与复算记录

详见 `checks/20261109_表面吸附任务_基线复算与模拟边界.md` 和 `figure_data/20261109_复算记录与科研边界.md`。

当前已完成：三组结构约束、几何外层代理、荷电基线逐列复现、自由能量级换算、干态输入合同、湿态启动门、图1/图2/图3预计算图。

当前停止：E:\\PE 和 GitHub 公共仓库均没有原子表面、辛酸/PE 力场、GROMACS/WHAM/MBAR；因此不生成 W(z)、K_dry、R_dry 或 R_wet。
""",
        "02_计算边界与输入缺口.txt": """02 计算边界与输入缺口

必须补齐后才启动干态 MD：
1. 182/357/604 三组各 3 份独立原子表面坐标和构模种子；
2. PE 力场、辛酸拓扑/电荷、晶体取向、周期盒和松弛协议；
3. 真实外表面面积算法、第一吸附层边界、伞形窗口与时间块设置；
4. 可运行的 GROMACS 与 WHAM/MBAR 工具链。

缺失项全部保持 unknown。ζ、接触角、AFM 粗糙度和表观表面能不能替代这些输入。
""",
        "00_导师先读.md": """# 32号导师先读：拉伸 PE 外表面小分子吸附

2026-11-09｜导师侧计算审阅包。先看本页，再看 `checks/20261109_复算记录与科研边界.md`、`figures/` 和 `figure_data/`。

本轮已经完成三组结构约束整理、复用已有纤维几何接口得到外层代理、随附九组荷电基线逐列复算、位点密度比与 `-RT ln(R)` 量级标尺。`face_charge.py` 的 9 行/20 个共同字段复现最大绝对差为 `1.42×10^-14`。

本轮没有生成分子吸附自由能。E:\\PE 中没有三组各三份原子表面、辛酸/辛酸根拓扑、PE 力场或 GROMACS/WHAM/MBAR，因此 `W(z)`、`K_dry`、`R_dry`、`R_wet` 均保持未产生。图2是预计算量级标尺，不是吸附曲线。

包内把实验/条件反演、任务书约束、模型计划和 unknown 分开存放。`1.84` 是 ζ 条件反演的有效位点密度比，`1.29` 是固定两相吸附性质的简单混合上限，均不能代替辛酸吸附能。取得原子结构和工具链后，按 `checks/20261109_干态吸附计算输入合同.md` 先跑干态；只有干态方向成立才启动湿态 HA/A^-。
""",
        "README.md": """# 32号导师侧计算审阅包导读

阅读顺序：

1. `00_导师先读.md`；
2. `checks/` 中的复算记录、干态输入合同和湿态启动门；
3. `figures/图1_结构约束与模拟边界.*`、`figures/图2_量级标尺与交付状态.*`、`figures/图3_几何外层代理.*`；
4. `figure_data/` 的轻量 CSV/JSON/复算输出；
5. `source_inputs/` 的任务书、荷电结论、九组 CSV 和原脚本；
6. `src/` 的复算与打包源码。

当前包是导师侧科研推进记录。它不把宏观结构表、几何外层代理、ζ 条件反演或计划中的伞形采样写成辛酸吸附实测。所有 `unknown` 字段保留，缺少原子结构/力场/MD 工具时不生成 `W(z)` 或 `R_dry`。
""",
    }


def build_package(destination_root: Path | None = None) -> Path:
    destination_root = OUTPUT_ROOT if destination_root is None else Path(destination_root)
    package = destination_root / PACKAGE_NAME
    if package.exists():
        raise FileExistsError(package)
    package.mkdir(parents=True)

    source_inputs = {
        WORK_ROOT / "20261109_珈绮计算任务_拉伸结构与小分子吸附.md": "source_inputs/20261109_珈绮计算任务_拉伸结构与小分子吸附.md",
        WORK_ROOT / "20261109_赵珈绮晨会提纲_阻盐机制与必要实验.md": "source_inputs/20261109_赵珈绮晨会提纲_阻盐机制与必要实验.md",
        WORK_ROOT / "20261109_膜外表面吸附荷电_第二步核算_v1.0.0/01_膜外表面吸附荷电_结论与依据.md": "source_inputs/01_膜外表面吸附荷电_结论与依据.md",
        WORK_ROOT / "20261109_膜外表面吸附荷电_第二步核算_v1.0.0/02_九组计算与实验.csv": "source_inputs/02_九组计算与实验.csv",
        WORK_ROOT / "20261109_膜外表面吸附荷电_第二步核算_v1.0.0/face_charge.py": "source_inputs/face_charge.py",
    }
    for source, rel in source_inputs.items():
        _copy(source, package / rel)

    if not DATA_ROOT.is_dir():
        raise FileNotFoundError(DATA_ROOT)
    for source in sorted(DATA_ROOT.iterdir()):
        if source.is_file() and source.name != "manifest.csv":
            folder = "figures" if source.suffix.lower() in {".png", ".svg", ".pdf"} else "figure_data"
            _copy(source, package / folder / source.name)

    geometry_input = ROOT / "WP3" / "v9_几何层表.csv"
    if geometry_input.is_file():
        _copy(geometry_input, package / "figure_data" / "20261109_几何层表_输入.csv")

    check_files = {
        ROOT / "阶段总结/20261109_表面吸附任务_基线复算与模拟边界.md": "checks/20261109_表面吸附任务_基线复算与模拟边界.md",
        ROOT / "data/processed/surface_adsorption/20261109_干态吸附计算输入合同.md": "checks/20261109_干态吸附计算输入合同.md",
        ROOT / "data/processed/surface_adsorption/20261109_湿态启动门与停止规则.md": "checks/20261109_湿态启动门与停止规则.md",
        ROOT / ".planning/2026-11-09-表面吸附任务/progress.md": "checks/20261109_进度检查点.md",
        ROOT / ".planning/2026-11-09-表面吸附任务/findings.md": "checks/20261109_发现与决策.md",
        ROOT / "项目进展计划与任务看板.md": "checks/项目进展计划与任务看板.md",
    }
    for source, rel in check_files.items():
        _copy(source, package / rel)

    source_files = [
        ROOT / "src/mechanism_surface_adsorption_baseline.py",
        ROOT / "src/surface_geometry_proxy.py",
        ROOT / "src/test_surface_adsorption_baseline.py",
        ROOT / "src/test_surface_geometry_proxy.py",
        ROOT / "src/build_adsorption_package.py",
        ROOT / "src/test_build_adsorption_package.py",
    ]
    for source in source_files:
        _copy(source, package / "src" / source.name)

    for name, content in _root_documents().items():
        _write(package / name, content)
    scope = {
        "package_number": 32,
        "package_name": PACKAGE_NAME,
        "date": "2026-11-09",
        "role": "ADVISOR_SIDE_COMPUTATIONAL_PROGRESS",
        "dry_md_status": "BLOCKED_UNTIL_ATOMIC_SURFACES_FORCE_FIELD_AND_MD_TOOLCHAIN",
        "wet_md_status": "NOT_STARTED_BY_DESIGN",
        "face_charge_reproduction": "PASS",
        "geometry_proxy_status": "COMPUTED_STRUCTURAL_SCENARIO_ONLY",
        "geometry_input_packaged": True,
        "source_inputs_unchanged": True,
        "experimental_claims_added": False,
        "synthetic_adsorption_curves_added": False,
    }
    _write(package / "scope.json", json.dumps(scope, ensure_ascii=False, indent=2) + "\n")

    def collect_manifest_rows() -> list[dict[str, str]]:
        rows = []
        for path in sorted(package.rglob("*")):
            if path.is_file() and path.name != "manifest.csv":
                rel = path.relative_to(package).as_posix()
                role = rel.split("/", 1)[0]
                rows.append({"relative_path": rel, "role": role, "status": "COPIED_OR_GENERATED"})
        return rows

    manifest_rows = collect_manifest_rows()
    _write(
        package / "03_文件清单.txt",
        "03 文件清单\n\n"
        + "\n".join([row["relative_path"] for row in manifest_rows] + ["03_文件清单.txt"])
        + "\n",
    )
    manifest_rows = collect_manifest_rows()
    with (package / "manifest.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["relative_path", "role", "status"])
        writer.writeheader()
        writer.writerows(manifest_rows)
    return package


if __name__ == "__main__":
    built = build_package()
    print(f"package={built}")
