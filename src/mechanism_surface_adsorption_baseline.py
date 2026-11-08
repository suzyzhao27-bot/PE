"""Reproduce the supplied face-charge reference and freeze adsorption-task inputs.

This module deliberately stops before molecular dynamics.  It computes only
quantities that are supported by the task brief and the supplied nine-row CSV:
structure constraints, charge-density ratios, the RT ln(R) scale, the simple
non-crystalline mixture bound, and a row-by-row reproduction check for the
existing ``face_charge.py`` calculation.
"""

from __future__ import annotations

import argparse
import csv
import importlib.metadata
import math
import shutil
import subprocess
import sys
from pathlib import Path


R = 8.31446261815324
T = 298.15

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _resolve_work_root() -> Path:
    """Locate supplied task materials in the common local layouts."""
    candidates = (PROJECT_ROOT.parent, PROJECT_ROOT.parent / "PE")
    marker = "20261109_膜外表面吸附荷电_第二步核算_v1.0.0"
    for candidate in candidates:
        if (candidate / marker).is_dir():
            return candidate
    return PROJECT_ROOT.parent


WORK_ROOT = _resolve_work_root()
SUPPLIED_CSV = (
    WORK_ROOT
    / "20261109_膜外表面吸附荷电_第二步核算_v1.0.0"
    / "02_九组计算与实验.csv"
)
FACE_CHARGE = (
    WORK_ROOT
    / "20261109_膜外表面吸附荷电_第二步核算_v1.0.0"
    / "face_charge.py"
)
DEFAULT_OUTPUT = PROJECT_ROOT / "data" / "processed" / "surface_adsorption"


def structure_rows() -> list[dict[str, object]]:
    """Return the three-group structural contract from the task brief.

    Values are constraints for constructing representative surfaces.  They are
    not atomic coordinates and are intentionally labelled as task-brief inputs.
    """

    return [
        {
            "draw": "182",
            "second_draw_ratio": 1.5,
            "thickness_nm": 770.0,
            "dsc_crystallinity_mass_frac": 0.5399,
            "xrd_crystallinity_frac": 0.549,
            "coherent_110_nm": 16.20,
            "coherent_200_nm": 14.79,
            "fiber_diameter_nm": 24.54,
            "afm_roughness_nm": 18.1,
            "surface_energy_mJ_m2": 42.3,
            "zeta_mV": -52.9,
            "site_density_per_nm2": 0.04318916524277331,
            "evidence_role": "task_brief_constraint",
        },
        {
            "draw": "357",
            "second_draw_ratio": 2.1,
            "thickness_nm": 370.0,
            "dsc_crystallinity_mass_frac": 0.5108,
            "xrd_crystallinity_frac": 0.511,
            "coherent_110_nm": 15.07,
            "coherent_200_nm": 13.68,
            "fiber_diameter_nm": 22.56,
            "afm_roughness_nm": 16.0,
            "surface_energy_mJ_m2": 43.3,
            "zeta_mV": -59.1,
            "site_density_per_nm2": 0.0509063375654427,
            "evidence_role": "task_brief_constraint",
        },
        {
            "draw": "604",
            "second_draw_ratio": 2.73,
            "thickness_nm": 230.0,
            "dsc_crystallinity_mass_frac": 0.4083,
            "xrd_crystallinity_frac": 0.406,
            "coherent_110_nm": 14.48,
            "coherent_200_nm": 12.85,
            "fiber_diameter_nm": 26.58,
            "afm_roughness_nm": 23.2,
            "surface_energy_mJ_m2": 43.4,
            "zeta_mV": -70.5,
            "site_density_per_nm2": 0.07935339114996252,
            "evidence_role": "task_brief_constraint",
        },
    ]


def compute_baseline() -> dict[str, float]:
    rows = structure_rows()
    by_draw = {row["draw"]: row for row in rows}
    thick = float(by_draw["182"]["site_density_per_nm2"])
    middle = float(by_draw["357"]["site_density_per_nm2"])
    thin = float(by_draw["604"]["site_density_per_nm2"])
    ratios = {
        "thin_over_thick": thin / thick,
        "middle_over_thick": middle / thick,
        "thin_over_middle": thin / middle,
        "noncrystalline_bound": (1.0 - 0.4083) / (1.0 - 0.5399),
    }
    ratios["thin_over_thick_delta_delta_g_kj_mol"] = -R * T * math.log(ratios["thin_over_thick"]) / 1000.0
    ratios["middle_over_thick_delta_delta_g_kj_mol"] = -R * T * math.log(ratios["middle_over_thick"]) / 1000.0
    ratios["thin_over_middle_delta_delta_g_kj_mol"] = -R * T * math.log(ratios["thin_over_middle"]) / 1000.0
    return ratios


def build_audit_rows() -> list[dict[str, str]]:
    return [
        {
            "item": "zeta_site_density_reference",
            "status": "REPRODUCED_FROM_SUPPLIED_CSV",
            "evidence": "pH~7 rows in 02_九组计算与实验.csv",
            "claim_limit": "conditional reference only; not measured adsorbate identity",
        },
        {
            "item": "dry_neutral_octanoic_acid_w_of_z",
            "status": "BLOCKED_MISSING_MD_INPUTS",
            "evidence": "no atomic surface coordinates, octanoic-acid topology, or GROMACS executable in E:\\PE",
            "claim_limit": "do not report K_dry, R_dry, or W(z)",
        },
        {
            "item": "wet_HA_and_A_minus",
            "status": "NOT_STARTED_BY_DESIGN",
            "evidence": "task brief gates wet calculation on a dry-state direction result",
            "claim_limit": "do not report K_wet or R_wet",
        },
        {
            "item": "structure_surface_area_and_ordered_fraction",
            "status": "BLOCKED_MISSING_ATOMIC_STRUCTURES",
            "evidence": "only macro constraints are supplied; no three independent surface realizations",
            "claim_limit": "no group-specific surface morphology claim",
        },
        {
            "item": "legacy_square_fiber_geometry_proxy",
            "status": "REPRODUCED_GEOMETRY_PROXY",
            "evidence": "WP3/v9_几何层表.csv; layer counts 11/6/4 for 182/357/604",
            "claim_limit": "structural scenario only; not atomic accessible area or adsorption capacity",
        },
        {
            "item": "simple_noncrystalline_mixture_bound",
            "status": "REPRODUCED_FROM_TASK_BRIEF",
            "evidence": "(1-0.4083)/(1-0.5399)",
            "claim_limit": "upper bound under fixed phase-specific adsorption; not a surface measurement",
        },
    ]


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"No rows for {path}")
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def _write_structure_table(output_dir: Path) -> Path:
    rows = []
    for row in structure_rows():
        rows.append({str(key): value for key, value in row.items()})
    path = output_dir / "20261109_结构约束表.csv"
    _write_csv(path, rows)
    return path


def _write_baseline_table(output_dir: Path) -> Path:
    values = compute_baseline()
    rows = [
        {
            "comparison": "thin_over_thick",
            "numerator": "604",
            "denominator": "182",
            "ratio": values["thin_over_thick"],
            "delta_delta_g_kj_mol": values["thin_over_thick_delta_delta_g_kj_mol"],
            "reference_or_role": "zeta-inverted effective site-density ratio",
            "formula": "Gamma_604/Gamma_182; -RT ln(R)",
        },
        {
            "comparison": "middle_over_thick",
            "numerator": "357",
            "denominator": "182",
            "ratio": values["middle_over_thick"],
            "delta_delta_g_kj_mol": values["middle_over_thick_delta_delta_g_kj_mol"],
            "reference_or_role": "zeta-inverted effective site-density ratio",
            "formula": "Gamma_357/Gamma_182; -RT ln(R)",
        },
        {
            "comparison": "thin_over_middle",
            "numerator": "604",
            "denominator": "357",
            "ratio": values["thin_over_middle"],
            "delta_delta_g_kj_mol": values["thin_over_middle_delta_delta_g_kj_mol"],
            "reference_or_role": "zeta-inverted effective site-density ratio",
            "formula": "Gamma_604/Gamma_357; -RT ln(R)",
        },
        {
            "comparison": "noncrystalline_bound",
            "numerator": "604 noncrystalline fraction",
            "denominator": "182 noncrystalline fraction",
            "ratio": values["noncrystalline_bound"],
            "delta_delta_g_kj_mol": "",
            "reference_or_role": "fixed phase-specific adsorption upper bound",
            "formula": "(1-Xc_604)/(1-Xc_182)",
        },
        {
            "comparison": "engineering_target_low",
            "numerator": "task brief",
            "denominator": "task brief",
            "ratio": 1.5,
            "delta_delta_g_kj_mol": "",
            "reference_or_role": "pre-declared interpretation scale",
            "formula": "target interval lower edge",
        },
        {
            "comparison": "engineering_target_high",
            "numerator": "task brief",
            "denominator": "task brief",
            "ratio": 2.3,
            "delta_delta_g_kj_mol": "",
            "reference_or_role": "pre-declared interpretation scale",
            "formula": "target interval upper edge",
        },
    ]
    path = output_dir / "20261109_量级基线.csv"
    _write_csv(path, rows)
    return path


def _write_audit_table(output_dir: Path) -> Path:
    path = output_dir / "20261109_计算状态审计.csv"
    _write_csv(path, build_audit_rows())
    return path


def _build_reproduction_input(supplied_csv: Path, output_dir: Path) -> Path:
    supplied = _read_csv(supplied_csv)
    rows = []
    for row in supplied:
        rows.append(
            {
                "draw": row["draw"],
                "h_group_nm": row["h_group_nm"],
                "pH_group": row["nominal_pH"],
                "pH": row["actual_pH"],
                "conditional_I_mM": row["conditional_I_mM"],
                "zeta_mV": row["measured_zeta_mV"],
            }
        )
    path = output_dir / "20261109_face_charge_reproduction_input.csv"
    _write_csv(path, rows)
    return path


def _compare_csvs(reference: Path, reproduced: Path, output_dir: Path) -> Path:
    reference_rows = _read_csv(reference)
    reproduced_rows = _read_csv(reproduced)
    if len(reference_rows) != len(reproduced_rows):
        raise ValueError(f"Row count differs: {len(reference_rows)} != {len(reproduced_rows)}")
    columns = [key for key in reference_rows[0] if key in reproduced_rows[0]]
    reports: list[dict[str, object]] = []
    max_all = 0.0
    for column in columns:
        max_diff = 0.0
        numeric = True
        for left, right in zip(reference_rows, reproduced_rows):
            if left[column] == right[column]:
                continue
            try:
                diff = abs(float(left[column]) - float(right[column]))
            except ValueError:
                numeric = False
                break
            max_diff = max(max_diff, diff)
        if numeric:
            max_all = max(max_all, max_diff)
        reports.append(
            {
                "column": column,
                "numeric": numeric,
                "max_absolute_difference": max_diff if numeric else "",
                "status": "PASS" if (numeric and max_diff <= 1e-10) or (not numeric and all(a[column] == b[column] for a, b in zip(reference_rows, reproduced_rows))) else "FAIL",
            }
        )
    if any(row["status"] == "FAIL" for row in reports):
        raise AssertionError("face_charge reproduction differs from supplied CSV")
    path = output_dir / "20261109_face_charge_reproduction_check.csv"
    _write_csv(path, reports)
    (output_dir / "20261109_face_charge_reproduction_summary.txt").write_text(
        f"rows={len(reference_rows)} columns={len(columns)} max_absolute_difference={max_all:.17g} status=PASS\n",
        encoding="utf-8",
    )
    return path


def _render_figures(output_dir: Path) -> tuple[Path, Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    plt.rcParams.update({"font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans"], "axes.unicode_minus": False})
    groups = ["182", "357", "604"]
    rows = structure_rows()
    thickness = [row["thickness_nm"] for row in rows]
    dsc = [100 * row["dsc_crystallinity_mass_frac"] for row in rows]
    xrd = [100 * row["xrd_crystallinity_frac"] for row in rows]
    coherent = [row["coherent_110_nm"] for row in rows]
    diameter = [row["fiber_diameter_nm"] for row in rows]
    zeta = [row["zeta_mV"] for row in rows]
    density = [row["site_density_per_nm2"] for row in rows]

    fig, axes = plt.subplots(2, 2, figsize=(11, 7.2), constrained_layout=True)
    ax = axes[0, 0]
    x = np.arange(3)
    ax.bar(x - 0.18, dsc, 0.36, label="DSC", color="#4C78A8")
    ax.bar(x + 0.18, xrd, 0.36, label="XRD", color="#F58518")
    ax.set_xticks(x, groups)
    ax.set_ylabel("结晶度 (%)")
    ax.set_title("体相结构约束（输入）")
    ax.legend(frameon=False)

    ax = axes[0, 1]
    ax.plot(groups, thickness, "o-", label="膜厚 (nm)", color="#54A24B")
    ax.plot(groups, coherent, "s--", label="XRD(110)相干尺寸 (nm)", color="#B279A2")
    ax.plot(groups, diameter, "^:", label="SEM纤维直径 (nm)", color="#E45756")
    ax.set_ylabel("nm")
    ax.set_title("尺度输入；不等于原子表面")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 0]
    ax2 = ax.twinx()
    ax.plot(groups, zeta, "o-", color="#4C78A8", label="ζ (mV)")
    ax2.plot(groups, density, "s--", color="#E45756", label="有效位点密度 (nm^-2)")
    ax.set_ylabel("ζ (mV)")
    ax2.set_ylabel("条件反演位点密度 (nm^-2)")
    ax.set_title("实测/条件反演对照；不用于调表面模型")
    ax.grid(axis="y", alpha=0.2)

    ax = axes[1, 1]
    ax.axis("off")
    ax.text(
        0.02,
        0.98,
        "本轮模拟边界\n\n"
        "[已] 三组宏观约束表已冻结\n"
        "[已] 荷电基线与 RT ln(R) 已复算\n"
        "[缺] 三组×3 原子表面结构\n"
        "[缺] 辛酸拓扑/PE 力场/GROMACS\n"
        "[未产生] W(z)、K干、R干\n"
        "[未启动] 湿态 HA/A-（等待干态方向）",
        va="top",
        fontsize=11,
        bbox={"boxstyle": "round", "facecolor": "#F2F2F2", "edgecolor": "#888888"},
    )
    fig.suptitle("图1 结构约束、实测对照与模拟边界（2026-11-09）", fontsize=14)
    png1 = output_dir / "图1_结构约束与模拟边界.png"
    svg1 = output_dir / "图1_结构约束与模拟边界.svg"
    fig.savefig(png1, dpi=300)
    fig.savefig(svg1)
    plt.close(fig)

    values = compute_baseline()
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.4), constrained_layout=True)
    ax = axes[0]
    labels = ["薄/厚\n位点比", "中/厚\n位点比", "薄/中\n位点比", "非晶相\n简单上限"]
    ratios = [values["thin_over_thick"], values["middle_over_thick"], values["thin_over_middle"], values["noncrystalline_bound"]]
    colors = ["#E45756", "#72B7B2", "#F2CF5B", "#B0B0B0"]
    ax.bar(labels, ratios, color=colors)
    ax.axhline(1.0, color="black", lw=1, label="R=1")
    ax.axhline(1.18, color="#888888", ls="--", lw=1, label="1.18")
    ax.axhline(1.84, color="#4C78A8", ls=":", lw=1.5, label="1.84")
    ax.set_ylabel("比值 R")
    ax.set_ylim(0, 2.1)
    ax.set_title("预计算量级标尺")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1]
    delta = [values["thin_over_thick_delta_delta_g_kj_mol"], values["middle_over_thick_delta_delta_g_kj_mol"], values["thin_over_middle_delta_delta_g_kj_mol"]]
    ax.bar(["薄/厚", "中/厚", "薄/中"], delta, color=["#E45756", "#72B7B2", "#F2CF5B"])
    ax.axhline(0.0, color="black", lw=1)
    ax.set_ylabel("−RT ln(R) (kJ/mol)")
    ax.set_title("相同活度/稀覆盖下的等效尺度")

    ax = axes[2]
    ax.axis("off")
    ax.text(
        0.02,
        0.98,
        "判读规则（预声明）\n\n"
        "• 1.5–2.3：工程解释尺度\n"
        "• 1.84：水中ζ条件反演参照\n"
        "• 1.29：固定两相性质的简单上限\n\n"
        "当前没有吸附曲线，不能给 R干/R湿。\n"
        "只有获得 W(z) 并完成面积归一后，\n"
        "才能把该标尺与分子模拟比较。",
        va="top",
        fontsize=11,
        bbox={"boxstyle": "round", "facecolor": "#F2F2F2", "edgecolor": "#888888"},
    )
    fig.suptitle("图2 量级标尺与交付状态（预计算；非吸附自由能曲线）", fontsize=14)
    png2 = output_dir / "图2_量级标尺与交付状态.png"
    svg2 = output_dir / "图2_量级标尺与交付状态.svg"
    fig.savefig(png2, dpi=300)
    fig.savefig(svg2)
    plt.close(fig)
    return png1, png2


def write_summary(output_dir: Path, supplied_csv: Path) -> Path:
    values = compute_baseline()
    package_versions = []
    for package_name in ("numpy", "scipy", "matplotlib", "pandas"):
        try:
            version = importlib.metadata.version(package_name)
        except importlib.metadata.PackageNotFoundError:
            version = "MISSING"
        package_versions.append(f"{package_name}={version}")
    command_status = ", ".join(
        f"{name}={'FOUND' if shutil.which(name) else 'MISSING'}"
        for name in ("gmx", "wham", "plumed", "obabel", "packmol")
    )
    text = f"""# 2026-11-09 表面吸附任务：可计算基线与停止边界

本记录只写本机实际复算和资料审计，不把计划中的分子模拟写成结果。

## 已完成

- 三组加工状态的宏观结构约束表：182/357/604，对应第二次双向拉伸倍率 1.5/2.1/2.73。
- 由随附九组 CSV 的 pH≈7 条件反演位点密度重算：薄/厚 = {values['thin_over_thick']:.12f}，中/厚 = {values['middle_over_thick']:.12f}，薄/中 = {values['thin_over_middle']:.12f}。
- 在 298.15 K、共同稀气相活度、稀覆盖和相同最大位点密度假设下，−RT ln(R) 分别为 {values['thin_over_thick_delta_delta_g_kj_mol']:.6f}、{values['middle_over_thick_delta_delta_g_kj_mol']:.6f}、{values['thin_over_middle_delta_delta_g_kj_mol']:.6f} kJ/mol。
- 固定两相性质的简单非晶相上限为 {values['noncrystalline_bound']:.12f}。
- 复用已跟踪的 `WP3/v9_几何层表.csv` 计算几何代理：182/357/604 的层数为 11/6/4，外层周长占比为 0.1818/0.3333/0.5000；这只是结构情景，不是原子表面或吸附容量。
- 将随附九组输出整理为 face_charge.py 的最小输入，并逐列复算，结果见 `20261109_face_charge_reproduction_summary.txt`。

## 尚未产生的量

- 没有原子坐标、半结晶 PE 表面、辛酸/辛酸根拓扑或兼容力场输入。
- 当前环境没有 GROMACS、WHAM/MBAR、MDTraj、ParmEd 或 OpenMM；因此没有生成 W(z)、K干、R干、R湿 或任何模拟误差条。
- 湿态 HA/A⁻ 按任务书依赖干态方向，当前未启动。

## 科学边界

`1.84` 是水中 ζ 条件反演的有效位点密度比，`1.29` 是固定两相吸附性质的简单混合上限；它们都不是辛酸吸附自由能。表面能、AFM 粗糙度、投影面积和单位质量比表面积没有被当作真实外表面位点密度。当前可交付的是“结构情景与目标幅度的可行性准备”，不能写成“薄组已由分子模拟证明更易吸附”。

输入来源：`{supplied_csv}`；脚本：`src/mechanism_surface_adsorption_baseline.py`；生成时间由文件系统记录。

## 本机环境

- Python：`{sys.version.split()[0]}`
- Python 包：`{', '.join(package_versions)}`
- 分子模拟工具：`{command_status}`
"""
    path = output_dir / "20261109_复算记录与科研边界.md"
    path.write_text(text, encoding="utf-8")
    return path


def run_reproduction(output_dir: Path, supplied_csv: Path, face_charge: Path) -> Path:
    input_csv = _build_reproduction_input(supplied_csv, output_dir)
    reproduced_csv = output_dir / "20261109_face_charge_reproduced.csv"
    result = subprocess.run(
        [sys.executable, str(face_charge), str(input_csv), str(reproduced_csv)],
        check=True,
        capture_output=True,
        text=True,
    )
    (output_dir / "20261109_face_charge_stdout.txt").write_text(result.stdout, encoding="utf-8")
    return _compare_csvs(supplied_csv, reproduced_csv, output_dir)


def write_manifest(output_dir: Path) -> Path:
    paths = sorted(path for path in output_dir.rglob("*") if path.is_file())
    rows = []
    for path in paths:
        rel = path.relative_to(output_dir).as_posix()
        if rel == "manifest.csv":
            continue
        role = "figure" if rel.startswith("图") or rel.startswith("figures/") else "processed_output"
        rows.append({"relative_path": rel, "role": role, "status": "GENERATED_OR_RECORDED"})
    path = output_dir / "manifest.csv"
    _write_csv(path, rows)
    return path


def build_outputs(
    output_dir: Path = DEFAULT_OUTPUT,
    supplied_csv: Path = SUPPLIED_CSV,
    face_charge: Path = FACE_CHARGE,
) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    if not supplied_csv.is_file():
        raise FileNotFoundError(supplied_csv)
    if not face_charge.is_file():
        raise FileNotFoundError(face_charge)
    _write_structure_table(output_dir)
    _write_baseline_table(output_dir)
    _write_audit_table(output_dir)
    run_reproduction(output_dir, supplied_csv, face_charge)
    _render_figures(output_dir)
    write_summary(output_dir, supplied_csv)
    write_manifest(output_dir)
    return output_dir


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--supplied-csv", type=Path, default=SUPPLIED_CSV)
    parser.add_argument("--face-charge", type=Path, default=FACE_CHARGE)
    args = parser.parse_args()
    output = build_outputs(args.output_dir, args.supplied_csv, args.face_charge)
    values = compute_baseline()
    print(f"output={output}")
    print(f"thin_over_thick={values['thin_over_thick']:.12f}")
    print(f"middle_over_thick={values['middle_over_thick']:.12f}")
    print(f"thin_over_middle={values['thin_over_middle']:.12f}")
    print(f"noncrystalline_bound={values['noncrystalline_bound']:.12f}")
    print("face_charge_reproduction=PASS")
    print("md_status=BLOCKED_MISSING_MD_INPUTS")


if __name__ == "__main__":
    main()
