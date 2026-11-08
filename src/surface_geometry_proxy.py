"""Summarise the project's existing square-fiber geometry as a non-atomic proxy.

This is a structural preparation step for the adsorption task.  It consumes the
already tracked v9 geometry table and reports layer count, perimeter/projection
ratio, outer-layer share, and a 1 nm cross-sectional shell fraction.  It never
computes an atomistic surface, a force field, W(z), K, or an adsorption ratio.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
GEOMETRY_TABLE = ROOT / "WP3" / "v9_几何层表.csv"
OUTPUT_DIR = ROOT / "data" / "processed" / "surface_adsorption"
SHELL_NM = 1.0


def _read_rows(path: Path = GEOMETRY_TABLE) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def _f(row: dict[str, str], key: str) -> float:
    return float(row[key])


def build_proxy_rows(path: Path = GEOMETRY_TABLE) -> list[dict[str, object]]:
    source_rows = _read_rows(path)
    output: list[dict[str, object]] = []
    for draw in ("182", "357", "604"):
        rows = [row for row in source_rows if row["case_id"] == f"family_{draw}_free"]
        if not rows:
            raise ValueError(f"Missing family_{draw}_free")
        first = rows[0]
        d_nm = _f(first, "d_nm")
        n = int(first["N"])
        perimeter_total = sum(_f(row, "wet_perimeter_nm") for row in rows)
        outer_perimeter = sum(_f(row, "wet_perimeter_nm") for row in rows if row["is_outer"] == "True")
        one_nm_shell = 1.0 - max(d_nm - 2 * SHELL_NM, 0.0) ** 2 / d_nm**2
        output.append(
            {
                "draw": draw,
                "thickness_nm": _f(first, "h_nm"),
                "N": n,
                "fiber_width_nm": d_nm,
                "period_nm": _f(first, "period_nm"),
                "wet_area_m2_per_m": _f(first, "wet_area_m2_per_m"),
                "projected_area_m2_per_m": _f(first, "projected_area_m2_per_m"),
                "wet_to_projected_area_ratio": _f(first, "wet_area_m2_per_m") / _f(first, "projected_area_m2_per_m"),
                "outer_layer_perimeter_fraction": outer_perimeter / perimeter_total,
                "one_nm_shell_fraction": one_nm_shell,
                "shell_thickness_nm": SHELL_NM,
                "evidence_role": "legacy_geometry_output",
                "model_status": "GEOMETRY_PROXY_NOT_ATOMIC_SURFACE",
                "adsorption_status": "NOT_COMPUTED",
            }
        )
    return output


def build_backed_rows(path: Path = GEOMETRY_TABLE) -> list[dict[str, object]]:
    source_rows = _read_rows(path)
    output: list[dict[str, object]] = []
    for draw in ("182", "357", "604"):
        rows = [row for row in source_rows if row["case_id"] == f"family_{draw}_backed"]
        if not rows:
            raise ValueError(f"Missing family_{draw}_backed")
        first = rows[0]
        total = sum(_f(row, "wet_perimeter_nm") for row in rows)
        outer = sum(_f(row, "wet_perimeter_nm") for row in rows if row["is_outer"] == "True")
        output.append(
            {
                "draw": draw,
                "N": int(first["N"]),
                "wet_area_m2_per_m": _f(first, "wet_area_m2_per_m"),
                "projected_area_m2_per_m": _f(first, "projected_area_m2_per_m"),
                "wet_to_projected_area_ratio": _f(first, "wet_area_m2_per_m") / _f(first, "projected_area_m2_per_m"),
                "outer_layer_perimeter_fraction": outer / total,
                "evidence_role": "legacy_geometry_output",
                "model_status": "GEOMETRY_PROXY_NOT_ATOMIC_SURFACE",
            }
        )
    return output


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_outputs(output_dir: Path = OUTPUT_DIR, geometry_table: Path = GEOMETRY_TABLE) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    free = build_proxy_rows(geometry_table)
    backed = build_backed_rows(geometry_table)
    _write_csv(output_dir / "20261109_几何外层代理_自由膜.csv", free)
    _write_csv(output_dir / "20261109_几何外层代理_背衬膜.csv", backed)
    contract = {
        "status": "GEOMETRY_PROXY_ONLY",
        "source_table": str(geometry_table),
        "shell_thickness_nm": SHELL_NM,
        "groups": [row["draw"] for row in free],
        "derived_quantities": [
            "layer_count_N",
            "wet_to_projected_area_ratio",
            "outer_layer_perimeter_fraction",
            "one_nm_shell_fraction",
        ],
        "not_computed": ["atomic_accessible_area", "W(z)", "K_dry", "R_dry", "R_wet"],
        "interpretation": "structural scenario proxy; not an atomistic surface or adsorption prediction",
    }
    (output_dir / "20261109_几何代理合同.json").write_text(json.dumps(contract, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    _write_markdown(output_dir / "20261109_几何代理说明.md", free, backed)
    _render_figure(output_dir, free)
    return output_dir


def _write_markdown(path: Path, free: list[dict[str, object]], backed: list[dict[str, object]]) -> None:
    lines = [
        "# 2026-11-09 现有纤维几何的外层代理",
        "",
        "本记录使用仓库已跟踪的 `WP3/v9_几何层表.csv`，只做结构接口复用。它不是三组原子表面，也不产生辛酸吸附能。",
        "",
        "| 组 | 层数 N | 纤维宽度 (nm) | 湿/投影面积比 | 外层周长占比 | 1 nm 截面壳层比例 |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in free:
        lines.append(
            f"| {row['draw']} | {row['N']} | {row['fiber_width_nm']:.5f} | "
            f"{row['wet_to_projected_area_ratio']:.5f} | {row['outer_layer_perimeter_fraction']:.5f} | "
            f"{row['one_nm_shell_fraction']:.5f} |"
        )
    lines += [
        "",
        "## 可以推进到哪里",
        "",
        "几何代理显示三组的层数和外层占比仍有可区分的结构情景：182/357/604 的层数为 11/6/4，自由膜外层周长占比分别为 0.1818/0.3333/0.5000。1 nm 截面壳层比例由纤维宽度得到，分别约为 0.1563/0.1695/0.1448。这个结果只说明宏观代理结构没有在输入阶段自动塌缩为同一模型。",
        "",
        "## 不能推出什么",
        "",
        "湿/投影面积比、外层周长占比和 1 nm 壳层比例都不是真实 PE 外表面原子面积、可达面积或单位面积吸附容量。它们不替代表面坐标、力场、辛酸拓扑和自由能采样，不能给出 `W(z)`、`K_dry`、`R_dry` 或 `R_wet`。",
        "",
        "## 背衬代理",
        "",
        "背衬几何的对应表保存在 `20261109_几何外层代理_背衬膜.csv`，底层不湿面的处理沿用已有 v9 几何接口。",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def _render_figure(output_dir: Path, rows: list[dict[str, object]]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans"],
        "axes.unicode_minus": False,
    })
    groups = [row["draw"] for row in rows]
    n_layers = [row["N"] for row in rows]
    outer = [row["outer_layer_perimeter_fraction"] for row in rows]
    shell = [row["one_nm_shell_fraction"] for row in rows]
    ratio = [row["wet_to_projected_area_ratio"] for row in rows]
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 4.2), constrained_layout=True)
    axes[0].bar(groups, n_layers, color="#4C78A8")
    axes[0].set_title("层数 N（已有几何接口）")
    axes[0].set_ylabel("层数")
    axes[1].plot(groups, outer, "o-", label="外层周长占比", color="#E45756")
    axes[1].plot(groups, shell, "s--", label="1 nm 截面壳层比例", color="#54A24B")
    axes[1].set_ylim(0, 0.6)
    axes[1].set_ylabel("比例")
    axes[1].set_title("外层代理；不是原子表面")
    axes[1].legend(frameon=False, fontsize=8)
    axes[2].bar(groups, ratio, color="#F2CF5B")
    axes[2].set_title("湿/投影面积比")
    axes[2].set_ylabel("几何比")
    fig.suptitle("图3 现有纤维几何的外层代理（结构情景，不是吸附结果）", fontsize=13)
    fig.savefig(output_dir / "图3_几何外层代理.png", dpi=300)
    fig.savefig(output_dir / "图3_几何外层代理.svg")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometry-table", type=Path, default=GEOMETRY_TABLE)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    args = parser.parse_args()
    output = write_outputs(args.output_dir, args.geometry_table)
    print(f"output={output}")
    print("status=GEOMETRY_PROXY_ONLY")


if __name__ == "__main__":
    main()
