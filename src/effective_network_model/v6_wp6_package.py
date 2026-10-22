"""Bundle the frozen v6 study with a path/size inventory (no new hashes)."""

from __future__ import annotations

import csv
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "WP6" / "reproduction"
MENTOR = Path("原始资料/赵珈绮/12_论文导向研究总方案_v6.0.0")
FIXED = [
    Path("项目进展计划与任务看板.md"),
    Path("阶段总结/20261021_v6_WP2数值验收.md"),
    Path("阶段总结/20261021_v6_WP3机制分工与M2资格.md"),
    Path("阶段总结/20261021_v6_WP4全曲线比较.md"),
    Path("阶段总结/20261022_v6_WP5预测辨识与稳健性.md"),
    Path("阶段总结/20261022_v6_WP6论文证据与交接.md"),
    Path("data/processed/salt_transport/salt_transport_tidy.csv"),
    Path("data/processed/origin_exports/Polymer_salt_2h/Book1__Sheet1.csv"),
    Path("原始资料/赵珈绮/01_研究框架与资料_v2/05_真实数据_v1核查/source_snapshots/Polymer_salt_2h.opju"),
    Path("src/data_ingest/build_wp1_manifest.py"),
]


def role(path: Path) -> str:
    parts = path.parts
    if parts[0] == "原始资料":
        return "read-only source or mentor specification"
    if parts[0] == "data":
        return "audited upstream input"
    if parts[0] == "src":
        return "versioned executable source"
    if parts[0] == "阶段总结" or path.name == "项目进展计划与任务看板.md":
        return "formal status or scientific stage conclusion"
    if parts[0] == "WP6":
        return "manuscript, figure or reproduction guide"
    return "frozen input or recorded output"


def main() -> None:
    files = set(FIXED)
    for root_name in ("WP1", "WP2", "WP3", "WP4", "WP5", "WP6"):
        for path in (ROOT / root_name).rglob("*"):
            if (path.is_file() and "__pycache__" not in path.parts and
                    path.name not in ("fontlist-v3.11.0.json", "v6_reproduction.zip",
                                      "inventory.csv")):
                files.add(path.relative_to(ROOT))
    for path in (ROOT / "src" / "effective_network_model").glob("v6_*.py"):
        files.add(path.relative_to(ROOT))
    for path in (ROOT / MENTOR).glob("0[0-5]_*.md"):
        files.add(path.relative_to(ROOT))
    missing = [str(path) for path in files if not (ROOT / path).is_file()]
    if missing:
        raise FileNotFoundError(missing)
    ordered = sorted(files, key=lambda p: p.as_posix())
    OUT.mkdir(exist_ok=True)
    with (OUT / "inventory.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("project_relative_path", "role", "bytes"))
        writer.writerows((p.as_posix(), role(p), (ROOT / p).stat().st_size) for p in ordered)
    with ZipFile(OUT / "v6_reproduction.zip", "w", compression=ZIP_DEFLATED) as archive:
        for path in ordered:
            archive.write(ROOT / path, arcname=path.as_posix())
        archive.write(OUT / "inventory.csv", arcname="WP6/reproduction/inventory.csv")
    print(f"{len(ordered)} files; archive bytes={(OUT / 'v6_reproduction.zip').stat().st_size}")


if __name__ == "__main__":
    main()
