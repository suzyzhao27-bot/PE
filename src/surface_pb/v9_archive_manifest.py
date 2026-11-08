"""Prepare reviewable 19-package copy plan without creating the directory."""

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / "原始资料/赵珈绮/18_薄化荷电主线任务书_v9.0.0"
DESTINATION = ROOT / "原始资料/赵珈绮/19_学生薄化荷电成果_v9.0.0"
OUT = ROOT / "WP6/v9_拟归档清单.csv"


def add(items, source, target, role):
    if not source.is_file():
        raise FileNotFoundError(source)
    items.append({"source_relative": source.relative_to(ROOT).as_posix(),
                  "archive_relative": target.as_posix(),
                  "role": role, "bytes": source.stat().st_size,
                  "status": "PREPARED_NOT_COPIED"})


def main():
    items = []
    for source in sorted(PACKAGE.rglob("*")):
        if source.is_file():
            add(items, source, Path("taskbook") / source.relative_to(PACKAGE),
                "FROZEN_ORIGINAL_18_COPY")
    for source in sorted((ROOT / "src/data_ingest").glob("v9_*.py")):
        add(items, source, Path("src/data_ingest") / source.name, "V9_CODE")
    for source in sorted((ROOT / "src/surface_pb").glob("v9_*.py")):
        add(items, source, Path("src/surface_pb") / source.name, "V9_CODE")
    old = ROOT / "src/surface_pb/fibre_pb.py"
    add(items, old, Path("src/surface_pb/fibre_pb.py"), "LEGACY_TECHNICAL_REFERENCE")
    for folder in ("WP1", "WP2", "WP3", "WP4", "WP5", "WP6"):
        for source in sorted((ROOT / folder).glob("v9_*")):
            if source.is_file() and source != OUT and source.name != "v9_拟归档README.md":
                add(items, source, Path(folder) / source.name, "V9_EVIDENCE")
    for source in sorted((ROOT / "阶段总结").glob("*v9*")):
        if source.is_file():
            add(items, source, Path("阶段总结") / source.name, "STAGE_CONCLUSION")
    for source in sorted((ROOT / "data/processed/figures").glob("v9_*")):
        if source.is_file():
            add(items, source, Path("data/processed/figures") / source.name, "V9_FIGURE")
    add(items, ROOT / "复现记录_v9_薄化荷电.md", Path("复现记录_v9_薄化荷电.md"),
        "REPRODUCTION_RECORD")
    add(items, ROOT / "WP6/v9_拟归档README.md", Path("README.md"), "ARCHIVE_README")
    add(items, ROOT / ".planning/2026-10-26-薄化荷电主线-v9/progress.md",
        Path("checkpoint_progress.md"), "CHECKPOINT_SNAPSHOT")
    items.append({"source_relative": "WP6/v9_拟归档清单.csv",
                  "archive_relative": "WP6/v9_拟归档清单.csv",
                  "role": "ARCHIVE_MANIFEST_SELF", "bytes": "",
                  "status": "PREPARED_NOT_COPIED"})
    if DESTINATION.is_dir():
        replay = DESTINATION / "WP6/v9_包内回放记录.json"
        status = "COPIED_AND_REPLAYED" if replay.is_file() else "COPIED_PENDING_REPLAY"
        for item in items:
            item["status"] = status
    targets = [r["archive_relative"] for r in items]
    if len(targets) != len(set(targets)):
        raise AssertionError("duplicate archive destination")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(items[0]))
        writer.writeheader()
        writer.writerows(items)
    print(f"listed {len(items)} files, {sum(r['bytes'] for r in items if isinstance(r['bytes'], int))} source bytes plus manifest; status {items[0]['status']}")


if __name__ == "__main__":
    main()
