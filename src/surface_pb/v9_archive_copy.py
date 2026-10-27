"""Copy the approved v9 evidence snapshot into the new 19 delivery folder."""

import csv
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DESTINATION = ROOT / "原始资料/赵珈绮/19_学生薄化荷电成果_v9.0.0"
MANIFEST = ROOT / "WP6/v9_拟归档清单.csv"


def main():
    with MANIFEST.open(encoding="utf-8-sig", newline="") as stream:
        records = list(csv.DictReader(stream))
    if not records:
        raise AssertionError("empty reviewable archive manifest")
    if DESTINATION.exists():
        raise FileExistsError(f"refusing to overwrite existing delivery: {DESTINATION}")
    root = ROOT.resolve()
    destination = DESTINATION.resolve()
    pairs = []
    targets = set()
    for record in records:
        source = (ROOT / record["source_relative"]).resolve()
        target = (DESTINATION / record["archive_relative"]).resolve()
        if not source.is_relative_to(root) or not source.is_file():
            raise AssertionError(f"source outside project or missing: {source}")
        if not target.is_relative_to(destination) or target in targets:
            raise AssertionError(f"unsafe or duplicate archive target: {target}")
        targets.add(target)
        pairs.append((source, target))
    DESTINATION.mkdir(parents=False, exist_ok=False)
    copied = 0
    for source, target in pairs:
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        copied += 1
    print(f"copied {copied} files to {DESTINATION}")


if __name__ == "__main__":
    main()
