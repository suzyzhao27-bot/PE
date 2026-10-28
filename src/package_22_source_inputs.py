"""Package frozen 22-series source inputs for local, unsent handoff."""

from __future__ import annotations

from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "原始资料/赵珈绮/22_机制预测与实验验证方案_v1.0.0"
BOOK2 = ROOT / "原始资料/赵珈绮/20_学生Stage2_B1成果_v1.0.0/WP1/book2_long.csv"
OUT = ROOT / "原始资料/赵珈绮/03_发送用ZIP/Stage2_22号源输入交接包_20261028_v1.0.0.zip"
NOTE = """22号源输入本地交接包（未送达）

本包仅含22号原始任务包和20号Book2长表，保留相对项目根目录路径。
不含本机新计算结果；代码、结果及封存JSON/PDF见Git的dev-zhao设备分支。
22号包内99_入口更新前备份属于历史快照，不覆盖项目根AGENTS.md与正式看板。
解包后先按freeze_predictions.json的source_input_manifest核对所用文件。
学生本人接收、运行和导师验收均未发生。
"""


def main() -> None:
    if not SOURCE.is_dir() or not BOOK2.is_file() or not OUT.parent.is_dir():
        raise FileNotFoundError("Source, Book2 or existing handoff directory missing")
    if OUT.exists():
        raise FileExistsError(f"Do not overwrite existing handoff: {OUT}")
    files = sorted(path for path in SOURCE.rglob("*") if path.is_file()) + [BOOK2]
    if len(files) != 82:
        raise ArithmeticError(f"Expected 81 package files plus Book2; got {len(files)}")
    expected = {path.relative_to(ROOT).as_posix() for path in files}
    with ZipFile(OUT, "x", compression=ZIP_DEFLATED, compresslevel=7) as archive:
        for path in files:
            archive.write(path, arcname=path.relative_to(ROOT).as_posix())
        archive.writestr("README_交接范围.txt", NOTE)
    with ZipFile(OUT) as archive:
        names = archive.namelist()
        if set(names) != expected | {"README_交接范围.txt"} or len(names) != 83:
            raise ArithmeticError("Archive member mismatch")
        if archive.testzip() is not None:
            raise ArithmeticError("Archive CRC failure")
    print(f"zip={OUT.name} source_files={len(files)} entries=83 bytes={OUT.stat().st_size}")


if __name__ == "__main__":
    main()
