"""Recalculate the four B1 I1/H1 critical boundaries from frozen Book2 input.

This is an erratum overlay. It does not edit the archived 20 package.
"""

from __future__ import annotations

import csv
from decimal import Decimal as D, localcontext
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
H18 = ROOT / "原始资料/赵珈绮/18_超薄PE膜荷电与阻盐研究资料/18H_Stage2_B1_学生任务书_20261027_v1.0.0"
P20 = ROOT / "原始资料/赵珈绮/20_学生Stage2_B1成果_v1.0.0"
OUT = ROOT / "data/processed/salt_transport/B1_I1_H1_临界勘误_20261028.csv"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def concentration(k: D) -> D:
    if k == 0:
        return D(0)
    if k < 10:
        return k * concentration(D(10)) / 10
    z = k.log10()
    return D(1000) / D("58.44") * D(10) ** (D(".0165") * z * z + D(".9398") * z - D("3.2793"))


def read_points() -> dict[tuple[str, D], D]:
    raw = rows(H18 / "inputs/origin/Book2__Sheet1.csv")
    mapping = rows(H18 / "inputs/曲线列映射.csv")
    points = {}
    for item in mapping:
        for line in raw:
            value = line[item["conductivity_column"]]
            if value:
                points[item["sample_id"], D(line[item["time_column"]])] = D(value)
    return points


def resistance(points: dict[tuple[str, D], D], sample: str, conversion: str) -> D:
    k0 = points[sample, D(0)]

    def increment(t: int) -> D:
        k = points[sample, D(t)]
        return concentration(k - k0) if conversion == "C1" else concentration(k) - concentration(k0)

    theta0 = D(1) - increment(0) / D(300)
    theta1 = D(1) - increment(20) / D(300)
    transfer = (theta0.ln() - theta1.ln()) / (D(20) * D(60))
    return 1 / transfer


def feasible(x: list[D], y: list[D], delta: D) -> bool:
    lower, upper = D(1), D("Infinity")
    for i in range(3):
        for k in range(3):
            coefficient = (1 - delta) * x[k] - (1 + delta) * x[i]
            rhs = (1 + delta) * y[k] - (1 - delta) * y[i]
            if coefficient > 0:
                upper = min(upper, rhs / coefficient)
            elif coefficient < 0:
                lower = max(lower, rhs / coefficient)
            elif rhs < 0:
                return False
    return upper >= lower


def boundary(x: list[D], y: list[D]) -> D:
    lower, upper = D(0), D(".999999")
    assert not feasible(x, y, lower) and feasible(x, y, upper)
    for _ in range(80):
        middle = (lower + upper) / 2
        if feasible(x, y, middle):
            upper = middle
        else:
            lower = middle
    return (lower + upper) / 2


def main() -> None:
    original = rows(P20 / "WP2/critical_tolerance.csv")
    submitted = {(r["scenario"], r["scope"], r["model"]): r for r in original}
    points = read_points()
    corrected = []
    with localcontext() as context:
        context.prec = 55
        for conversion in ("C1", "C2"):
            x = [resistance(points, f"R{h}", conversion) for h in (65, 700, 2600)]
            y = [resistance(points, f"B{h}", conversion) for h in (65, 700, 2600)]
            value = boundary(x, y)
            for timing in ("T85", "T80"):
                scenario = f"{conversion}_{timing}"
                old = submitted[scenario, "I1", "H1"]
                assert not D(old["delta_lower"]) <= value <= D(old["delta_upper"])
                corrected.append({"scenario": scenario, "scope": "I1", "model": "H1",
                                  "independent_boundary": str(value),
                                  "submitted_lower": old["delta_lower"],
                                  "submitted_upper": old["delta_upper"],
                                  "original_bracket_contains_boundary": "False",
                                  "boundary_type": "infimum; finite q not identified"})
    if not OUT.parent.is_dir():
        raise FileNotFoundError(OUT.parent)
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(corrected[0]))
        writer.writeheader()
        writer.writerows(corrected)
    print(f"Wrote {len(corrected)} corrected rows: {OUT}")


if __name__ == "__main__":
    main()
