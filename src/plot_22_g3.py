"""Render G3 Book2 salt observations and conditional urea targets."""

from __future__ import annotations

import csv
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
BOOK = ROOT / "原始资料/赵珈绮/20_学生Stage2_B1成果_v1.0.0/WP1/book2_long.csv"
UREA = ROOT / "data/processed/salt_transport/22_legacy_urea_conditions_20261028.csv"
OUT = ROOT / "data/processed/salt_transport"
FONT = Path("C:/Windows/Fonts/msyh.ttc")
W, H = 1400, 980
INK = "#243144"
GRID = "#d7dde6"
BLUE = "#2b6cb0"
ORANGE = "#df6c28"
BG = "#ffffff"


def read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def font(size: int) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(str(FONT), size)


def plot_panel(draw: ImageDraw.ImageDraw, area: tuple[int, int, int, int],
               title: str, ylabel: str, ymax: float, yticks: list[float],
               series: list[tuple[str, list[tuple[float, float]], str]],
               bands: list[tuple[list[tuple[float, float]], list[tuple[float, float]], str]] | None = None) -> None:
    left, top, right, bottom = area
    plot_left, plot_right = left + 120, right - 45
    plot_top, plot_bottom = top + 62, bottom - 62
    draw.text((left, top), title, fill=INK, font=font(28))
    draw.text((left, top + 36), ylabel, fill=INK, font=font(19))
    for tick in yticks:
        y = plot_bottom - (tick / ymax) * (plot_bottom - plot_top)
        draw.line((plot_left, y, plot_right, y), fill=GRID, width=2)
        draw.text((left + 20, y - 13), f"{tick:g}", fill=INK, font=font(17))
    for t in (0, 20, 40, 60, 80, 100, 120):
        x = plot_left + t / 120 * (plot_right - plot_left)
        draw.line((x, plot_top, x, plot_bottom), fill=GRID, width=1)
        draw.text((x - 14, plot_bottom + 9), str(t), fill=INK, font=font(17))
    draw.line((plot_left, plot_top, plot_left, plot_bottom), fill=INK, width=3)
    draw.line((plot_left, plot_bottom, plot_right, plot_bottom), fill=INK, width=3)
    draw.text((plot_right - 112, plot_bottom + 30), "时间 / min", fill=INK, font=font(18))
    for low, high, fill in bands or []:
        polygon = ([(plot_left + t / 120 * (plot_right - plot_left),
                     plot_bottom - value / ymax * (plot_bottom - plot_top)) for t, value in low]
                   + [(plot_left + t / 120 * (plot_right - plot_left),
                       plot_bottom - value / ymax * (plot_bottom - plot_top)) for t, value in reversed(high)])
        draw.polygon(polygon, fill=fill)
    for index, (label, points, color) in enumerate(series):
        coords = [(plot_left + t / 120 * (plot_right - plot_left),
                   plot_bottom - value / ymax * (plot_bottom - plot_top)) for t, value in points]
        if len(coords) > 1:
            draw.line(coords, fill=color, width=5, joint="curve")
        for x, y in coords:
            draw.ellipse((x - 5, y - 5, x + 5, y + 5), fill=color)
        legend_x = plot_left + 20 + index * 440
        legend_y = plot_top + 8
        draw.line((legend_x, legend_y + 11, legend_x + 40, legend_y + 11), fill=color, width=5)
        draw.text((legend_x + 50, legend_y), label, fill=INK, font=font(18))


def main() -> None:
    if not OUT.is_dir():
        raise FileNotFoundError(OUT)
    book = read(BOOK)
    all_urea = read(UREA)
    urea = [r for r in all_urea if r["scenario"] == "C1_T85"]
    urea_c2 = [r for r in all_urea if r["scenario"] == "C2_T85"]
    assert len(book) == 41 and len(urea) == len(urea_c2) == 41
    written = []
    for h in (65, 700, 2600):
        image = Image.new("RGB", (W, H), BG)
        draw = ImageDraw.Draw(image)
        draw.text((65, 28), f"G3 旧盐曲线与尿素条件靶值  |  标称厚度 {h} nm",
                  fill=INK, font=font(36))
        draw.text((65, 79), "Book2既有观测 + 22号条件计算；不是新实验或机制验证",
                  fill="#5f6671", font=font(20))
        salt_lines = []
        urea_lines = []
        urea_bands = []
        for state, color in (("B", BLUE), ("R", ORANGE)):
            sid = f"{state}{h}"
            points = sorted((float(r["time_min"]), float(r["kappa_uScm"])) for r in book
                            if r["sample_id"] == sid)
            initial = points[0][1]
            salt_lines.append((f"{state} 原电导减初值", [(t, k - initial) for t, k in points], color))
            target = sorted((float(r["source_time_min"]), float(r["urea_receiver_W_mol_m3"]))
                            for r in urea if r["sample_id"] == sid)
            alt = sorted((float(r["source_time_min"]), float(r["urea_receiver_W_mol_m3"]))
                         for r in urea_c2 if r["sample_id"] == sid)
            assert [t for t, _ in target] == [t for t, _ in alt]
            low = [(t, min(v1, v2)) for (t, v1), (_, v2) in zip(target, alt)]
            high = [(t, max(v1, v2)) for (t, v1), (_, v2) in zip(target, alt)]
            urea_bands.append((low, high, "#d8e8f6" if state == "B" else "#fbe2d1"))
            label = "B / C锚点（C1线）" if state == "B" else "R若仅W（C1线）"
            urea_lines.append((label, target, color))
        plot_panel(draw, (65, 135, 1330, 510), "既有盐电导记录", "接收侧电导增量 / µS·cm⁻¹",
                   2200, [0, 500, 1000, 1500, 2000], salt_lines)
        plot_panel(draw, (65, 535, 1330, 910), "尿素条件靶值", "接收侧尿素 / mM（条件计算）",
                   0.5, [0, 0.1, 0.2, 0.3, 0.4, 0.5], urea_lines, urea_bands)
        draw.text((65, 925), "线为C1_T85，浅色带为C1–C2；原第5点85 min，T80见CSV。幂式要求共享几何与外阻按扩散比例缩放。",
                  fill="#5f6671", font=font(18))
        if h == 2600:
            draw.text((65, 951), "R2600原表缺120 min，图中未补点。外阻未知时不使用此靶值作唯一预测。",
                      fill="#5f6671", font=font(18))
        path = OUT / f"22_G3_旧盐与尿素条件靶值_{h}nm_20261028.png"
        image.save(path)
        written.append(path)
    print(f"G3 images={len(written)}: " + ", ".join(p.name for p in written))


if __name__ == "__main__":
    main()
