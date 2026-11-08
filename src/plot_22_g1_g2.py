"""Render G1 fixed-thickness resistance and G2 paired-probe illustration."""

from __future__ import annotations

import csv
import math
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/salt_transport"
FONT = Path("C:/Windows/Fonts/msyh.ttc")
INK, GRID, BLUE, ORANGE = "#243144", "#d7dde6", "#2b6cb0", "#df6c28"


def read(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def font(size: int) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(str(FONT), size)


def g1_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], title: str,
             rows: list[dict[str, str]]) -> None:
    left, top, right, bottom = box
    xl, xr, yt, yb = left + 105, right - 25, top + 62, bottom - 68
    draw.text((left, top), title, fill=INK, font=font(27))
    for exponent in range(-8, -2):
        y = yb - (exponent + 8) / 5 * (yb - yt)
        draw.line((xl, y, xr, y), fill=GRID, width=2)
        draw.text((left + 8, y - 12), f"10^{exponent}", fill=INK, font=font(17))
    for h in (65, 700, 2600):
        x = xl + (math.log10(h) - math.log10(65)) / (math.log10(2600) - math.log10(65)) * (xr - xl)
        draw.line((x, yt, x, yb), fill=GRID, width=1)
        draw.text((x - 22, yb + 14), str(h), fill=INK, font=font(17))
    draw.line((xl, yt, xl, yb), fill=INK, width=3)
    draw.line((xl, yb, xr, yb), fill=INK, width=3)
    draw.text((xr - 100, yb + 37), "h / nm", fill=INK, font=font(18))
    legend = []
    for b, probe, color, label in (
        (.0001, "salt_zero_charge", BLUE, "盐 b=0.0001"),
        (.0001, "urea", "#73a2d1", "尿素 b=0.0001"),
        (.01, "salt_zero_charge", ORANGE, "盐 b=0.01"),
        (.01, "urea", "#e8a46f", "尿素 b=0.01"),
    ):
        points = sorted((float(r["h_nm"]), float(r["permeability_total_m_s"])) for r in rows
                        if float(r["b"]) == b and r["probe"] == probe
                        and float(r["salt_feed_mM"]) == 1 and float(r["time_min"]) == 0)
        assert len(points) == 3
        coords = [(xl + (math.log10(h) - math.log10(65)) /
                   (math.log10(2600) - math.log10(65)) * (xr - xl),
                   yb - (math.log10(p) + 8) / 5 * (yb - yt)) for h, p in points]
        draw.line(coords, fill=color, width=5)
        for x, y in coords:
            draw.ellipse((x - 6, y - 6, x + 6, y + 6), fill=color)
        legend.append((label, color))
    for i, (label, color) in enumerate(legend):
        x = xl + (i % 2) * 285
        y = yt + 14 + (i // 2) * 30
        draw.line((x, y + 10, x + 34, y + 10), fill=color, width=5)
        draw.text((x + 43, y), label, fill=INK, font=font(17))


def render_g1() -> Path:
    rows = read("22_neutral_baseline_20261028.csv")
    image = Image.new("RGB", (1450, 780), "white")
    draw = ImageDraw.Draw(image)
    draw.text((55, 28), "G1 固定结构的厚度与外阻效应", fill=INK, font=font(36))
    draw.text((55, 79), "零电荷解析基线；纵轴为总体系P / (m/s)，对数；b为指定情景，不是实测b(h)",
              fill="#5f6671", font=font(20))
    zero = [r for r in rows if r["external"] == "zero"]
    reference = [r for r in rows if r["external"] == "reference_2600_neutral"]
    g1_panel(draw, (55, 140, 700, 690), "外阻 = 0", zero)
    g1_panel(draw, (750, 140, 1390, 690), "固定参考外阻", reference)
    draw.text((55, 726), "固定b时增厚使P下降；旧厚端较快不能仅由h解释。孔口及外阻为22号条件近似。",
              fill="#5f6671", font=font(20))
    path = DATA / "22_G1_厚度外阻条件图_20261028.png"
    image.save(path)
    return path


def g2_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], title: str,
             rows: list[dict[str, str]], ymax: float) -> None:
    left, top, right, bottom = box
    xl, xr, yt, yb = left + 52, right - 20, top + 52, bottom - 51
    draw.text((left, top), title, fill=INK, font=font(24))
    for fraction in (0, .25, .5, .75, 1):
        y = yb - fraction * (yb - yt)
        draw.line((xl, y, xr, y), fill=GRID, width=1)
        draw.text((left + 5, y - 11), f"{ymax*fraction:g}", fill=INK, font=font(15))
    for minute in (0, 60, 120, 180, 240):
        x = xl + minute / 240 * (xr - xl)
        draw.line((x, yt, x, yb), fill=GRID, width=1)
        draw.text((x - 13, yb + 10), str(minute), fill=INK, font=font(15))
    draw.line((xl, yt, xl, yb), fill=INK, width=2)
    draw.line((xl, yb, xr, yb), fill=INK, width=2)
    for state, color in (("B", BLUE), ("E", ORANGE)):
        points = sorted((float(r["time_min"]), float(r["receiver_mol_m3"])) for r in rows
                        if r["state"] == state)
        assert len(points) == 4
        coords = [(xl + t / 240 * (xr - xl), yb - v / ymax * (yb - yt)) for t, v in points]
        draw.line(coords, fill=color, width=4)
        for x, y in coords:
            draw.ellipse((x - 4, y - 4, x + 4, y + 4), fill=color)
        legend_x = xl + (0 if state == "B" else 130)
        draw.text((legend_x, yt + 8), state, fill=color, font=font(18))


def render_g2() -> Path:
    curves = read("22_paired_illustrations_20261028.csv")
    effects = read("22_paired_effects_20261028.csv")
    image = Image.new("RGB", (1500, 950), "white")
    draw = ImageDraw.Draw(image)
    draw.text((52, 25), "G2 两探针成对示例：700 nm，1 mM，零外阻", fill=INK, font=font(35))
    draw.text((52, 77), "接收侧浓度；盐初始1 mM，尿素初始10 mM；W/C/WC为22号固定情景，非材料拟合",
              fill="#5f6671", font=font(19))
    for col, mechanism in enumerate(("W", "C", "WC")):
        selected = [r for r in curves if int(float(r["h_nm"])) == 700 and
                    int(float(r["c_feed_mM"])) == 1 and r["external"] == "zero" and
                    r["mechanism"] == mechanism]
        effect = next(r for r in effects if int(float(r["h_nm"])) == 700 and
                      int(float(r["c_feed_mM"])) == 1 and r["external"] == "zero" and
                      r["mechanism"] == mechanism)
        left = 52 + col * 485
        right = left + 455
        g2_panel(draw, (left, 145, right, 450), f"{mechanism}：盐 / mM",
                 [r for r in selected if r["probe"] == "salt"], 0.5)
        g2_panel(draw, (left, 495, right, 800), f"{mechanism}：尿素 / mM",
                 [r for r in selected if r["probe"] == "urea"], 5)
        draw.text((left + 7, 823),
                  f"Es={float(effect['E_s_initial']):.3g}  Eu={float(effect['E_u_initial']):.3g}  Pi={float(effect['Pi_initial']):.3g}",
                  fill=INK, font=font(18))
    draw.text((52, 881), "Es/Eu/Pi为初始膜内有效P之比，非终点浓度比。蓝=B，橙=E；60/120/240 min取样。",
              fill="#5f6671", font=font(18))
    draw.text((52, 912), "此图仅展示一个固定情景；其余53对及10 mM数值未决标记见CSV。",
              fill="#5f6671", font=font(18))
    path = DATA / "22_G2_两探针情景_700nm_1mM_20261028.png"
    image.save(path)
    return path


def main() -> None:
    if not DATA.is_dir():
        raise FileNotFoundError(DATA)
    a, b = render_g1(), render_g2()
    print(f"G1={a.name}; G2={b.name}")


if __name__ == "__main__":
    main()
