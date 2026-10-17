"""Render the existing Markdown working draft as an internal review PDF."""

import argparse
import html
import re
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    CondPageBreak,
    Image,
    KeepTogether,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)
from PIL import Image as PILImage


FONT = Path("C:/Windows/Fonts/simhei.ttf")
LINK = re.compile(r"!?\[([^\]]+)\]\((?:<)?[^)]+(?:>)?\)")
FIGURE = re.compile(r"!\[([^\]]+)\]\(([^)]+)\)")


def plain_markdown(value):
    value = LINK.sub(r"\1", value)
    value = re.sub(r"(?<!\w)\*([^*]+)\*", r"\1", value.replace("**", ""))
    value = value.translate(str.maketrans({
        "²": "^2", "ŷ": "y_hat", "ʰ": "^h", "ⁱ": "^i", "⁶": "^6",
        "⁻": "^-", "₀": "_0", "−": "-",
    }))
    return html.escape(value.replace("`", ""))


def parse_table(lines, styles, width):
    cells = [[part.strip() for part in line.strip().strip("|").split("|")] for line in lines]
    cells = [cells[0]] + cells[2:]
    columns = len(cells[0])
    if any(len(row) != columns for row in cells):
        raise ValueError("Inconsistent Markdown table columns")
    rows = [[Paragraph(plain_markdown(cell), styles["table"]) for cell in row] for row in cells]
    table = Table(rows, colWidths=[width / columns] * columns, repeatRows=1, hAlign="LEFT")
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#e7eef2")),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f7f9fa")]),
        ("GRID", (0, 0), (-1, -1), 0.3, colors.HexColor("#bdc8ce")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ("RIGHTPADDING", (0, 0), (-1, -1), 6),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
    ]))
    return table


def page_decoration(canvas, document):
    canvas.saveState()
    canvas.setFont("SimHei", 8)
    canvas.setFillColor(colors.HexColor("#65717a"))
    canvas.drawString(54, 32, "内部审阅稿 · 非学校模板 · 非学生正式提交")
    canvas.drawRightString(A4[0] - 54, 32, str(document.page))
    canvas.restoreState()


def render(source_path, output_path):
    if not FONT.is_file():
        raise FileNotFoundError(f"Chinese font unavailable: {FONT}")
    if not output_path.parent.is_dir():
        raise FileNotFoundError(f"Output directory does not exist: {output_path.parent}")
    pdfmetrics.registerFont(TTFont("SimHei", str(FONT)))
    styles = {
        "title": ParagraphStyle("title", fontName="SimHei", fontSize=17, leading=25,
                                textColor=colors.HexColor("#173c53"), spaceAfter=15),
        "h1": ParagraphStyle("h1", fontName="SimHei", fontSize=13, leading=21,
                             textColor=colors.HexColor("#17658d"), spaceBefore=16, spaceAfter=7),
        "h2": ParagraphStyle("h2", fontName="SimHei", fontSize=10.5, leading=18,
                             textColor=colors.HexColor("#173c53"), spaceBefore=11, spaceAfter=5),
        "body": ParagraphStyle("body", fontName="SimHei", fontSize=9.1, leading=16,
                               wordWrap="CJK", spaceAfter=8, alignment=TA_LEFT),
        "list": ParagraphStyle("list", fontName="SimHei", fontSize=9.1, leading=16,
                               wordWrap="CJK", leftIndent=12, firstLineIndent=-12,
                               spaceAfter=8),
        "table": ParagraphStyle("table", fontName="SimHei", fontSize=8.2, leading=13,
                                wordWrap="CJK"),
    }
    width = A4[0] - 108
    story = []
    lines = source_path.read_text(encoding="utf-8").splitlines()
    paragraph = []

    def flush():
        if paragraph:
            story.append(Paragraph(plain_markdown(" ".join(paragraph)), styles["body"]))
            paragraph.clear()

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            flush()
            i += 1
            continue
        figure = FIGURE.fullmatch(line)
        if figure:
            flush()
            image_path = (source_path.parent / figure.group(2)).resolve()
            if not image_path.is_file():
                raise FileNotFoundError(f"Figure missing: {image_path}")
            with PILImage.open(image_path) as bitmap:
                original_width, original_height = bitmap.size
            height = width * original_height / original_width
            story.append(CondPageBreak(height + 26))
            story.append(Image(str(image_path), width=width, height=height))
            story.append(Spacer(1, 11))
            i += 1
            continue
        if line.startswith("|"):
            flush()
            block = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                block.append(lines[i])
                i += 1
            story.append(KeepTogether([parse_table(block, styles, width), Spacer(1, 10)]))
            continue
        if line.startswith("# "):
            flush()
            story.append(Paragraph(plain_markdown(line[2:]), styles["title"]))
        elif line.startswith("## "):
            flush()
            story.append(CondPageBreak(255 if line.startswith("## 4.") else 90))
            story.append(Paragraph(plain_markdown(line[3:]), styles["h1"]))
        elif line.startswith("### "):
            flush()
            story.append(CondPageBreak(70))
            story.append(Paragraph(plain_markdown(line[4:]), styles["h2"]))
        elif re.match(r"(?:\d+\.|-)\s", line):
            flush()
            story.append(Paragraph(plain_markdown(line), styles["list"]))
        else:
            paragraph.append(line)
        i += 1
    flush()

    document = SimpleDocTemplate(str(output_path), pagesize=A4, leftMargin=54,
                                 rightMargin=54, topMargin=52, bottomMargin=57,
                                 title="PE v5 学生报告工作稿 - 内部审阅预览",
                                 author="PE research workspace")
    document.build(story, onFirstPage=page_decoration, onLaterPages=page_decoration)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    render(args.source.resolve(), args.output.resolve())


if __name__ == "__main__":
    main()
