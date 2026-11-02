"""Render a human-readable hash-bearing summary of the immutable JSON freeze."""

from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/processed/salt_transport/freeze_predictions.json"
OUT = ROOT / "paper/20261028_22号实验前预测封存摘要.pdf"
FONT = "C:/Windows/Fonts/msyh.ttc"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fmt_range(values: list[float]) -> str:
    low, high = min(values), max(values)
    return f"{low:.3g}" if abs(high - low) < 1e-12 else f"{low:.3g}-{high:.3g}"


def main() -> None:
    if OUT.exists():
        raise FileExistsError(f"Frozen PDF exists; create a new version instead: {OUT}")
    data = json.loads(SOURCE.read_text(encoding="utf-8"))
    if len(data["prediction_groups"]) != 54 or len(data["experiment_distinguishability_rows"]) != 54:
        raise ArithmeticError("Incomplete freeze")
    pdfmetrics.registerFont(TTFont("YaHei", FONT, subfontIndex=0))
    title = ParagraphStyle("title", fontName="YaHei", fontSize=18, leading=25, textColor=colors.HexColor("#23364d"), spaceAfter=12)
    h2 = ParagraphStyle("h2", fontName="YaHei", fontSize=12, leading=17, textColor=colors.HexColor("#23364d"), spaceBefore=10, spaceAfter=6)
    body = ParagraphStyle("body", fontName="YaHei", fontSize=9, leading=15, spaceAfter=6)
    note = ParagraphStyle("note", fontName="YaHei", fontSize=8, leading=12, textColor=colors.HexColor("#4f5d6d"), spaceAfter=7)
    story = [Paragraph("22号实验前条件预测封存摘要", title),
             Paragraph(f"封存时间：{data['created_at_Asia_Shanghai']}（上海）", body),
             Paragraph("实际运行者：本机 Codex。学生本人未运行、未签署；导师未验收。", body),
             Paragraph("本 PDF 摘要对应下列 JSON 哈希。完整54对时序、54实验行、代码、输入及图源的 SHA256 清单均在 JSON；本 PDF 不替代 JSON。", note)]
    hash_rows = [["对象", "SHA256"]]
    hash_rows.append(["freeze_predictions.json", digest(SOURCE)])
    paths = {m["path"]: m["sha256"] for m in data["source_input_manifest"] + data["code_manifest"]}
    for label, suffix in (("成对矩阵", "matrices/paired_predictions.json"),
                          ("PB主矩阵", "matrices/pb_main.json"),
                          ("输运代码", "src/mechanism_transport.py"),
                          ("PB代码", "src/mechanism_radial_pb.py")):
        hit = [v for k, v in paths.items() if k.endswith(suffix)]
        if len(hit) != 1:
            raise ArithmeticError(f"Missing hash for {suffix}")
        hash_rows.append([label, hit[0]])
    hashes = Table(hash_rows, colWidths=[112, 394], repeatRows=1)
    hashes.setStyle(TableStyle([("FONTNAME", (0, 0), (-1, -1), "YaHei"),
                                ("FONTSIZE", (0, 0), (-1, -1), 7.2),
                                ("LEADING", (0, 0), (-1, -1), 11),
                                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#e8eef5")),
                                ("GRID", (0, 0), (-1, -1), .3, colors.HexColor("#c5cfd9")),
                                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                                ("LEFTPADDING", (0, 0), (-1, -1), 6),
                                ("RIGHTPADDING", (0, 0), (-1, -1), 4)]))
    story += [hashes, Spacer(1, 9), Paragraph("固定模型与验收边界", h2)]
    pb = data["pb_budget"]
    story += [Paragraph(f"PB 实际调用 {pb['actual_calls']}/{pb['limit']}：{pb['pass']} PASS、{pb['failed']} FAILED。"
                        "T12 有10/36个目标的31/16节点积分差仍超过2%，相关10 mM带电结果保留数值未决。", body),
              Paragraph("1与10 mM低盐候选同时封存。实际选择待考试后方法校准，按固定5倍LOQ规则执行；新实验结果不得回填本版。", body),
              Paragraph("以下数值为统一孔半径25 nm、迂曲度2、指定 b/σ/外阻的条件示例。E_s、E_u、Π是初始膜内有效透过系数比，非实测材料增益。", note),
              Paragraph("54对初始效应的分组范围", h2)]
    grouped: dict[tuple[str, float, str], list[dict[str, object]]] = defaultdict(list)
    for row in data["prediction_groups"]:
        grouped[(row["mechanism"], float(row["salt_feed_mM"]), row["external"])].append(row)
    rows = [["机制", "盐mM", "外阻", "E_s", "E_u", "Π", "状态"]]
    for (mechanism, feed, external), group in sorted(grouped.items()):
        status = "T12未决" if any(r["numeric_status"] == "NUMERIC_UNRESOLVED_T12" for r in group) else "条件值"
        rows.append([mechanism, f"{feed:g}", "0" if external == "zero" else "参考",
                     fmt_range([float(r["initial_membrane_E_s"]) for r in group]),
                     fmt_range([float(r["initial_membrane_E_u"]) for r in group]),
                     fmt_range([float(r["initial_membrane_Pi"]) for r in group]), status])
    table = Table(rows, colWidths=[48, 49, 53, 101, 72, 111, 72], repeatRows=1)
    table.setStyle(TableStyle([("FONTNAME", (0, 0), (-1, -1), "YaHei"),
                              ("FONTSIZE", (0, 0), (-1, -1), 8),
                              ("LEADING", (0, 0), (-1, -1), 11),
                              ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#e8eef5")),
                              ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f6f8fb")]),
                              ("GRID", (0, 0), (-1, -1), .3, colors.HexColor("#c5cfd9")),
                              ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                              ("LEFTPADDING", (0, 0), (-1, -1), 4),
                              ("RIGHTPADDING", (0, 0), (-1, -1), 4)]))
    story += [table, Spacer(1, 8),
              Paragraph("G1/G2/G3图源及G4逐行判别表均见 JSON 哈希清单。当前LOQ、净流和材料参数未测；54对示例不能替代按实际装置条件的预测。", note)]

    def footer(canvas, document):
        canvas.saveState()
        canvas.setFont("YaHei", 8)
        canvas.setFillColor(colors.HexColor("#667487"))
        canvas.drawString(42, 27, "22号实验前条件预测 - Codex技术封存摘要")
        canvas.drawRightString(A4[0] - 42, 27, f"第 {document.page} 页")
        canvas.restoreState()

    doc = SimpleDocTemplate(str(OUT), pagesize=A4, leftMargin=42, rightMargin=42,
                            topMargin=42, bottomMargin=46, title="22号实验前条件预测封存摘要")
    doc.build(story, onFirstPage=footer, onLaterPages=footer)
    print(f"pdf={OUT.name} json_sha256={digest(SOURCE)}")


if __name__ == "__main__":
    main()
