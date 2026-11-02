"""Audit the immutable 22-series technical freeze against its saved sources."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import pdfplumber


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/salt_transport"
FREEZE = DATA / "freeze_predictions.json"
PDF = ROOT / "paper/20261028_22号实验前预测封存摘要.pdf"
OUT = DATA / "22_freeze_review_20261028.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def rows(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    frozen = json.loads(FREEZE.read_text(encoding="utf-8"))
    issues: list[str] = []
    checked_files = 0
    for category in ("source_input_manifest", "code_manifest", "derived_artifact_manifest"):
        for entry in frozen[category]:
            path = ROOT / entry["path"]
            checked_files += 1
            if not path.is_file():
                issues.append(f"missing:{entry['path']}")
            elif sha256(path) != entry["sha256"]:
                issues.append(f"hash_mismatch:{entry['path']}")
    figures = {row["path"]: row["sha256"] for row in frozen["figure_sources"]}
    artifacts = {row["path"]: row["sha256"] for row in frozen["derived_artifact_manifest"]}
    if len(figures) != 5 or any(artifacts.get(path) != digest for path, digest in figures.items()):
        issues.append("figure_manifest_not_subset")

    matrix_path = ROOT / "原始资料/赵珈绮/22_机制预测与实验验证方案_v1.0.0/matrices/paired_predictions.json"
    matrix = {row["pair_id"]: row for row in json.loads(matrix_path.read_text(encoding="utf-8-sig"))}
    effects = {row["pair_id"]: row for row in rows("22_paired_effects_20261028.csv")}
    curves = defaultdict(dict)
    for row in rows("22_paired_illustrations_20261028.csv"):
        if float(row["time_min"]) == 0:
            continue
        key = (row["state"], row["probe"], int(float(row["time_min"])))
        curves[row["pair_id"]][key] = row
    groups = frozen["prediction_groups"]
    if len(groups) != 54 or len({g["pair_id"] for g in groups}) != 54:
        issues.append("prediction_group_count_or_identity")
    for group in groups:
        pid = group["pair_id"]
        if pid not in matrix or pid not in effects:
            issues.append(f"pair_missing_source:{pid}")
            continue
        source = matrix[pid]
        for field in ("h_nm", "c_feed_mM", "mechanism", "external", "b_B", "b_E",
                      "sigma_B_C_m2", "sigma_E_C_m2"):
            stored = group["salt_feed_mM"] if field == "c_feed_mM" else group[field]
            if stored != source[field]:
                issues.append(f"pair_field_mismatch:{pid}:{field}")
        for saved, source_field in (("initial_membrane_E_s", "E_s_initial"),
                                    ("initial_membrane_E_u", "E_u_initial"),
                                    ("initial_membrane_Pi", "Pi_initial")):
            if group[saved] != float(effects[pid][source_field]):
                issues.append(f"effect_mismatch:{pid}:{saved}")
        if group["numeric_status"] != effects[pid]["numeric_status"]:
            issues.append(f"numeric_status_mismatch:{pid}")
        seen = set()
        for state in group["pool_states"]:
            key = (state["state"], state["probe"], state["time_min"])
            seen.add(key)
            source_state = curves[pid].get(key)
            if source_state is None or state["donor_mM"] != float(source_state["donor_mol_m3"]) or \
                    state["receiver_mM"] != float(source_state["receiver_mol_m3"]) or \
                    state["numeric_status"] != source_state["numeric_status"]:
                issues.append(f"pool_state_mismatch:{pid}:{key}")
        if len(seen) != 12 or len(curves[pid]) != 12:
            issues.append(f"pool_state_count:{pid}")

    experiment_rows = rows("22_G4_实验行判别地图_20261028.csv")
    if frozen["experiment_distinguishability_rows"] != experiment_rows or len(experiment_rows) != 54:
        issues.append("experiment_map_mismatch")
    if {row.get("run_id") for row in experiment_rows} != {f"H{i:03d}" for i in range(1, 55)}:
        issues.append("experiment_run_ids_incomplete")
    sensitivity = rows("22_algebraic_sensitivity_20261028.csv")
    sensitivity_keys = {(r["pair_id"], r["state"], r["probe"], r["time_min"]) for r in sensitivity}
    if len(sensitivity) != 648 or len(sensitivity_keys) != 648:
        issues.append("sensitivity_rows_or_keys")
    if {r["pair_id"] for r in sensitivity} != set(matrix):
        issues.append("sensitivity_pair_ids")
    if any(not (float(r["min_J_mol_m2_s"]) <= float(r["J_base_mol_m2_s"])
                <= float(r["max_J_mol_m2_s"])) for r in sensitivity):
        issues.append("sensitivity_envelope_excludes_base")
    ledger = rows("22_pb_call_ledger_20261028.csv")
    counts = Counter(row["status"] for row in ledger)
    if (len(ledger), counts["PASS"], counts["FAILED"]) != (437, 432, 5):
        issues.append("pb_ledger_changed")
    if frozen["pb_budget"] != {"actual_calls": 437, "limit": 444, "pass": 432, "failed": 5}:
        issues.append("pb_budget_mismatch")
    target = rows("22_pb_target_recheck_20261028.csv")
    if sum(row["status"] == "NUMERIC_UNRESOLVED_OVER_2_PERCENT" for row in target) != 10:
        issues.append("T12_unresolved_changed")
    if frozen["low_salt_candidates_mM"] != [1, 10]:
        issues.append("low_salt_candidates_changed")
    with pdfplumber.open(PDF) as document:
        text = "\n".join(page.extract_text() or "" for page in document.pages)
        pages = len(document.pages)
    if pages != 1 or sha256(FREEZE) not in text or "T12未决" not in text:
        issues.append("pdf_json_or_status_mismatch")

    result = {"result": "PASS" if not issues else "FAIL", "checked_manifest_files": checked_files,
              "prediction_groups": len(groups), "experiment_rows": len(experiment_rows),
              "sensitivity_rows": len(sensitivity),
              "pdf_pages": pages, "issues": issues,
              "scope": "Manifest, source-row, pair identity, pool state, experiment map, PB budget and PDF linkage; no model rerun."}
    with OUT.open("w", encoding="utf-8") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2)
    print(json.dumps(result, ensure_ascii=False))
    if issues:
        raise ArithmeticError(f"Freeze review found {len(issues)} issue(s)")


if __name__ == "__main__":
    main()
