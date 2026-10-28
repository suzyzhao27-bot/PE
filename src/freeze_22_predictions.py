"""Create immutable Codex technical freeze from fixed pre-experiment artifacts."""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/salt_transport"
PACK = ROOT / "原始资料/赵珈绮/22_机制预测与实验验证方案_v1.0.0"
OUT = DATA / "freeze_predictions.json"


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def manifest(paths: list[Path]) -> list[dict[str, str]]:
    result = []
    for path in sorted(set(paths)):
        if not path.is_file():
            raise FileNotFoundError(path)
        result.append({"path": path.relative_to(ROOT).as_posix(), "sha256": digest(path)})
    return result


def read(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def main() -> None:
    if OUT.exists():
        raise FileExistsError(f"Frozen artifact exists; create a new version instead: {OUT}")
    pair_path = PACK / "matrices/paired_predictions.json"
    pairs = json.loads(pair_path.read_text(encoding="utf-8-sig"))
    effects = {r["pair_id"]: r for r in read("22_paired_effects_20261028.csv")}
    curves = read("22_paired_illustrations_20261028.csv")
    sensitivity = read("22_algebraic_sensitivity_20261028.csv")
    experiments = read("22_G4_实验行判别地图_20261028.csv")
    ledger = read("22_pb_call_ledger_20261028.csv")
    target = read("22_pb_target_recheck_20261028.csv")
    if (len(pairs), len(effects), len(curves), len(sensitivity), len(experiments), len(ledger)) != (
            54, 54, 864, 648, 54, 437):
        raise ArithmeticError("Frozen matrix or result count changed")
    if sum(r["status"] == "NUMERIC_UNRESOLVED_OVER_2_PERCENT" for r in target) != 10:
        raise ArithmeticError("T12 unresolved count changed")
    by_pair: dict[str, list[dict[str, str]]] = defaultdict(list)
    by_sensitivity: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in curves:
        by_pair[row["pair_id"]].append(row)
    for row in sensitivity:
        by_sensitivity[row["pair_id"]].append(row)
    prediction_groups = []
    for pair in pairs:
        pid = pair["pair_id"]
        if len(by_pair[pid]) != 16 or len(by_sensitivity[pid]) != 12:
            raise ArithmeticError(f"Incomplete pair {pid}")
        times = []
        for row in by_pair[pid]:
            if float(row["time_min"]) == 0:
                continue
            times.append({"state": row["state"], "probe": row["probe"],
                          "time_min": int(float(row["time_min"])),
                          "donor_mM": float(row["donor_mol_m3"]),
                          "receiver_mM": float(row["receiver_mol_m3"]),
                          "numeric_status": row["numeric_status"]})
        prediction_groups.append({
            "pair_id": pid, "h_nm": pair["h_nm"], "salt_feed_mM": pair["c_feed_mM"],
            "mechanism": pair["mechanism"], "external": pair["external"],
            "b_B": pair["b_B"], "b_E": pair["b_E"],
            "sigma_B_C_m2": pair["sigma_B_C_m2"], "sigma_E_C_m2": pair["sigma_E_C_m2"],
            "initial_membrane_E_s": float(effects[pid]["E_s_initial"]),
            "initial_membrane_E_u": float(effects[pid]["E_u_initial"]),
            "initial_membrane_Pi": float(effects[pid]["Pi_initial"]),
            "numeric_status": effects[pid]["numeric_status"],
            "pool_states": times,
            "instantaneous_flux_envelope_source": "22_algebraic_sensitivity_20261028.csv",
        })
    code_paths = list((ROOT / "src").glob("mechanism_*.py"))
    code_paths += list((ROOT / "src").glob("plot_22_*.py"))
    code_paths += [ROOT / "src/freeze_22_predictions.py", ROOT / "src/render_22_freeze_pdf.py"]
    input_paths = [
        PACK / "02_方程参数与实现规格.md", PACK / "03_计算矩阵与验收.md",
        PACK / "04_考试后实验与判别.md", PACK / "05_复用清单与数据接口.md",
        PACK / "matrices/pb_main.json", PACK / "matrices/pb_refinement.json",
        pair_path, PACK / "matrices/experiment_runs.json",
        ROOT / "原始资料/赵珈绮/20_学生Stage2_B1成果_v1.0.0/WP1/book2_long.csv",
    ]
    result_paths = [p for p in DATA.glob("22_*") if p.is_file() and p != OUT]
    result_paths += list((ROOT / "阶段总结").glob("20261028_22号*.md"))
    artifact_manifest = manifest(result_paths)
    now = datetime.now(ZoneInfo("Asia/Shanghai")).isoformat(timespec="seconds")
    frozen = {
        "schema": "22_pre_experiment_conditional_predictions_v1",
        "status": "FROZEN_CODEX_TECHNICAL_STUDENT_UNREVIEWED",
        "created_at_Asia_Shanghai": now,
        "actual_operator": "Codex on local dev-zhao workspace; student did not run or sign",
        "code_version": {"branch": git("branch", "--show-current"), "base_commit": git("rev-parse", "HEAD"),
                         "working_tree": "UNCOMMITTED_SNAPSHOT_SHA256_PER_FILE"},
        "pb_budget": {"actual_calls": len(ledger), "limit": 444,
                      "pass": sum(r["status"] == "PASS" for r in ledger),
                      "failed": sum(r["status"] == "FAILED" for r in ledger)},
        "T12": {"target_count": 36, "unresolved_count": 10,
                "status": "NUMERIC_UNRESOLVED_OVER_2_PERCENT",
                "source": "data/processed/salt_transport/22_pb_target_recheck_20261028.csv"},
        "low_salt_candidates_mM": [1, 10],
        "low_salt_selection": "PENDING_POST_EXAM_METHOD_CALIBRATION_5X_LOQ_RULE",
        "interpretation": "Uniform illustrative model only; not observed membrane properties, material fit, or experimental validation.",
        "source_input_manifest": manifest(input_paths),
        "code_manifest": manifest(code_paths),
        "derived_artifact_manifest": artifact_manifest,
        "figure_sources": [r for r in artifact_manifest if r["path"].lower().endswith(".png")],
        "prediction_groups": prediction_groups,
        "experiment_distinguishability_rows": experiments,
    }
    if len(frozen["figure_sources"]) != 5:
        raise ArithmeticError("Expected G1/G2/G3 figures")
    with OUT.open("x", encoding="utf-8") as stream:
        json.dump(frozen, stream, ensure_ascii=False, indent=2)
    print(f"frozen_pairs={len(prediction_groups)} experiment_rows={len(experiments)} "
          f"code={len(frozen['code_manifest'])} inputs={len(frozen['source_input_manifest'])} "
          f"artifacts={len(frozen['derived_artifact_manifest'])}")


if __name__ == "__main__":
    main()
