"""G4 conditional detection and flow limits plus 54-run experiment crosswalk."""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "原始资料/赵珈绮/22_机制预测与实验验证方案_v1.0.0"
CURVES = ROOT / "data/processed/salt_transport/22_paired_illustrations_20261028.csv"
MODEL_OUT = ROOT / "data/processed/salt_transport/22_G4_模型检测需求_20261028.csv"
RUN_OUT = ROOT / "data/processed/salt_transport/22_G4_实验行判别地图_20261028.csv"
MODEL_FIELDS = ["pair_id", "h_nm", "salt_feed_mM", "mechanism", "external", "state",
                "salt_receiver_60_uM", "urea_receiver_60_uM", "salt_LOQ_at_most_uM",
                "urea_LOQ_at_most_uM", "salt_equilibrium_fraction_60",
                "urea_equilibrium_fraction_60", "flow10_0_60_uL", "flow30_0_60_uL",
                "flow10_60_120_uL", "flow30_60_120_uL",
                "flow10_120_240_uL", "flow30_120_240_uL",
                "numeric_status", "actual_LOQ_status"]
RUN_FIELDS = ["run_id", "group", "thickness_level", "parent_film_id", "state",
              "salt_feed", "urea_feed_mM", "role", "purpose", "illustration_link",
              "conditional_salt_60_uM_range_1mM", "conditional_salt_60_uM_range_10mM",
              "conditional_salt_60_uM_range_600mM", "conditional_urea_60_uM_range",
              "actual_LOQ_status", "actual_flow_status", "identifiability_status"]


def csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def flow_threshold(series: dict[float, dict[str, str]], t0: float, t1: float) -> float:
    start, end = series[t0], series[t1]
    volume_m3 = {0.0: 1e-5, 60.0: 9.8e-6, 120.0: 9.6e-6}[t0]
    amount = volume_m3 * (float(end["receiver_mol_m3"]) - float(start["receiver_mol_m3"]))
    reference_c = max(float(start["donor_mol_m3"]), float(start["receiver_mol_m3"]),
                      float(end["donor_mol_m3"]), float(end["receiver_mol_m3"]))
    if amount < -1e-12 or reference_c <= 0:
        raise ArithmeticError("Invalid model transfer or concentration")
    return max(0.0, amount) / reference_c * 1e9  # microliters, before 10% factor


def model_map(curves: list[dict[str, str]], pairs: list[dict]) -> list[dict]:
    grouped = defaultdict(list)
    for row in curves:
        grouped[row["pair_id"], row["state"], row["probe"]].append(row)
    output = []
    for pair in pairs:
        for state in ("B", "E"):
            salt = {float(r["time_min"]): r for r in grouped[pair["pair_id"], state, "salt"]}
            urea = {float(r["time_min"]): r for r in grouped[pair["pair_id"], state, "urea"]}
            assert set(salt) == set(urea) == {0.0, 60.0, 120.0, 240.0}
            s60 = float(salt[60.0]["receiver_mol_m3"]) * 1000  # mM = mol/m^3; to micromolar
            u60 = float(urea[60.0]["receiver_mol_m3"]) * 1000
            limits = []
            for probe, series in (("salt", salt), ("urea", urea)):
                for t0, t1 in ((0.0, 60.0), (60.0, 120.0), (120.0, 240.0)):
                    limits.append((0.1 * flow_threshold(series, t0, t1), probe, f"{t0:g}-{t1:g}"))
            by_window = {w: min(v for v, _, ww in limits if ww == w)
                         for w in ("0-60", "60-120", "120-240")}
            output.append({
                "pair_id": pair["pair_id"], "h_nm": pair["h_nm"],
                "salt_feed_mM": pair["c_feed_mM"], "mechanism": pair["mechanism"],
                "external": pair["external"], "state": state,
                "salt_receiver_60_uM": s60, "urea_receiver_60_uM": u60,
                "salt_LOQ_at_most_uM": s60 / 5, "urea_LOQ_at_most_uM": u60 / 5,
                "salt_equilibrium_fraction_60": 2 * s60 / (1000 * float(pair["c_feed_mM"])),
                "urea_equilibrium_fraction_60": 2 * u60 / 10000,
                "flow10_0_60_uL": by_window["0-60"],
                "flow30_0_60_uL": 3 * by_window["0-60"],
                "flow10_60_120_uL": by_window["60-120"],
                "flow30_60_120_uL": 3 * by_window["60-120"],
                "flow10_120_240_uL": by_window["120-240"],
                "flow30_120_240_uL": 3 * by_window["120-240"],
                "numeric_status": ("NUMERIC_UNRESOLVED_T12" if pair["c_feed_mM"] == 10
                                   and pair[f"sigma_{state}_C_m2"] != 0
                                   else "MODEL_ONLY_TARGET_CHECK"),
                "actual_LOQ_status": "NOT_MEASURED",
            })
    return output


PURPOSE = {
    "E1": "thin/thick paired B/E and low/high salt; evaluate E_s,E_u,Pi",
    "E2": "middle high-salt held-out thickness response",
    "E3": "water-wait control; not paired-model W mechanism",
    "E4": "three-layer serial-path control; not a separate paired illustration",
    "C1": "support-only external resistance reference",
    "C2": "impermeable seal leakage/background",
    "C3": "zero-solute PE leachate and method blank",
    "C4": "urea-free salt perturbation check",
}


def experiment_map(runs: list[dict], model_rows: list[dict]) -> list[dict]:
    nominal = {"thin": 65, "middle": 700, "thick": 2600}
    output = []
    for run in runs:
        group, level, state = run["group"], run["thickness_level"], run["PE_state"]
        model_link = group in ("E1", "E2") and level in nominal and state in ("B", "E")
        candidates = [r for r in model_rows if model_link and int(r["h_nm"]) == nominal[level]
                      and r["state"] == state]
        ranges = {}
        for feed in (1, 10, 600):
            subset = [r for r in candidates if int(r["salt_feed_mM"]) == feed]
            ranges[feed] = (f"{min(r['salt_receiver_60_uM'] for r in subset):.6g}.."
                            f"{max(r['salt_receiver_60_uM'] for r in subset):.6g}") if subset else "NA"
        urea_range = (f"{min(r['urea_receiver_60_uM'] for r in candidates):.6g}.."
                      f"{max(r['urea_receiver_60_uM'] for r in candidates):.6g}") if candidates else "NA"
        output.append({
            "run_id": run["run_id"], "group": group, "thickness_level": level,
            "parent_film_id": run["parent_film_id"], "state": state,
            "salt_feed": run["salt_feed_mM"], "urea_feed_mM": run["urea_feed_mM"],
            "role": run["role"], "purpose": PURPOSE[group],
            "illustration_link": ("CONDITIONAL_NOMINAL_H_NOT_MATERIAL_FORECAST" if model_link
                                  else "NO_DIRECT_54_PAIR_ILLUSTRATION"),
            "conditional_salt_60_uM_range_1mM": ranges[1] if run["salt_level"] == "low" else "NA",
            "conditional_salt_60_uM_range_10mM": ranges[10] if run["salt_level"] == "low" else "NA",
            "conditional_salt_60_uM_range_600mM": ranges[600] if run["salt_level"] == "high" else "NA",
            "conditional_urea_60_uM_range": urea_range,
            "actual_LOQ_status": "NOT_MEASURED",
            "actual_flow_status": "NOT_MEASURED",
            "identifiability_status": "PRE_EXPERIMENT_UNDETERMINED",
        })
    return output


def main() -> None:
    curves = csv_rows(CURVES)
    with (PACKAGE / "matrices/paired_predictions.json").open(encoding="utf-8-sig") as stream:
        pairs = json.load(stream)
    with (PACKAGE / "matrices/experiment_runs.json").open(encoding="utf-8-sig") as stream:
        runs = json.load(stream)
    assert len(curves) == 864 and len(pairs) == len(runs) == 54
    model = model_map(curves, pairs)
    experiments = experiment_map(runs, model)
    assert len(model) == 108 and len(experiments) == 54
    write_csv(MODEL_OUT, model, MODEL_FIELDS)
    write_csv(RUN_OUT, experiments, RUN_FIELDS)
    print(f"model_rows={len(model)} experiment_rows={len(experiments)} "
          f"min_urea_LOQ_requirement_uM={min(r['urea_LOQ_at_most_uM'] for r in model):.6g} "
          f"min_flow10_0_60_uL={min(r['flow10_0_60_uL'] for r in model):.6g}")


if __name__ == "__main__":
    main()
