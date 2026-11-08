"""Repeat prespecified time and 2600-nm holdouts for two worst WP5 scenarios."""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "effective_network_model"))
from v6_fit_wp4 import csv_write, load_films, unpack  # noqa: E402
from v6_wp5_holdouts import MANIFEST_SHA256, trimmed  # noqa: E402
from v6_wp5_sensitivity import (SCENARIOS, directed_films, fit, perturbed_films,
                                scenario_prediction)  # noqa: E402

OUT = ROOT / "WP5"


def reconstruct(name, films):
    fixed = dict(SCENARIOS)
    if name in fixed:
        return fixed[name], films
    if name.startswith("reading_"):
        suffix = name.removeprefix("reading_")
        epsilon = float(suffix.split("_")[0])
        if "directed200" in suffix:
            return {}, directed_films(films, epsilon)
        sample_id = int(suffix.rsplit("seed", 1)[1])
        return {}, perturbed_films(films, epsilon, sample_id)
    raise ValueError(name)


def transformed(p):
    return np.array([p["ws_ref"], math.log10(p["ell_nm"]), math.log10(p["dm_ratio"]),
                     math.log10(p["gamma_per_m"]), math.log10(p["q_conductivity"])])


def main():
    films = load_films()
    base = json.loads((ROOT / "WP4" / "M1_拟合摘要.json").read_text(encoding="utf-8"))
    base_p = base["parameters"]
    with (OUT / "07_固定情景与读数扰动.csv").open(encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.DictReader(fh))
    ranked = []
    widths = np.array([20, 3, 12, 8, 8])
    for row in rows:
        if not row["parameters"]:
            continue
        p = json.loads(row["parameters"])
        drift = float(np.max(np.abs((transformed(p) - transformed(base_p)) / widths)))
        loss_change = abs(float(row["loss_bdf"]) - base["best_loss"]) / base["best_loss"]
        ranked.append((max(drift, loss_change), row["scenario_id"], drift, loss_change))
    ranked.sort(reverse=True)
    selected = ranked[:2]
    result = []
    details = {}
    for rank, (score, name, drift, loss_change) in enumerate(selected, 1):
        scenario, data = reconstruct(name, films)
        for split, train, target in (
            ("train_0_60", [trimmed(f) for f in data], data),
            ("hold_2600nm", [f for f in data if f.h_nm != 2600],
             [next(f for f in data if f.h_nm == 2600)]),
        ):
            best, runs = fit(train, scenario, max_nfev=120)
            details[f"{name}/{split}"] = {"best": best, "runs": runs}
            if best is None:
                result.append({"worst_rank": rank, "scenario_id": name, "split_id": split,
                               "data_manifest_hash": MANIFEST_SHA256,
                               "selection_score": score, "parameter_drift": drift,
                               "loss_change_fraction": loss_change, "model_id": "M1",
                               "E_macro": "", "E_2600": "", "endpoint_error_factor_2600": "",
                               "optimizer_status": "all starts failed",
                               "success_starts": 0, "failed_starts": 12})
                continue
            p = unpack(np.array(best["vector"]), "M1")
            per_film = []
            endpoint_factor = ""
            for film in target:
                y, _ = scenario_prediction(film, p, scenario, fast=False)
                mask = film.times_min > 60 if split == "train_0_60" else film.times_min > 0
                err = float(np.sqrt(np.mean((np.arcsinh(y[mask]) -
                                             np.arcsinh(film.increment[mask])) ** 2)))
                per_film.append((film.h_nm, err))
                if film.h_nm == 2600:
                    endpoint_factor = max(y[-1] / film.increment[-1],
                                          film.increment[-1] / max(y[-1], 1e-300))
            result.append({"worst_rank": rank, "scenario_id": name, "split_id": split,
                           "data_manifest_hash": MANIFEST_SHA256,
                           "selection_score": score, "parameter_drift": drift,
                           "loss_change_fraction": loss_change, "model_id": "M1",
                           "E_macro": float(np.mean([x[1] for x in per_film])),
                           "E_2600": next((x[1] for x in per_film if x[0] == 2600), ""),
                           "endpoint_error_factor_2600": endpoint_factor,
                           "optimizer_status": best["message"],
                           "success_starts": sum(r["success"] for r in runs),
                           "failed_starts": sum(not r["success"] for r in runs)})
            print(name, split, result[-1]["E_macro"], flush=True)
    csv_write(OUT / "10_最不利情景重复留出.csv", result)
    (OUT / "11_最不利情景搜索记录.json").write_text(json.dumps(
        {"selection_rule": "largest max of parameter-coordinate drift normalized by frozen box width and absolute full-data loss change fraction; uses all 44 prescribed scenarios, no independent validation claim",
         "ranked": ranked, "selected": selected, "details": details},
        ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
