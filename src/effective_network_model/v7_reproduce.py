"""Representative v7 replay from packaged inputs and archived fit parameters.

This checks forward calculations and figure sources. It does not repeat optimization.
"""

from __future__ import annotations

import csv
import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np

from v7_d1_objective import ROOT, curves, graphs, predict
from v7_network_transport import SPEC_PATH
from v7_straight_chain import straight_channel, closed_two_baths
from v7_paper_figures import figure1

TESTS = ROOT / "tests"
TASK = ROOT / "原始资料/赵珈绮/14_论文补强学生任务书_v7.0.0"


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    started = time.perf_counter()
    TESTS.mkdir(exist_ok=True)
    spec = json.loads(SPEC_PATH.read_text(encoding="utf-8-sig"))
    data, nets = curves(), graphs()
    inputs = {"frozen_spec": SPEC_PATH, "salt_tidy": TASK / "inputs/salt/salt_transport_tidy.csv",
              "curve_roles": TASK / "inputs/salt/curve_roles_v7.csv",
              "networks": ROOT / "WP2/v7_network_tables.jsonl.gz"}
    report = {"mode": "representative_replay_not_refit", "environment": {"python": sys.version.split()[0],
              "platform": platform.platform(), "numpy": np.__version__},
              "input_sha256": {name: sha(path) for name, path in inputs.items()}, "checks": {}}

    t0 = time.perf_counter()
    physics = spec["physics"]
    d = 2*physics["D_Na_m2_s"]*physics["D_Cl_m2_s"]/(physics["D_Na_m2_s"]+physics["D_Cl_m2_s"])
    p, lag, p_series = straight_channel(200e-9, 1.0, 128, False, d, physics["a_ref_m"])
    ref = next(r for r in read_csv(ROOT / "WP2/v7_straight_chain_benchmark.csv") if
               r["h_nm"] == "200" and r["K"] == "1.0" and r["access"] == "false" and r["cells"] == "128")
    bath = closed_two_baths(200e-9, 0.1, 32, d, physics["a_ref_m"], 1e3)
    report["checks"]["straight_chain"] = {"elapsed_s": time.perf_counter()-t0,
              "P_m_s": p, "P_relative_to_archived": abs(p-float(ref["P_m_s"]))/float(ref["P_m_s"]),
              "lag_relative_to_archived": abs(lag-float(ref["t_lag_s"]))/float(ref["t_lag_s"]),
              "P_relative_to_series": abs(p-p_series)/p_series,
              "bath_max_mass_relative_error": bath["max_mass_relative_error"]}

    params = read_csv(ROOT / "WP4/parameters.csv")
    for name, task, held in (("NS_six_film_forward", "fullB", ""), ("NS_2600_film_holdout_replay", "film", "2600")):
        t0 = time.perf_counter()
        par = next(r for r in params if r["model"] == "NS" and r["task"] == task and r["held_out_h_nm"] == held)
        predictions, diag = predict(data, nets, spec, float(par["alpha_m_inv"]), float(par["beta"]), float(par["Q_uS_cm_per_unit_uR"]))
        saved = read_csv(ROOT / ("WP4/v7_NS_fullB_points.csv" if task == "fullB" else "WP4/v7_NS_film_points.csv"))
        if task == "film":
            saved = [r for r in saved if r["hold_h_nm"] == "2600"]
        diffs = [abs(float(r["predicted_increment_uS_cm"]) - predictions[int(r["h_nm"] if task == "fullB" else r["hold_h_nm"])][int(r["point_index"])]) for r in saved]
        report["checks"][name] = {"elapsed_s": time.perf_counter()-t0,
                 "points": len(saved), "max_abs_difference_uS_cm": max(diffs),
                 "max_mass_relative_error": max(v["max_mass_relative_error"] for v in diag.values()),
                 "minimum_u": min(v["minimum_u"] for v in diag.values()),
                 "fitted_parameters_from": "WP4/parameters.csv", "refit_performed": False}
        if task == "film":
            report["checks"][name]["endpoint_pred_over_observed"] = predictions[2600][-1]/(data[2600][-1][1]-data[2600][0][1])

    t0 = time.perf_counter()
    figure1()
    generated = ROOT / "data/processed/figures/v7_Fig1_material_observation.csv"
    ref_lines = (TESTS / "reference_Fig1_material_observation.csv").read_text(encoding="utf-8-sig").splitlines()
    actual_lines = generated.read_text(encoding="utf-8-sig").splitlines()
    report["checks"]["figure1"] = {"elapsed_s": time.perf_counter()-t0,
              "source_rows": len(actual_lines)-1, "source_exact_match": actual_lines == ref_lines,
              "svg_present": generated.with_suffix(".svg").is_file(),
              "png_present": generated.with_suffix(".png").is_file()}

    report["elapsed_total_s"] = time.perf_counter()-started
    report["pass"] = (report["checks"]["straight_chain"]["P_relative_to_archived"] < 1e-9
              and report["checks"]["straight_chain"]["lag_relative_to_archived"] < 1e-9
              and report["checks"]["straight_chain"]["bath_max_mass_relative_error"] < 1e-8
              and report["checks"]["NS_six_film_forward"]["max_abs_difference_uS_cm"] < 1e-6
              and report["checks"]["NS_2600_film_holdout_replay"]["max_abs_difference_uS_cm"] < 1e-6
              and report["checks"]["figure1"]["source_exact_match"])
    (TESTS / "representative_reproduction.json").write_text(json.dumps(report, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"pass": report["pass"], "elapsed_total_s": report["elapsed_total_s"],
                      "checks": report["checks"]}, ensure_ascii=False))
    if not report["pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
