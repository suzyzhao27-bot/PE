"""Frozen F S01-S17 one-factor sensitivity, three NS full-B starts each."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from v7_d1_objective import H_VALUES, ROOT, curves, graphs, starts
from v7_e2_inverse_start import residuals
from v7_network_geometry import generate, layer_count, morphology
from v7_network_transport import SPEC_PATH, closed_baths, fixed_bath_p, system


OUT = ROOT / "WP5"
CASE_ID = "S01"
CASE_KEY = "q"
CASE_VALUE = 0.25
GRAPHS = AUDIT = MAIN_FORWARD = RUNS = EXT_RUNS = REPORT = POINTS = COMPARISON = None
_DATA = None
_NETS = None
_SPEC = None


def configure(case_id):
    global CASE_ID, CASE_KEY, CASE_VALUE, GRAPHS, AUDIT, MAIN_FORWARD, RUNS, EXT_RUNS, REPORT, POINTS, COMPARISON
    if case_id not in ("S01", "S02", "S03", "S04", "S05", "S06", "S07", "S08", "S09", "S10", "S11", "S12", "S13", "S14", "S15", "S16", "S17"):
        raise ValueError("Only frozen cases S01-S17 are implemented")
    CASE_ID = case_id
    CASE_KEY, CASE_VALUE = {
        "S01": ("q", 0.25), "S02": ("q", 1),
        "S03": ("p_z", 0.35), "S04": ("p_z", 1),
        "S05": ("sigma_ln", 0), "S06": ("sigma_ln", 0.7),
        "S07": ("kappa_V", 0.25),
        "S08": ("extrapolation", "linear_in_log_h"),
        "S09": ("plot_anchor_joint_endpoint", "low"),
        "S10": ("plot_anchor_joint_endpoint", "high"),
        "S11": ("spacing_factor", 0.75),
        "S12": ("layers_rule", "round_half_up"),
        "S13": ("D_b_multiplier", 0.7),
        "S14": ("D_b_multiplier", 1.3),
        "S15": ("reservoir_ratio", 2),
        "S16": ("preload_membrane_u", 0.25),
        "S17": ("external_layer_m", 1e-5)}[case_id]
    prefix = f"v7_F_{case_id}_"
    GRAPHS = OUT / (prefix + "networks.jsonl.gz")
    AUDIT = OUT / (prefix + "geometry_audit.csv")
    MAIN_FORWARD = OUT / (prefix + "main_forward.json")
    RUNS = OUT / (prefix + "fit_runs.jsonl")
    EXT_RUNS = OUT / (prefix + "fit_extensions.jsonl")
    REPORT = OUT / (prefix + "report.json")
    POINTS = OUT / (prefix + "point_predictions.csv")
    COMPARISON = OUT / (prefix + "comparison.json")


class BudgetReached(Exception):
    pass


def case_spec():
    spec = json.loads(SPEC_PATH.read_text(encoding="utf-8-sig"))
    case = next(c for c in spec["matrices"]["F"]["cases"] if c["id"] == CASE_ID)
    assert case == {"id": CASE_ID, "key": CASE_KEY, "value": CASE_VALUE,
                    "conditional": None}
    if CASE_KEY == "D_b_multiplier":
        assert "D_b_multiplier" not in spec["physics"]
        spec["physics"]["D_b_multiplier"] = CASE_VALUE
    elif CASE_KEY == "reservoir_ratio":
        assert spec["physics"]["reservoir_ratio"] == 1
        spec["physics"]["reservoir_ratio"] = CASE_VALUE
    elif CASE_KEY == "preload_membrane_u":
        assert spec["physics"]["initial_membrane_u"] == 0
        spec["physics"]["initial_membrane_u"] = CASE_VALUE
    elif CASE_KEY == "external_layer_m":
        assert "external_layer_m" not in spec["physics"]
        spec["physics"]["external_layer_m"] = CASE_VALUE
    elif CASE_KEY == "plot_anchor_joint_endpoint":
        anchors = spec["geometry"]["anchors"]
        plot_anchors = [a for a in anchors if a["source_type"] == "plot_estimate"]
        assert {a["draw"] for a in plot_anchors} == {81, 182, 357}
        for anchor in plot_anchors:
            for value_key, interval_key in (("dp_nm", "dp_interval_nm"),
                                            ("phi_A", "phi_interval")):
                assert anchor[interval_key][0] <= anchor[value_key] <= anchor[interval_key][1]
                anchor[value_key] = anchor[interval_key][0 if CASE_VALUE == "low" else 1]
    else:
        geometry_key = {"extrapolation": "S_extrapolation",
                        "layers_rule": "nz_rule"}.get(CASE_KEY, CASE_KEY)
        assert spec["geometry"][geometry_key] != case["value"]
        spec["geometry"][geometry_key] = case["value"]
    return spec


def observations():
    raw = curves()
    return {h: {"times": np.array([p[0] for p in raw[h]]),
                "observed": np.array([p[1] - raw[h][0][1] for p in raw[h]])}
            for h in H_VALUES}


def reference_parameters():
    with (OUT / "v7_E1_equivalent_selection.csv").open(
            encoding="utf-8-sig", newline="") as stream:
        reference = next(r for r in csv.DictReader(stream)
                         if r["label"] == "fullB_best")
    return np.array([float(reference[k])
                     for k in ("log10_alpha", "log10_Q", "beta")])


def load_nets():
    nets = {h: {} for h in H_VALUES}
    with gzip.open(GRAPHS, "rt", encoding="utf-8") as stream:
        for line in stream:
            g = json.loads(line)
            assert g["series"] == "S"
            nets[g["h_nm"]][g["seed"]] = g
    assert all(set(nets[h]) == set(range(7001, 7009)) for h in H_VALUES)
    return nets


def prepare():
    if GRAPHS.exists() or AUDIT.exists() or MAIN_FORWARD.exists():
        raise RuntimeError(f"{CASE_ID} preparation already recorded")
    spec = case_spec()
    geom = spec["geometry"]
    max_nz = max(layer_count(h, geom["spacing_factor"] *
                             sum(morphology(h, "S", geom["anchors"],
                                            geom["S_extrapolation"])[:2]),
                             geom["nz_rule"])
                 for h in H_VALUES)
    audits = []
    with gzip.open(GRAPHS, "wt", encoding="utf-8", newline="\n") as stream:
        for h in H_VALUES:
            for seed in range(7001, 7009):
                g, a = generate(h, "S", seed, spec, max_nz)
                assert a["volume_relative_error"] < 1e-12
                assert a["left_port_area_fraction"] <= 1 + 1e-12
                assert a["right_port_area_fraction"] <= 1 + 1e-12
                stream.write(json.dumps(g, ensure_ascii=False) + "\n")
                audits.append(a)
    with AUDIT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(audits[0]))
        writer.writeheader()
        writer.writerows(audits)
    assert len(audits) == 48
    print(f"F/{CASE_ID} prepared 48 S graphs; max volume error "
          f"{max(a['volume_relative_error'] for a in audits):.3g}")
    main_forward()


def main_forward():
    if MAIN_FORWARD.exists():
        raise RuntimeError(f"{CASE_ID} main forward already recorded")
    if not GRAPHS.exists() or not AUDIT.exists():
        raise RuntimeError(f"{CASE_ID} geometry preparation is incomplete")
    data = observations()
    x_ref = reference_parameters()
    with ProcessPoolExecutor(max_workers=6, initializer=init_worker,
                             initargs=(CASE_ID,)) as pool:
        predicted, mass, minimum = predict(pool, x_ref)
    r = residuals(data, predicted)
    forward = {"case": CASE_ID, "model": "NS",
               "changed": {CASE_KEY: CASE_VALUE},
               "main_transformed": x_ref.tolist(), "L": float(r @ r),
               "max_mass_relative_error": mass, "minimum_u": minimum,
               "predicted_increment_by_h":
                   {str(h): predicted[h].tolist() for h in H_VALUES}}
    MAIN_FORWARD.write_text(json.dumps(forward, ensure_ascii=False, indent=2) + "\n",
                            encoding="utf-8")
    print(f"F/{CASE_ID} main-parameter forward L={forward['L']:.8g}; mass={mass:.3g}")


def init_worker(case_id):
    global _DATA, _NETS, _SPEC
    configure(case_id)
    _DATA = observations()
    _NETS = load_nets()
    _SPEC = case_spec()


def forward_one(task):
    h, seed, alpha, beta = task
    cap, edges, _ = system(_NETS[h][seed], beta, alpha, _SPEC)
    initial = initial_state(cap, _SPEC)
    result = closed_baths(cap, edges, _DATA[h]["times"] * 60, initial=initial)
    return h, seed, result["u_R"], result["mass_relative_error"], result["minimum_u"]


def initial_state(cap, spec):
    preload = spec["physics"]["initial_membrane_u"]
    if preload == 0:
        return None
    initial = np.zeros(len(cap))
    initial[2:] = preload
    initial[0] = 1 - preload * float(np.sum(cap[2:])) / cap[0]
    if initial[0] < 0:
        raise ValueError("Preload requires negative donor concentration")
    assert math.isclose(float(cap @ initial), float(cap[0]), rel_tol=1e-12)
    return initial


def predict(pool, x):
    alpha, q, beta = 10 ** x[0], 10 ** x[1], x[2]
    group = {h: [] for h in H_VALUES}
    tasks = [(h, seed, alpha, beta) for h in H_VALUES
             for seed in range(7001, 7009)]
    for h, seed, right, mass, minimum in pool.map(forward_one, tasks):
        group[h].append((seed, right, mass, minimum))
    predicted = {}
    for h in H_VALUES:
        rows = sorted(group[h])
        assert [r[0] for r in rows] == list(range(7001, 7009))
        assert max(r[2] for r in rows) <= 1e-8
        assert min(r[3] for r in rows) >= -1e-12
        predicted[h] = q * np.mean(np.vstack([r[1] for r in rows]), axis=0)
    return predicted, max(r[2] for rows in group.values() for r in rows), min(
        r[3] for rows in group.values() for r in rows)


def fit(index, extended=False):
    if index not in range(3):
        raise ValueError("F uses three frozen starts")
    if not MAIN_FORWARD.exists():
        raise RuntimeError("Run prepare first")
    original = ([json.loads(line) for line in RUNS.read_text(
                encoding="utf-8").splitlines()] if RUNS.exists() else [])
    if extended:
        if len(original) != 3:
            raise RuntimeError("All three original starts must finish before extension")
        eligible = sorted((r for r in original
                           if r["status"] in ("max_eval", "max_eval_or_solver_stop")
                           and r["best_L"] is not None),
                          key=lambda r: (r["best_L"], r["start_index"]))[:2]
        if index not in [r["start_index"] for r in eligible]:
            raise RuntimeError("Only two lowest-loss unfinished starts may extend")
        if EXT_RUNS.exists() and any(
                json.loads(line)["start_index"] == index
                for line in EXT_RUNS.read_text(encoding="utf-8").splitlines()):
            raise RuntimeError("Frozen extension already recorded")
        prior = next(r for r in original if r["start_index"] == index)
        x0 = np.array(prior["best_transformed"])
        out = EXT_RUNS
    else:
        if any(r["start_index"] == index for r in original):
            raise RuntimeError(f"{CASE_ID} start already recorded")
        prior = None
        x0 = starts("NS")[index]
        out = RUNS
    data = observations()
    count = 0
    best = {"L": math.inf, "x": None}
    started = time.perf_counter()
    pool = ProcessPoolExecutor(max_workers=6, initializer=init_worker,
                               initargs=(CASE_ID,))

    def objective(x):
        nonlocal count
        if count >= 200:
            raise BudgetReached()
        count += 1
        predicted, _, _ = predict(pool, x)
        r = residuals(data, predicted)
        loss = float(r @ r)
        if loss < best["L"]:
            best.update(L=loss, x=x.tolist())
        if count % 10 == 0:
            phase = "extension" if extended else "original"
            print(f"F/{CASE_ID} {phase} start{index}: calls {count}/200; best L={best['L']:.6g}",
                  flush=True)
        return r

    status, message, optimizer = "", "", {}
    try:
        result = least_squares(objective, x0, bounds=(np.zeros(3), [5, 6, 6]),
                               method="trf", jac="2-point", x_scale=1.0,
                               ftol=1e-8, xtol=1e-8, gtol=1e-8, max_nfev=200)
        optimizer = {"scipy_success": bool(result.success),
                     "scipy_status": int(result.status),
                     "scipy_nfev": int(result.nfev),
                     "message": str(result.message)}
        status = "solver_success" if result.success else "max_eval_or_solver_stop"
        message = str(result.message)
    except BudgetReached:
        status, message = "max_eval", "200 actual objective calls reached"
    except Exception as exc:
        status, message = "FAIL", repr(exc)
    finally:
        pool.shutdown(wait=True)
    record = {"case": CASE_ID, "model": "NS", "start_index": index,
              "phase": "extension" if extended else "original",
              "prior_actual_objective_calls": prior["actual_objective_calls"] if prior else 0,
              "start_transformed": x0.tolist(),
              "bounds_transformed": {"lower": [0, 0, 0], "upper": [5, 6, 6]},
              "actual_objective_call_limit": 200, "actual_objective_calls": count,
              "status": status, "message": message, "optimizer": optimizer,
              "best_L": best["L"] if math.isfinite(best["L"]) else None,
              "best_transformed": best["x"],
              "elapsed_s": time.perf_counter() - started}
    record["cumulative_actual_objective_calls"] = (
        record["prior_actual_objective_calls"] + count)
    with out.open("a", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(record, ensure_ascii=False) + "\n")
    print(f"F/{CASE_ID} {record['phase']} start{index}: {status}; "
          f"calls={count}; best L={record['best_L']}")
    if status == "FAIL":
        raise RuntimeError(message)


def final_runs():
    original = [json.loads(line) for line in RUNS.read_text(encoding="utf-8").splitlines()]
    assert len(original) == 3 and sorted(r["start_index"] for r in original) == [0, 1, 2]
    extensions = ([json.loads(line) for line in EXT_RUNS.read_text(
                  encoding="utf-8").splitlines()] if EXT_RUNS.exists() else [])
    assert len(extensions) <= 2
    assert len({r["start_index"] for r in extensions}) == len(extensions)
    eligible = sorted((r for r in original
                       if r["status"] in ("max_eval", "max_eval_or_solver_stop")
                       and r["best_L"] is not None),
                      key=lambda r: (r["best_L"], r["start_index"]))[:2]
    assert {r["start_index"] for r in extensions} <= {
        r["start_index"] for r in eligible}
    by_index = {r["start_index"]: r for r in original}
    for r in extensions:
        prior = by_index[r["start_index"]]
        assert prior["status"] in ("max_eval", "max_eval_or_solver_stop")
        assert r["prior_actual_objective_calls"] == prior["actual_objective_calls"]
        assert r["cumulative_actual_objective_calls"] <= 400
        by_index[r["start_index"]] = r
    return original, extensions, [by_index[i] for i in range(3)]


def report():
    if REPORT.exists() or POINTS.exists():
        raise RuntimeError(f"{CASE_ID} report already recorded")
    original, extensions, runs = final_runs()
    good = [r for r in runs if r["status"] == "solver_success"]
    assert good
    best = min(good, key=lambda r: r["best_L"])
    if any(r["status"] != "solver_success" and r["best_L"] is not None
           and r["best_L"] < best["best_L"] for r in runs):
        raise RuntimeError("Lower unfinished run requires frozen extension review")
    x_ref = reference_parameters()
    x_fit = np.array(best["best_transformed"])
    data = observations()
    with ProcessPoolExecutor(max_workers=6, initializer=init_worker,
                             initargs=(CASE_ID,)) as pool:
        base, mass_base, min_base = predict(pool, x_ref)
        fitted, mass_fit, min_fit = predict(pool, x_fit)
    losses = {}
    for label, pred in (("main_parameters", base), ("refit", fitted)):
        r = residuals(data, pred)
        losses[label] = float(r @ r)
    assert math.isclose(losses["refit"], best["best_L"], rel_tol=1e-6,
                        abs_tol=1e-8)
    nets, spec = load_nets(), case_spec()
    p_by_h = {}
    points = []
    for h in H_VALUES:
        p_by_h[str(h)] = {}
        for label, x in (("main_parameters", x_ref), ("refit", x_fit)):
            values = []
            for seed in range(7001, 7009):
                cap, edges, area = system(nets[h][seed], x[2], 10 ** x[0], spec)
                values.append(fixed_bath_p(cap, edges, area))
            p_by_h[str(h)][label] = float(np.mean(values))
        for i, t in enumerate(data[h]["times"]):
            points.append({"case": CASE_ID, "h_nm": h, "time_min": t,
                           "observed_increment_uS_cm": data[h]["observed"][i],
                           "main_parameters_increment_uS_cm": base[h][i],
                           "refit_increment_uS_cm": fitted[h][i]})
    with POINTS.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(points[0]))
        writer.writeheader()
        writer.writerows(points)
    report_data = {"case": CASE_ID, "model": "NS",
                   "changed": {CASE_KEY: CASE_VALUE},
                   "status_by_start": {str(r["start_index"]): r["status"] for r in runs},
                   "actual_objective_calls": sum(r["actual_objective_calls"]
                                                for r in original + extensions),
                   "extension_start_indices": [r["start_index"] for r in extensions],
                   "best_start_index": best["start_index"],
                   "main_transformed": x_ref.tolist(),
                   "fit_transformed": x_fit.tolist(),
                   "L": losses, "P_m_s_by_h": p_by_h,
                   "max_mass_relative_error": max(mass_base, mass_fit),
                   "minimum_u": min(min_base, min_fit),
                   "point_count": len(points)}
    REPORT.write_text(json.dumps(report_data, ensure_ascii=False, indent=2) + "\n",
                      encoding="utf-8")
    print(f"F/{CASE_ID} reported: L={losses}; points={len(points)}")


def compare():
    if COMPARISON.exists():
        raise RuntimeError(f"{CASE_ID} comparison already recorded")
    report_data = json.loads(REPORT.read_text(encoding="utf-8"))
    x_ref = reference_parameters()
    spec = json.loads(SPEC_PATH.read_text(encoding="utf-8-sig"))
    base_nets = graphs()
    base_p = {}
    for h in H_VALUES:
        values = []
        for graph in base_nets[h]:
            cap, edges, area = system(graph, x_ref[2], 10 ** x_ref[0], spec)
            values.append(fixed_bath_p(cap, edges, area))
        base_p[str(h)] = float(np.mean(values))
    with (ROOT / "WP4/all_predictions.csv").open(
            encoding="utf-8-sig", newline="") as stream:
        base_points = [r for r in csv.DictReader(stream)
                       if r["model"] == "NS" and r["task"] == "fullB"
                       and float(r["time_min"]) == 120]
    assert len(base_points) == 6
    base_120 = {r["h_nm"]: float(r["predicted_increment_uS_cm"])
                for r in base_points}
    with POINTS.open(encoding="utf-8-sig", newline="") as stream:
        case_points = [r for r in csv.DictReader(stream)
                       if float(r["time_min"]) == 120]
    assert len(case_points) == 6
    by_h = {}
    for row in case_points:
        h = row["h_nm"]
        p = report_data["P_m_s_by_h"][h]
        by_h[h] = {
            "P_main_m_s": base_p[h],
            f"P_{CASE_ID}_main_over_main": p["main_parameters"] / base_p[h],
            f"P_{CASE_ID}_refit_over_main": p["refit"] / base_p[h],
            "y120_main_uS_cm": base_120[h],
            f"y120_{CASE_ID}_main_uS_cm": float(row["main_parameters_increment_uS_cm"]),
            f"y120_{CASE_ID}_refit_uS_cm": float(row["refit_increment_uS_cm"]),
            "y120_observed_uS_cm": float(row["observed_increment_uS_cm"]),
        }
    ratios = {}
    for label, key in (("main", "y120_main_uS_cm"),
                       (f"{CASE_ID}_main", f"y120_{CASE_ID}_main_uS_cm"),
                       (f"{CASE_ID}_refit", f"y120_{CASE_ID}_refit_uS_cm"),
                       ("observed", "y120_observed_uS_cm")):
        ratios[label] = by_h["200"][key] / by_h["350"][key]
    _, _, runs = final_runs()
    good = [r for r in runs if r["status"] == "solver_success"]
    best = min(good, key=lambda r: r["best_L"])
    stable = any(r is not best and abs(r["best_L"] - best["best_L"]) /
                 max(best["best_L"], np.finfo(float).tiny) <= 0.01
                 for r in good)
    result = {"case": CASE_ID, "model": "NS",
              "status": "multi_start_stable" if stable else "local_optimum_unstable",
              "boundary_within_one_percent": [
                  name for i, name in enumerate(("log10_alpha", "log10_Q", "beta"))
                  if min(report_data["fit_transformed"][i],
                         [5, 6, 6][i] - report_data["fit_transformed"][i])
                  < 0.01 * [5, 6, 6][i]],
              "base_main_L": float(reference_L()),
              f"{CASE_ID}_L": report_data["L"], "by_h": by_h,
              "ratio_200_over_350_at_120": ratios}
    COMPARISON.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n",
                          encoding="utf-8")
    print(f"F/{CASE_ID} comparison: {result['status']}; 200/350={ratios}")


def reference_L():
    with (OUT / "v7_E1_equivalent_selection.csv").open(
            encoding="utf-8-sig", newline="") as stream:
        return next(float(r["L_recomputed"]) for r in csv.DictReader(stream)
                    if r["label"] == "fullB_best")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", choices=("S01", "S02", "S03", "S04", "S05", "S06", "S07", "S08", "S09", "S10", "S11", "S12", "S13", "S14", "S15", "S16", "S17"),
                        default="S01")
    parser.add_argument("action", choices=("prepare", "forward", "fit", "extend", "report", "compare"))
    parser.add_argument("start_index", type=int, nargs="?")
    args = parser.parse_args()
    configure(args.case)
    if args.action == "prepare":
        prepare()
    elif args.action == "forward":
        main_forward()
    elif args.action in ("fit", "extend"):
        if args.start_index is None:
            parser.error("fit requires start_index")
        fit(args.start_index, extended=args.action == "extend")
    elif args.action == "report":
        report()
    else:
        compare()


if __name__ == "__main__":
    main()
