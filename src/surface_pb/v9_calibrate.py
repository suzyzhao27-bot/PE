"""Frozen seven-Gamma, two-start 182-only v9 calibration with actual-call ledger."""

import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from v9_geometry import family_geometry
from v9_observation import solve_velocity
from v9_pb_1d import Chemical, EPS0, F, Z
from v9_pb_2d import solve_equilibrium

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / "原始资料/赵珈绮/18_薄化荷电主线任务书_v9.0.0"
if not PACKAGE.is_dir():
    PACKAGE = ROOT / "taskbook"
OUT = ROOT / "WP4"
CALLS = OUT / "v9_校准目标调用.jsonl"
STARTS = OUT / "v9_校准双起点.jsonl"


class CallCap(Exception):
    pass


def append(path, data):
    with path.open("a", encoding="utf-8") as stream:
        stream.write(json.dumps(data, ensure_ascii=False) + "\n")


def ledger(path):
    if not path.exists():
        return []
    with path.open(encoding="utf-8") as stream:
        return [json.loads(line) for line in stream if line.strip()]


def project_gradient(y, gradient):
    value = gradient.copy()
    value[(y <= 1e-8) & (value > 0)] = 0
    value[(y >= 1 - 1e-8) & (value < 0)] = 0
    return float(np.max(np.abs(value)))


def main():
    if json.loads((ROOT / "WP3/v9_技术验收.json").read_text(encoding="utf-8"))["overall"] != "TECHNICAL_PASS":
        raise RuntimeError("frozen ten technical gates must pass before calibration")
    spec = json.loads((PACKAGE / "configs/frozen_spec_v9.json").read_text(encoding="utf-8"))
    chem_spec = spec["chemistry"]
    with (ROOT / "WP1/v9_九组观测与条件盐度.csv").open(encoding="utf-8-sig", newline="") as stream:
        observations = {r["case_id"]: r for r in csv.DictReader(stream)}
    rows = [observations[f"182_pH{pH}"] for pH in (3, 7, 10)]
    geometry = family_geometry(spec, 182, "backed")
    velocity_cache = {}
    limits = np.asarray([[2.0, -4.0], [6.0, 3.0]])
    span = limits[1] - limits[0]
    completed = {(r["gamma_nm2"], r["start_index"]) for r in ledger(STARTS)}
    prior_calls = ledger(CALLS)
    cap = int(chem_spec["max_actual_objective_calls_per_start"])

    def objective(y, gamma, start_index):
        run_id = f"{gamma:g}_{start_index}"
        old = [r for r in prior_calls if r["run_id"] == run_id]
        count = len(old) + sum(1 for r in active_calls if r["run_id"] == run_id)
        if count >= cap:
            raise CallCap(f"actual objective call cap {cap} reached")
        physical = limits[0] + span * y
        chemical = Chemical(float(gamma), float(physical[0]), float(physical[1]))
        predictions = []
        diagnostics = []
        try:
            for row in rows:
                case_id = row["case_id"]
                field = solve_equilibrium(geometry, float(row["conditional_NaCl_mM"]),
                                          float(row["pH_measured_mean"]), chemical=chemical)
                if case_id not in velocity_cache:
                    velocity_cache[case_id] = solve_velocity(
                        geometry, field, float(row["H_full_gap_um_mean_assumed"]))
                flow = velocity_cache[case_id]
                rho_mobile = F * (Z @ field.concentrations_mol_m3)
                current = float((rho_mobile * flow.velocity_m_s) @ field.mesh.volumes_m2)
                zeta = -flow.H_um**2 * 1e-12 / (12 * EPS0 * 78.3) * current / flow.Q_half_m3_s * 1000
                predictions.append(zeta)
                diagnostics.append({"case_id": case_id, "Gauss_relative": field.Gauss_relative,
                                    "max_species_M": field.max_species_M,
                                    "physical_PB_warning": field.physical_PB_warning})
            residual = [(pred - float(row["zeta_reported_mean_mV"])) / float(row["comparison_scale_mV"])
                        for pred, row in zip(predictions, rows)]
            record = {"run_id": run_id, "gamma_nm2": gamma, "start_index": start_index,
                      "call_number": count + 1, "scaled_y": list(map(float, y)),
                      "pI": float(physical[0]), "b": float(physical[1]),
                      "zeta_pred_mV": predictions, "standardized_residual": residual,
                      "E_cal": float(np.sqrt(np.mean(np.square(residual)))),
                      "numeric_status": "VALID" if all(d["Gauss_relative"] <= 1e-6
                                                        and not d["physical_PB_warning"] for d in diagnostics)
                                        else "PHYSICAL_OR_GAUSS_RISK",
                      "diagnostics": diagnostics}
            append(CALLS, record)
            active_calls.append(record)
            return np.asarray(residual)
        except Exception as exc:
            record = {"run_id": run_id, "gamma_nm2": gamma, "start_index": start_index,
                      "call_number": count + 1, "scaled_y": list(map(float, y)),
                      "pI": float(physical[0]), "b": float(physical[1]),
                      "numeric_status": "NUMERICALLY_UNRESOLVED",
                      "error": f"{type(exc).__name__}: {exc}"}
            append(CALLS, record)
            active_calls.append(record)
            raise

    for gamma in chem_spec["Gamma_nm2_grid"]:
        for start_index, start in enumerate(chem_spec["starts"], 1):
            if (gamma, start_index) in completed:
                continue
            run_id = f"{gamma:g}_{start_index}"
            old = [r for r in prior_calls if r["run_id"] == run_id]
            valid_old = [r for r in old if r.get("numeric_status") == "VALID"]
            begin = min(valid_old, key=lambda r: r["E_cal"])["scaled_y"] if valid_old else (np.asarray(start) - limits[0]) / span
            active_calls = []
            scipy_status = "NOT_RETURNED"
            gradient = None
            reason = None
            try:
                result = least_squares(lambda y: objective(y, gamma, start_index), begin,
                                       bounds=(np.zeros(2), np.ones(2)), method="trf",
                                       jac="2-point", diff_step=1e-3, x_scale=1.0,
                                       ftol=1e-10, xtol=1e-10, gtol=1e-10, max_nfev=cap)
                scipy_status = str(result.message)
                gradient = project_gradient(result.x, result.jac.T @ result.fun / 3)
            except CallCap as exc:
                reason = str(exc)
            except Exception as exc:
                reason = f"{type(exc).__name__}: {exc}"
            all_calls = old + active_calls
            valid = [r for r in all_calls if r.get("numeric_status") == "VALID"]
            best = min(valid, key=lambda r: r["E_cal"]) if valid else None
            status = ("STABLE_LOCAL" if gradient is not None and gradient <= 1e-3 and best
                      else "OPTIMIZATION_UNRESOLVED")
            if reason and not reason.startswith("actual objective call cap"):
                status = "NUMERICALLY_UNRESOLVED"
            summary = {"gamma_nm2": gamma, "start_index": start_index, "declared_start": start,
                       "resumed_from_prior_calls": bool(old), "actual_objective_calls": len(all_calls),
                       "PB_field_upper_bound_used": 3 * len(all_calls), "numeric_diff_step_scaled": 1e-3,
                       "scaled_projected_gradient_inf": gradient, "status": status,
                       "termination_reason": reason or scipy_status,
                       "best_valid": {k: best[k] for k in ("pI", "b", "E_cal", "zeta_pred_mV", "scaled_y")}
                                     if best else None,
                       "best_call_number": best["call_number"] if best else None,
                       "boundary_hit": bool(best and (min(best["scaled_y"]) <= 1e-6
                                                      or max(best["scaled_y"]) >= 1 - 1e-6))}
            append(STARTS, summary)
            completed.add((gamma, start_index))
            prior_calls.extend(active_calls)
            print(f"Γ={gamma:g} start={start_index}: {status}, calls={len(all_calls)}, "
                  f"best_E={best['E_cal'] if best else 'NA'}", flush=True)
    print("CALIBRATION_STARTS_RECORDED")


if __name__ == "__main__":
    main()
