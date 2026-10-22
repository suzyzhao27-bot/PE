"""Run v6 WP2 technical checks N1–N9 on synthetic, dimensionless states."""

from __future__ import annotations

import csv
import hashlib
import json
import platform
import sys
from dataclasses import replace
from datetime import datetime, timezone
from math import exp, pi
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import least_squares

from transport import (Config, D_B_REF, build_operator, field, m0_right_fraction,
                       minutes_to_seconds, solve, steady_permeance)


ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = ROOT / "WP2/configs/technical_validation.json"
TABLE = ROOT / "WP2/02_数值验证表.csv"
RECORD = ROOT / "WP2/03_运行记录.json"
BOOK1 = ROOT / "data/processed/origin_exports/Polymer_salt_2h/Book1__Sheet1.csv"
RAW_OPJU = ROOT / "原始资料/赵珈绮/01_研究框架与资料_v2/05_真实数据_v1核查/source_snapshots/Polymer_salt_2h.opju"
ROWS: list[dict[str, object]] = []


def record(case: str, metric: str, measured: float, limit: float,
           *, comparator: str = "<=", note: str = "", diagnostic: bool = False) -> None:
    passed = (measured <= limit) if comparator == "<=" else (measured >= limit)
    ROWS.append(dict(case=case, metric=metric, measured=f"{measured:.12g}",
                     limit=f"{limit:.12g}", comparator=comparator,
                     status="OBSERVED" if diagnostic else ("PASS" if passed else "FAIL"), note=note))
    if not passed and not diagnostic:
        print(f"FAIL {case} {metric}: {measured} {comparator} {limit}")


def right_flux(tr) -> np.ndarray:
    return np.array([tr.operator.fluxes(tr.states[:, i])[1]
                     for i in range(tr.states.shape[1])])


def breakthrough_ratio(t: float, h: float, dm: float) -> float:
    n = np.arange(1, 100)
    return float(1 + 2 * np.sum((-1.) ** n * np.exp(-n * n * pi * pi * dm * t / (h * h))))


def lag_from_trace(tr, cfg: Config) -> float:
    return float(tr.times_s[-1] - np.trapezoid(right_flux(tr), tr.times_s) / steady_permeance(cfg))


def n1_n2_n3() -> None:
    base = Config(h_m=200e-9, ws_ref=0, dm_m2_s=1e-10,
                  access=False, fixed_pools=True)
    t = np.linspace(0, 0.004, 2001)
    tr = solve(base, t, 100)
    expected_lag = base.h_m ** 2 / (6 * base.dm_m2_s)
    record("N1", "relative_time_lag_error", abs(lag_from_trace(tr, base) / expected_lag - 1), 0.001)
    for time_s in (0.0001, 0.0002, 0.0004):
        idx = int(round(time_s / t[1]))
        flux_ratio = right_flux(tr)[idx] / steady_permeance(base)
        exact = breakthrough_ratio(time_s, base.h_m, base.dm_m2_s)
        record("N1", f"flux_relative_error_{time_s:g}s", abs(flux_ratio / exact - 1), 0.001)

    # ell >> h makes w(x) constant, with actual w0 = 2*ws_ref for r=0.
    ratios = []
    permeances = []
    for w0 in (0., 3., 7.):
        cfg = replace(base, ell_m=1e20, ws_ref=w0 / 2)
        curve = solve(cfg, np.array([0., 0.0001, 0.0004, 0.004]), 100)
        ratios.append(right_flux(curve)[1:3] / steady_permeance(cfg))
        permeances.append(steady_permeance(cfg))
    for w0, shape, p in zip((3., 7.), ratios[1:], permeances[1:]):
        record("N2", f"breakthrough_shape_w{w0:g}", float(np.max(np.abs(shape - ratios[0]))), 0.001)
        record("N2", f"permeance_ratio_w{w0:g}", abs(p / permeances[0] / exp(-w0) - 1), 0.001)

    cfg = replace(base, h_m=700e-9, ws_ref=7, ell_m=900e-9, access=True)
    curve = solve(cfg, np.array([0., 0.005, 0.05]), 100)
    record("N3", "spatial_field_steady_flux_relative_error",
           abs(right_flux(curve)[-1] / steady_permeance(cfg) - 1), 0.001)
    # When membrane storage is negligible against bath volume, M1 recovers M0.
    closed = replace(base, gamma_per_m=1e3, fixed_pools=False)
    times = np.array([0., 1., 10., 100.])
    dynamic = solve(closed, times, 100).right_fraction
    quasi = m0_right_fraction(times, closed)
    record("N3", "M0_M1_small_storage_relative_error",
           float(np.max(np.abs(dynamic[1:] - quasi[1:]) / quasi[1:])), 0.001)


def n4_n5_n6() -> None:
    cfg = Config(h_m=200e-9, ws_ref=3, dm_m2_s=1e-12,
                 gamma_per_m=1e6, chi=2, tau_s=0.01)
    op = build_operator(cfg, 50)
    mobile_eq = 0.7 * np.exp(-op.field_centers)
    y_eq = op.initial(0.7, 0.7, mobile_eq, cfg.chi * mobile_eq)
    flux_scale = cfg.dm_m2_s / cfg.h_m
    record("N4", "equilibrium_flux_dimensionless",
           max(abs(x) for x in op.fluxes(y_eq)) / flux_scale, 1e-8)
    record("N4", "equilibrium_derivative_absolute_per_s",
           float(np.max(np.abs(op.matrix @ y_eq))), 1e-8)
    mobile = 0.1 * np.exp(-op.field_centers)
    t = np.array([0., 0.001, 0.01, 0.1, 1., 10., 100.])
    tr = solve(cfg, t, 50, left=1, right=0.1, mobile=mobile, storage=cfg.chi * mobile)
    record("N4", "closed_mass_relative_error", tr.mass_error, 1e-8)
    record("N4", "minimum_concentration", float(np.min(tr.states)), -1e-11, comparator=">=")

    f = [tr.operator.free_energy(tr.states[:, i]) for i in range(len(t))]
    record("N6", "max_positive_free_energy_step",
           max(0., float(np.max(np.diff(f)))), 1e-8)
    record("N6", "equilibrium_total_free_energy_defined", abs(op.free_energy(y_eq)), 1e6)

    base = Config(h_m=200e-9, ws_ref=0, dm_m2_s=1e-10,
                  access=False, fixed_pools=True)
    times = np.linspace(0, 0.004, 2001)
    m1 = solve(base, times, 100)
    zero = solve(replace(base, chi=0, tau_s=1e-8), times, 100)
    record("N5", "chi_zero_exact_state_difference",
           float(np.max(np.abs(zero.states - m1.states))), 1e-12)
    fast_cfg = replace(base, chi=2, tau_s=1e-8)
    fast = solve(fast_cfg, times, 100)
    expected = 3 * base.h_m ** 2 / (6 * base.dm_m2_s)
    record("N5", "fast_exchange_lag_relative_error", abs(lag_from_trace(fast, fast_cfg) / expected - 1), 0.001)
    slow = solve(replace(base, chi=2, tau_s=1e6), times, 100)
    record("N5", "slow_exchange_lag_relative_to_M1",
           abs(lag_from_trace(slow, base) / lag_from_trace(m1, base) - 1), 0.001)


def n7_n8() -> None:
    base = Config(h_m=200e-9, ws_ref=2, dm_m2_s=D_B_REF * 1e-4,
                  gamma_per_m=1e5, q_conductivity=200)
    for h_nm in (65., 200., 2600.):
        cfg = base.with_nm(h_nm)
        t_char = cfg.h_m ** 2 / cfg.dm_m2_s
        values = []
        tight = None
        for n in (50, 100, 200):
            sample_times = np.array([0., 0.25 * t_char, t_char, 10 * t_char])
            tr = solve(cfg, sample_times, n)
            record("N7", f"mass_h{h_nm:g}_nx{n}", tr.mass_error, 1e-8)
            values.append(tr.right_fraction[1:])
            if n == 200:
                tight = solve(cfg, sample_times, n, rtol=1e-10, atol=1e-13)
        scale = np.maximum(np.abs(values[-1]), 1e-8)
        record("N7", f"grid_relative_h{h_nm:g}",
               float(np.max(np.abs(values[0] - values[-1]) / scale)), 0.01)
        record("N7", f"time_tolerance_relative_h{h_nm:g}",
               float(np.max(np.abs(values[-1] - tight.right_fraction[1:]) / scale)), 0.01)

    m2cfg = replace(base, chi=10, tau_s=100.)
    values = []
    for n in (50, 100, 200):
        tr = solve(m2cfg, np.array([0., 1., 10., 100.]), n)
        values.append(tr.right_fraction[1:])
        record("N7", f"M2_mass_nx{n}", tr.mass_error, 1e-8)
    scale = np.maximum(np.abs(values[-1]), 1e-8)
    record("N7", "M2_grid_relative", float(np.max(np.abs(values[0] - values[-1]) / scale)), 0.01)
    tight = solve(m2cfg, np.array([0., 1., 10., 100.]), 200,
                  rtol=1e-10, atol=1e-13)
    record("N7", "M2_time_tolerance_relative",
           float(np.max(np.abs(values[-1] - tight.right_fraction[1:]) / scale)), 0.01)

    record("N8", "nm_to_m_error", abs(base.with_nm(200).h_m - 200e-9), 1e-20)
    record("N8", "minute_to_second_error",
           float(np.max(np.abs(minutes_to_seconds([0, 20, 120]) - [0, 1200, 7200]))), 1e-12)
    t = solve(base, np.array([0., 60.]), 50)
    record("N8", "observation_mapping_error",
           float(np.max(np.abs(t.conductivity_increment - base.q_conductivity * t.right_fraction))), 1e-12,
           note="Q is conditional raw conductivity scale; no measured molar calibration")


def n9() -> None:
    with BOOK1.open(encoding="utf-8-sig", newline="") as stream:
        source_rows = list(csv.DictReader(stream))
    sparse_min = np.array([float(row["col_3_C"]) for row in source_rows], dtype=float)
    assert np.array_equal(sparse_min, [0, 20, 40, 60, 80, 100, 120])
    dense_min = np.arange(0, 121, 5, dtype=float)
    base = Config(h_m=200e-9, ws_ref=3, dm_m2_s=1e-17,
                  gamma_per_m=1e5, q_conductivity=100)
    for name, cfg, symbols, initial, true in (
        ("M1", base, ("log10_dm", "log10_gamma"), (-16.5, 4.5), (-17., 5.)),
        ("M2", replace(base, chi=10, tau_s=1200.),
         ("log10_chi", "log10_tau"), (0.7, 3.4), (1., np.log10(1200.))),
    ):
        target_dense = solve(cfg, minutes_to_seconds(dense_min), 50).conductivity_increment
        target_sparse = solve(cfg, minutes_to_seconds(sparse_min), 50).conductivity_increment

        def predict(theta: np.ndarray, times: np.ndarray) -> np.ndarray:
            if name == "M1":
                trial = replace(cfg, dm_m2_s=10 ** theta[0], gamma_per_m=10 ** theta[1])
            else:
                trial = replace(cfg, chi=10 ** theta[0], tau_s=10 ** theta[1])
            return solve(trial, minutes_to_seconds(times), 50).conductivity_increment

        fit = least_squares(lambda x: predict(x, dense_min)[1:] - target_dense[1:],
                            np.array(initial), diff_step=1e-3,
                            max_nfev=60, xtol=1e-7, ftol=1e-7, gtol=1e-7)
        dense_error = float(np.max(np.abs(predict(fit.x, dense_min) - target_dense)))
        sparse_error = float(np.max(np.abs(predict(fit.x, sparse_min) - target_sparse)))
        record("N9", f"{name}_dense_synthetic_rebuild_absolute", dense_error, 5e-5,
               note=f"fixed other parameters; recovered {symbols}={fit.x.tolist()} true={true}")
        record("N9", f"{name}_real_time_grid_rebuild_absolute", sparse_error, 5e-5,
               note="sampled at Book1 non-65-min time grid; no independent replicate")
        # Sensitivity singular values indicate whether the sparse seven-point grid
        # can distinguish the two varied parameters. This is a local diagnostic.
        x = np.asarray(true)
        delta = 1e-3
        jac = np.column_stack([(predict(x + delta * np.eye(2)[k], sparse_min)[1:]
                                - predict(x - delta * np.eye(2)[k], sparse_min)[1:]) / (2 * delta)
                               for k in range(2)])
        singular = np.linalg.svd(jac, compute_uv=False)
        record("N9", f"{name}_sparse_min_singular_value", float(singular[-1]), 0., comparator=">=",
               note=f"singular values={singular.tolist()}; absolute noise scale not known", diagnostic=True)

        # All nominal shared parameters: a six-point nonzero curve alone cannot
        # uniquely recover seven M2 quantities. Save the local null direction.
        full_names = ["ws", "log10_ell", "log10_dm", "log10_gamma", "log10_Q"]
        full_base = np.array([cfg.ws_ref, np.log10(cfg.ell_m), np.log10(cfg.dm_m2_s),
                              np.log10(cfg.gamma_per_m), np.log10(cfg.q_conductivity)])
        if name == "M2":
            full_names += ["log10_chi", "log10_tau"]
            full_base = np.r_[full_base, np.log10(cfg.chi), np.log10(cfg.tau_s)]

        def full_predict(v: np.ndarray) -> np.ndarray:
            changes = dict(ws_ref=v[0], ell_m=10 ** v[1], dm_m2_s=10 ** v[2],
                           gamma_per_m=10 ** v[3], q_conductivity=10 ** v[4])
            if name == "M2":
                changes.update(chi=10 ** v[5], tau_s=10 ** v[6])
            trial = replace(cfg, **changes)
            y = solve(trial, minutes_to_seconds(sparse_min), 50).conductivity_increment[1:]
            return np.arcsinh(y)

        eye = np.eye(len(full_base))
        full_jac = np.column_stack([
            (full_predict(full_base + delta * eye[k]) - full_predict(full_base - delta * eye[k])) / (2 * delta)
            for k in range(len(full_base))
        ])
        _u, sv, vh = np.linalg.svd(full_jac, full_matrices=True)
        record("N9", f"{name}_full_sparse_sensitivity_singular_min", float(sv[-1]), 0., comparator=">=",
               note=f"names={full_names}; singular_values={sv.tolist()}; weakest_direction={vh[-1].tolist()}; local analysis only", diagnostic=True)


def main() -> int:
    specification = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    assert abs(specification["reference"]["db_m2_s"] - D_B_REF) < 1e-23
    if RECORD.is_file():
        previous = json.loads(RECORD.read_text(encoding="utf-8"))
        raw_source_sha256 = previous["raw_source_sha256"]
        raw_source_hash_status = "reused recorded provenance; not recomputed in this run"
    else:
        raw_source_sha256 = hashlib.sha256(RAW_OPJU.read_bytes()).hexdigest()
        raw_source_hash_status = "one-time source provenance check"
    started = datetime.now(timezone.utc).isoformat()
    n1_n2_n3()
    n4_n5_n6()
    n7_n8()
    n9()
    with TABLE.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(ROWS[0]))
        writer.writeheader()
        writer.writerows(ROWS)
    summary = {status: sum(x["status"] == status for x in ROWS)
               for status in ("PASS", "OBSERVED", "FAIL")}
    record_data = {
        "started_utc": started,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "actual_executor": "Codex on local workspace; not student independent execution",
        "command": ".venv/Scripts/python.exe WP2/src/validate.py",
        "python": sys.version,
        "numpy": np.__version__, "scipy": scipy.__version__,
        "configuration": str(CONFIG_PATH.relative_to(ROOT)),
        "sample_time_source": str(BOOK1.relative_to(ROOT)),
        "raw_source": str(RAW_OPJU.relative_to(ROOT)),
        "raw_source_sha256": raw_source_sha256,
        "raw_source_hash_status": raw_source_hash_status,
        "validation_table": str(TABLE.relative_to(ROOT)),
        "summary": summary,
        "scope": "synthetic technical validation; no PE fit or inference",
        "limitations": [
            "N7 <1% grid/tolerance criterion checked; smaller-than-model-improvement criterion awaits WP4 fitted comparison",
            "N9 recovers two selected parameters with all other values fixed; full five/seven-parameter identifiability belongs to WP5",
        ],
    }
    RECORD.write_text(json.dumps(record_data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"WP2 technical checks: {summary}")
    return 0 if summary["FAIL"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
