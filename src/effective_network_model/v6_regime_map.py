"""WP3 fixed v6 mechanism/regime matrix; synthetic conditional calculations only."""

from __future__ import annotations

import csv
import json
import math
import sys
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "WP2/src"))
from transport import Config, D_B_REF, field, m0_right_fraction, solve, steady_permeance  # noqa: E402


OUTPUT = ROOT / "WP3/02_消融及适用域结果.csv"
RUN = ROOT / "WP3/04_运行记录.json"
ELL_M = 900e-9
A_M = 25e-9
GAMMA_PER_M = 1e-4 / ELL_M
OMEGA = np.logspace(-6, 6, 41)
ETA = np.logspace(math.log10(.05), math.log10(20), 41)
THICKNESSES_NM = (65, 120, 200, 350, 700, 2600)
GRID_CELLS = 50


def make_config(h_m: float, dm_ratio: float, ws: float, chi: float = 0) -> Config:
    dm = D_B_REF * dm_ratio
    return Config(h_m=h_m, a_m=A_M, ell_m=ELL_M, ws_ref=ws,
                  dm_m2_s=dm, db_m2_s=D_B_REF, gamma_per_m=GAMMA_PER_M,
                  volume_ratio=1., robin_r=0., chi=chi,
                  tau_s=ELL_M * ELL_M / dm)


def equilibrium_right(cfg: Config, model: str, n: int = GRID_CELLS) -> float:
    if model == "M0":
        return cfg.volume_ratio / (1 + cfg.volume_ratio)
    x = (np.arange(n) + .5) * cfg.h_m / n
    partition_length = float(np.sum(np.exp(-field(x, cfg))) * cfg.h_m / n)
    return cfg.volume_ratio / (cfg.volume_ratio + 1 + cfg.gamma_per_m * (1 + cfg.chi) * partition_length)


def threshold_time(times: np.ndarray, response: np.ndarray, limit: float) -> float | None:
    good = np.where(response >= limit)[0]
    if len(good) == 0:
        return None
    idx = int(good[0])
    if idx == 0:
        return float(times[0])
    y0, y1 = response[idx - 1], response[idx]
    if y1 == y0:
        return float(times[idx])
    return float(times[idx - 1] + (limit - y0) * (times[idx] - times[idx - 1]) / (y1 - y0))


def one_trajectory(model: str, cfg: Config, times: np.ndarray):
    if model == "M0":
        return m0_right_fraction(times, cfg), None
    tr = solve(cfg, times, GRID_CELLS, rtol=1e-8, atol=1e-11)
    return tr.right_fraction, tr


def one_row(kind: str, model: str, cfg: Config, eta: float, omega: float,
            response: float, tr, idx: int, t50: float | None, eq_right: float) -> dict[str, object]:
    mobile = slow = 0.
    left = cfg.volume_ratio - response if model == "M0" else 0.
    mass_error = 0.
    min_state = 0.
    if tr is not None:
        y = tr.states[:, idx]
        left = cfg.volume_ratio * y[0]
        mobile = cfg.gamma_per_m * tr.operator.dx_m * float(np.sum(y[1:tr.operator.n_cells + 1]))
        if tr.operator.has_storage:
            slow = cfg.gamma_per_m * tr.operator.dx_m * float(np.sum(y[tr.operator.n_cells + 2:]))
        mass_error = tr.mass_error
        min_state = float(np.min(tr.states))
    t = omega * ELL_M * ELL_M / cfg.dm_m2_s
    p = steady_permeance(cfg)
    b_factor = response / (cfg.gamma_per_m * p * t) if response > 0 and t > 0 else None
    return dict(record_type=kind, model_id=model, dm_over_db=cfg.dm_m2_s / D_B_REF,
                ws_ref=cfg.ws_ref, chi=cfg.chi, tau_s=cfg.tau_s if cfg.chi else "",
                h_nm=cfg.h_m * 1e9, eta=eta, omega=omega, time_s=t,
                P_connected_m_s=p, receiver_fraction=response,
                receiver_equilibrium_fraction=eq_right, B_factor=b_factor,
                left_pool_salt_star=left, mobile_salt_star=mobile,
                slow_salt_star=slow, receiver_salt_star=response,
                t50_own_equilibrium_s=t50 if t50 is not None else "not_reached",
                S_ss="", S_T="", dlnB_dlnh="",
                mass_relative_error=mass_error, minimum_state=min_state,
                note="synthetic fixed matrix; not fitted PE parameters")


def add_slopes(rows: list[dict[str, object]]) -> None:
    # For each fixed physical T (fixed Ω within a scenario), eta varies with h.
    for omega_idx in range(len(OMEGA)):
        cohort = rows[omega_idx::len(OMEGA)]
        logh = np.log([float(r["h_nm"]) for r in cohort])
        p = np.array([float(r["P_connected_m_s"]) for r in cohort])
        y = np.array([float(r["receiver_fraction"]) for r in cohort])
        ss = np.gradient(np.log(p), logh)
        # Near-zero receiver values have unstable logarithmic derivatives;
        # retain raw responses but mark slopes unresolved below this numerical
        # reporting floor (not a measurement threshold).
        eq = np.array([float(r["receiver_equilibrium_fraction"]) for r in cohort])
        ok = y > np.maximum(1e-12, 1e-8 * eq)
        st = np.full(len(y), np.nan)
        for i in range(len(y)):
            if not ok[i]:
                continue
            if i == 0 and ok[1]:
                st[i] = (np.log(y[1]) - np.log(y[0])) / (logh[1] - logh[0])
            elif i == len(y) - 1 and ok[-2]:
                st[i] = (np.log(y[-1]) - np.log(y[-2])) / (logh[-1] - logh[-2])
            elif 0 < i < len(y) - 1 and ok[i - 1] and ok[i + 1]:
                st[i] = (np.log(y[i + 1]) - np.log(y[i - 1])) / (logh[i + 1] - logh[i - 1])
        for i, r in enumerate(cohort):
            r["S_ss"] = float(ss[i])
            if np.isfinite(st[i]):
                r["S_T"] = float(st[i])
                r["dlnB_dlnh"] = float(st[i] - ss[i])


def scenario_rows(model: str, dm_ratio: float, chi: float) -> list[dict[str, object]]:
    rows = []
    for eta in ETA:
        cfg = make_config(float(eta * ELL_M), dm_ratio, 7., chi)
        times = np.r_[0., OMEGA * ELL_M * ELL_M / cfg.dm_m2_s]
        response, tr = one_trajectory(model, cfg, times)
        eq_right = equilibrium_right(cfg, model)
        t50 = threshold_time(times, response, .5 * eq_right)
        for j, omega in enumerate(OMEGA):
            rows.append(one_row("regime_map", model, cfg, float(eta), float(omega),
                                float(response[j + 1]), tr, j + 1, t50, eq_right))
    add_slopes(rows)
    return rows


def ablation_rows() -> list[dict[str, object]]:
    rows = []
    for dm_ratio in (1., 1e-8):
        for h_nm in THICKNESSES_NM:
            for ws in (0., 7.):
                for model, chi in (("M0", 0.), ("M1", 0.), ("M2", 10.), ("M2", 1000.)):
                    cfg = make_config(h_nm * 1e-9, dm_ratio, ws, chi)
                    times = np.array([0., 7200.])
                    response, tr = one_trajectory(model, cfg, times)
                    eq_right = equilibrium_right(cfg, model)
                    omega = 7200 * cfg.dm_m2_s / (ELL_M * ELL_M)
                    rows.append(one_row("ablation_120min", model, cfg, h_nm * 1e-9 / ELL_M,
                                        omega, float(response[-1]), tr, 1, None, eq_right))
    return rows


def main() -> None:
    started = datetime.now(timezone.utc).isoformat()
    all_rows = []
    for ratio in (1., 1e-8):
        for model, chi in (("M0", 0.), ("M1", 0.), ("M2", 10.), ("M2", 1000.)):
            print(f"regime {model} Dm/Db={ratio:g} chi={chi:g}", flush=True)
            all_rows.extend(scenario_rows(model, ratio, chi))
    all_rows.extend(ablation_rows())
    assert len(all_rows) == 8 * 41 * 41 + 2 * 6 * 2 * 4
    assert max(float(x["mass_relative_error"]) for x in all_rows) < 1e-8
    assert min(float(x["minimum_state"]) for x in all_rows) >= -1e-10
    with OUTPUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(all_rows[0]))
        writer.writeheader()
        writer.writerows(all_rows)
    RUN.write_text(json.dumps({
        "started_utc": started, "finished_utc": datetime.now(timezone.utc).isoformat(),
        "actual_executor": "Codex on local workspace; not student independent execution",
        "command": ".venv/Scripts/python.exe src/effective_network_model/v6_regime_map.py",
        "companion_artifacts_present": [str(p.relative_to(ROOT)) for p in (
            ROOT / "WP3/量级反算_20与60min.csv",
            ROOT / "WP3/figures/有限时间厚度斜率适用图.png") if p.is_file()],
        "solver": "WP2/src/transport.py", "source_specification": "导师 v6.0.0 04 section 5",
        "grid_eta": [float(ETA[0]), float(ETA[-1]), len(ETA)],
        "grid_omega": [float(OMEGA[0]), float(OMEGA[-1]), len(OMEGA)],
        "grid_cells": GRID_CELLS, "rows": len(all_rows),
        "maximum_mass_relative_error": max(float(x["mass_relative_error"]) for x in all_rows),
        "minimum_state": min(float(x["minimum_state"]) for x in all_rows),
        "scope": "synthetic fixed matrix; no PE fit or measured material map",
    }, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"{OUTPUT}: {len(all_rows)} rows")


if __name__ == "__main__":
    main()
