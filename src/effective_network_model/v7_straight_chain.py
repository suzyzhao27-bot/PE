"""v7 T4 straight-channel reference and one closed two-bath conservation case.

This is a reference benchmark, not the N0/NS random network solver.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve


ROOT = Path(__file__).resolve().parents[2]
SPEC = ROOT / "原始资料/赵珈绮/14_论文补强学生任务书_v7.0.0/configs/frozen_spec.json"
OUT = ROOT / "WP2"


def straight_channel(h: float, k: float, n: int, access: bool, d: float, radius: float) -> tuple[float, float, float]:
    """Return P, fixed-bath cumulative-flux lag, and exact series P.

    The finite-volume capacity and edge partition are both K. The lag is
    evaluated from the first time moment of the linear transient.
    """
    area = np.pi * radius**2
    dx = h / n
    capacity = k * area * dx
    internal_g = d * k * area / dx
    port_resistance = dx / (2 * d * k * area)
    if access:
        port_resistance += 1 / (4 * d * radius)
    port_g = 1 / port_resistance
    diagonal = np.full(n, -2 * internal_g / capacity)
    diagonal[0] = diagonal[-1] = -(internal_g + port_g) / capacity
    if n == 1:
        diagonal[0] = -2 * port_g / capacity
        sub = np.array([])
    else:
        sub = np.full(n - 1, internal_g / capacity)
    generator = diags([sub, diagonal, sub], [-1, 0, 1], format="csc")
    source = np.zeros(n)
    source[0] = port_g / capacity
    u_steady = spsolve(-generator, source)
    assert np.all(u_steady >= 0) and np.all(u_steady <= 1 + 1e-12)
    j_steady = port_g * u_steady[-1]
    integrated_deficit = spsolve(-generator, u_steady)
    lag = port_g * integrated_deficit[-1] / j_steady
    exact_p = 1 / (h / (d * k) + (2 * area / (4 * d * radius) if access else 0))
    return j_steady / area, lag, exact_p


def closed_two_baths(h: float, k: float, n: int, d: float, radius: float, alpha: float) -> dict:
    area = np.pi * radius**2
    dx = h / n
    left_volume = right_volume = area / alpha
    capacities = np.array([left_volume, *([k * area * dx] * n), right_volume])
    conductances = np.array([2 * d * k * area / dx, *([d * k * area / dx] * (n - 1)), 2 * d * k * area / dx])
    operator = np.zeros((n + 2, n + 2))
    for edge, g in enumerate(conductances):
        operator[edge, edge] -= g / capacities[edge]
        operator[edge, edge + 1] += g / capacities[edge]
        operator[edge + 1, edge] += g / capacities[edge + 1]
        operator[edge + 1, edge + 1] -= g / capacities[edge + 1]
    initial = np.zeros(n + 2)
    initial[0] = 1
    times = np.array([0, 20, 60, 120]) * 60.0
    result = solve_ivp(lambda _t, u: operator @ u, (0, times[-1]), initial, method="BDF",
                       jac=operator, t_eval=times, rtol=1e-9, atol=1e-13)
    assert result.success, result.message
    salt = capacities @ result.y
    relative_error = np.max(np.abs(salt - salt[0])) / salt[0]
    min_concentration = float(result.y.min())
    assert relative_error <= 1e-8, relative_error
    assert min_concentration >= -1e-12, min_concentration
    return {"h_nm": h * 1e9, "K": k, "cells": n, "alpha_m_inv": alpha,
            "times_min": [0, 20, 60, 120], "right_activity": result.y[-1].tolist(),
            "max_mass_relative_error": float(relative_error), "min_u": min_concentration,
            "equilibrium_u": float(salt[0] / sum(capacities)), "solver": "scipy.solve_ivp.BDF"}


def main() -> None:
    spec = json.loads(SPEC.read_text(encoding="utf-8-sig"))
    physics = spec["physics"]
    d_na, d_cl = physics["D_Na_m2_s"], physics["D_Cl_m2_s"]
    d = 2 * d_na * d_cl / (d_na + d_cl)
    radius = physics["a_ref_m"]
    config_hash = hashlib.sha256(SPEC.read_bytes()).hexdigest()
    print(f"frozen config sha256={config_hash}; D_b={d:.12g} m2/s; a={radius:.12g} m; K=1,0.1; h=200,2600 nm; cells=32,64,128")
    rows = []
    for h_nm in (200, 2600):
        h = h_nm * 1e-9
        for k in (1.0, 0.1):
            for access in (False, True):
                for n in spec["numerics"]["straight_chain_cells"]:
                    p, lag, series_p = straight_channel(h, k, n, access, d, radius)
                    p_rel = abs(p - series_p) / series_p
                    assert p_rel < 1e-6, (h_nm, k, access, n, p_rel)
                    continuum_lag = h * h / (6 * d) if not access else None
                    rows.append({"h_nm": h_nm, "K": k, "access": str(access).lower(), "cells": n,
                                 "config_sha256": config_hash,
                                 "D_b_m2_s": d, "P_m_s": p, "series_P_m_s": series_p,
                                 "P_relative_error": p_rel, "t_lag_s": lag,
                                 "continuum_no_access_t_lag_s": continuum_lag if continuum_lag is not None else "",
                                 "t_lag_relative_error": abs(lag - continuum_lag) / continuum_lag if continuum_lag is not None else ""})
    for h_nm in (200, 2600):
        for k in (1.0, 0.1):
            for access in ("false", "true"):
                selected = [r for r in rows if r["h_nm"] == h_nm and r["K"] == k and r["access"] == access]
                selected.sort(key=lambda r: r["cells"])
                if access == "false":
                    assert selected[-1]["t_lag_relative_error"] <= 0.01
                assert abs(selected[-1]["t_lag_s"] - selected[-2]["t_lag_s"]) / selected[-1]["t_lag_s"] <= 0.01
    with (OUT / "v7_straight_chain_benchmark.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    closed = closed_two_baths(200e-9, 0.1, 32, d, radius, 1e3)
    (OUT / "v7_two_bath_conservation.json").write_text(json.dumps(closed, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"T4 reference {len(rows)} cases PASS; max P error {max(r['P_relative_error'] for r in rows):.3g}; "
          f"max finest no-access lag error {max(r['t_lag_relative_error'] for r in rows if r['cells']==128 and r['access']=='false'):.3g}")
    print(f"Closed baths: max relative salt error {closed['max_mass_relative_error']:.3g}; min u {closed['min_u']:.3g}")


if __name__ == "__main__":
    main()
