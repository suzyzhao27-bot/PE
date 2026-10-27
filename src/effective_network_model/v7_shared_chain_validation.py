"""Run the prescribed straight-channel T4 cases through the shared solver."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import spsolve

from v7_network_transport import SPEC_PATH, fixed_bath_p, generator
from v7_straight_chain import straight_channel


ROOT = Path(__file__).resolve().parents[2]


def chain(h: float, k: float, n: int, access: bool, d: float, radius: float):
    area = np.pi * radius**2
    dx = h / n
    capacity = k * area * dx
    bath_capacity = area / 1000.0  # Unused by fixed-bath P or lag.
    capacities = np.array([bath_capacity, bath_capacity] + [capacity] * n)
    internal_g = d * k * area / dx
    port_resistance = dx / (2 * d * k * area)
    if access:
        port_resistance += 1 / (4 * d * radius)
    port_g = 1 / port_resistance
    edges = [(0, 2, port_g)]
    edges.extend((i + 2, i + 3, internal_g) for i in range(n - 1))
    edges.append((n + 1, 1, port_g))
    return capacities, edges, area, port_g


def shared_lag(capacities, edges, port_g: float) -> float:
    # The membrane block receives u_L=1, u_R=0. Its transient deficit has
    # first moment (-A)^-1 u_steady, where A is the membrane generator.
    op = generator(capacities, edges)
    membrane = op[2:, 2:]
    source = np.zeros(len(capacities) - 2)
    source[0] = port_g / capacities[2]
    steady = spsolve(-membrane, source)
    deficit_integral = spsolve(-membrane, steady)
    return float(deficit_integral[-1] / steady[-1])


def main() -> None:
    spec = json.loads(SPEC_PATH.read_text(encoding="utf-8-sig"))
    physics = spec["physics"]
    d_na, d_cl = physics["D_Na_m2_s"], physics["D_Cl_m2_s"]
    d = 2 * d_na * d_cl / (d_na + d_cl)
    radius = physics["a_ref_m"]
    rows = []
    for h_nm in (200, 2600):
        h = h_nm * 1e-9
        for k in (1.0, 0.1):
            for access in (False, True):
                for n in spec["numerics"]["straight_chain_cells"]:
                    capacities, edges, area, port_g = chain(h, k, n, access, d, radius)
                    p = fixed_bath_p(capacities, edges, area)
                    lag = shared_lag(capacities, edges, port_g)
                    ref_p, ref_lag, exact_p = straight_channel(h, k, n, access, d, radius)
                    continuum_lag = h*h/(6*d) if not access else None
                    rows.append({"h_nm": h_nm, "K": k, "access": access, "cells": n,
                                 "shared_P_m_s": p, "series_P_m_s": exact_p,
                                 "shared_lag_s": lag,
                                 "P_rel_series": abs(p-exact_p)/exact_p,
                                 "P_rel_reference": abs(p-ref_p)/ref_p,
                                 "lag_rel_reference": abs(lag-ref_lag)/ref_lag,
                                 "lag_rel_continuum": abs(lag-continuum_lag)/continuum_lag
                                 if continuum_lag is not None else ""})
    assert len(rows) == 24
    assert max(row["P_rel_series"] for row in rows) < 1e-6
    assert max(row["P_rel_reference"] for row in rows) < 1e-6
    assert max(row["lag_rel_reference"] for row in rows) < 1e-6
    assert max(row["lag_rel_continuum"] for row in rows
               if row["cells"] == 128 and not row["access"]) <= 0.01
    for h_nm in (200, 2600):
        for k in (1.0, 0.1):
            for access in (False, True):
                selected = [r for r in rows if (r["h_nm"], r["K"], r["access"])
                            == (h_nm, k, access)]
                selected.sort(key=lambda r: r["cells"])
                assert abs(selected[-1]["shared_lag_s"]-selected[-2]["shared_lag_s"]) / selected[-1]["shared_lag_s"] <= 0.01
                assert abs(selected[-1]["shared_P_m_s"]-selected[-2]["shared_P_m_s"]) / selected[-1]["shared_P_m_s"] <= 0.01
    out = ROOT / "WP2/v7_shared_chain_validation.csv"
    with out.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"T4 shared solver: {len(rows)} cases PASS; max P-series relative error "
          f"{max(r['P_rel_series'] for r in rows):.3g}; max lag-reference relative error "
          f"{max(r['lag_rel_reference'] for r in rows):.3g}")


if __name__ == "__main__":
    main()
