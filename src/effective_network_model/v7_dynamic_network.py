"""Frozen ND conductance activation on the same fixed-capacity NS network."""

from __future__ import annotations

import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import coo_matrix, diags

from v7_network_transport import closed_baths, system


def dynamic_system(graph, beta, alpha, tau0_s, b, spec, access=True):
    if tau0_s < 0 or b < 0:
        raise ValueError("tau0_s and b must be nonnegative")
    capacities, edges, area = system(graph, beta, alpha, spec, access=access)
    if len(edges) != len(graph["edges"]):
        raise RuntimeError("static and dynamic edge lists differ")
    a_ref = spec["physics"]["a_ref_m"]
    radii = np.array([edge["a_m"] for edge in graph["edges"]], dtype=float)
    tau = tau0_s * np.exp(b * (a_ref / radii - 1.0))
    i = np.array([edge[0] for edge in edges], dtype=int)
    j = np.array([edge[1] for edge in edges], dtype=int)
    g_inf = np.array([edge[2] for edge in edges], dtype=float)
    if np.any(radii <= 0) or np.any(~np.isfinite(tau)) or np.any(tau < 0):
        raise ValueError("invalid dynamic edge radius or timescale")
    return capacities, edges, area, (i, j, g_inf, tau)


def activation(t_s, tau_s, f0=0.0):
    if not 0 <= f0 <= 1 or t_s < 0:
        raise ValueError("f0 must be in [0,1] and t_s nonnegative")
    tau = np.asarray(tau_s, dtype=float)
    if np.any(tau < 0):
        raise ValueError("tau_s must be nonnegative")
    f = np.ones_like(tau)
    positive = tau > 0
    f[positive] = f0 + (1.0 - f0) * (-np.expm1(-t_s / tau[positive]))
    return f


def effective_time_b0(times_s, tau0_s, f0=0.0):
    times = np.asarray(times_s, dtype=float)
    if tau0_s == 0 or f0 == 1:
        return times.copy()
    return times - (1.0 - f0) * tau0_s * (-np.expm1(-times / tau0_s))


def dynamic_closed_baths(capacities, edges, dynamic_edges, times_s,
                         f0=0.0, initial=None, rtol=1e-7, atol=1e-11):
    times = np.asarray(times_s, dtype=float)
    if times.ndim != 1 or len(times) < 2 or times[0] != 0 or np.any(np.diff(times) <= 0):
        raise ValueError("times_s must start at zero and strictly increase")
    if np.any(capacities <= 0):
        raise ValueError("capacities must be positive")
    i, j, g_inf, tau = dynamic_edges
    n = len(capacities)
    if initial is None:
        initial = np.array([1.0, 0.0] + [0.0] * (n - 2))
    else:
        initial = np.asarray(initial, dtype=float)
    if len(initial) != n or np.any(initial < 0):
        raise ValueError("invalid initial activities")
    if f0 == 1 or np.all(tau == 0):
        return closed_baths(capacities, edges, times, initial=initial,
                            rtol=rtol, atol=atol)
    inv_cap = 1.0 / capacities
    edge_count = len(g_inf)
    rows = np.concatenate((i, i, j, j))
    cols = np.concatenate((i, j, i, j))

    def coefficients(t):
        return g_inf * activation(t, tau, f0)

    def rhs(t, u):
        g = coefficients(t)
        flux = g * (u[j] - u[i])
        change = np.bincount(i, weights=flux, minlength=n) - np.bincount(
            j, weights=flux, minlength=n)
        return change * inv_cap

    def jac(t, _u):
        g = coefficients(t)
        values = np.concatenate((-g * inv_cap[i], g * inv_cap[i],
                                 g * inv_cap[j], -g * inv_cap[j]))
        assert len(values) == 4 * edge_count
        return coo_matrix((values, (rows, cols)), shape=(n, n)).tocsc()

    result = solve_ivp(rhs, (0.0, float(times[-1])), initial, t_eval=times,
                       method="BDF", jac=jac, rtol=rtol, atol=atol)
    if not result.success:
        raise RuntimeError(result.message)
    salt = capacities @ result.y
    if salt[0] <= 0:
        raise ValueError("initial salt must be positive")
    return {"u_R": result.y[1],
            "mass_relative_error": float(np.max(np.abs(salt - salt[0])) / salt[0]),
            "minimum_u": float(np.min(result.y)), "all_u": result.y}


def free_energy(capacities, all_u):
    u = np.asarray(all_u)
    if np.min(u) < -1e-12:
        raise ValueError("negative activities exceed rounding allowance")
    positive = np.maximum(u, 0.0)
    terms = np.zeros_like(positive)
    mask = positive > 0
    terms[mask] = positive[mask] * np.log(positive[mask]) - positive[mask]
    return capacities @ terms
