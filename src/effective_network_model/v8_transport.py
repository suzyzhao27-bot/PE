"""Frozen v8.0.1 local salt edges, steady network and conservative baths."""

import math
from collections import deque

import numpy as np
from scipy.linalg import expm
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import cg, spsolve

from v8_prescreen import DB, RS, conductance


def weighted_graph(graph, geom, q, model, beta=0.0, r_s=RS):
    nodes = graph["nodes"]
    n = len(nodes)
    s = geom["s_nm"] * 1e-9
    dz = graph["h_nm"] * 1e-9 / geom["nz"]
    epsilon = geom["epsilon"]
    nxy = geom["nx"] * geom["ny"]
    weighted = []
    for e in graph["edges"]:
        i, j = e["i"], e["j"]
        port = e["port_flag"]
        if port:
            node = j if i == "L" else i
            raw = q * nodes[node]["r_m"]
            face = s * s
        else:
            raw = q * min(nodes[i]["r_m"], nodes[j]["r_m"])
            face = s * s if i // nxy != j // nxy else s * dz
        a = min(raw, math.sqrt(epsilon * face / math.pi))
        if model == "Tcal":
            g = conductance(a, e["L_m"], port, r_s)
        elif model in ("Ncal", "Scal"):
            k = math.exp(-beta * 25e-9 / a)
            throat = e["L_m"] / (DB * math.pi * a * a * k)
            g = 1 / (throat + (1 / (4 * DB * a) if port else 0))
        else:
            raise ValueError(model)
        if g > 0:
            ii = n if i == "L" else n + 1 if i == "R" else i
            jj = n if j == "L" else n + 1 if j == "R" else j
            weighted.append((ii, jj, g))
    area = (geom["nx"] * s) * (geom["ny"] * s)
    return weighted, area


def steady(n, weighted, area, rtol=1e-10, solver="direct"):
    """Return conductance and active-node solution. Reservoirs are n and n+1."""
    left, right = n, n + 1
    neighbors = [[] for _ in range(n + 2)]
    for i, j, g in weighted:
        neighbors[i].append(j)
        neighbors[j].append(i)
    def reach(start):
        found = {start}
        queue = deque([start])
        while queue:
            for j in neighbors[queue.popleft()]:
                if j not in found:
                    found.add(j)
                    queue.append(j)
        return found
    connected = reach(left)
    if right not in connected:
        return {"G_m3s": 0.0, "P_m_s": 0.0, "potentials": {},
                "residual_relative": 0.0, "conservation_relative": 0.0, "spanning_nodes": 0}
    active = sorted((connected & reach(right)) - {left, right})
    remap = {node: pos for pos, node in enumerate(active)}
    rows, cols, values = [], [], []
    rhs = np.zeros(len(active))
    for i, j, g in weighted:
        for node, other in ((i, j), (j, i)):
            if node not in remap:
                continue
            k = remap[node]
            rows.append(k); cols.append(k); values.append(g)
            if other in remap:
                rows.append(k); cols.append(remap[other]); values.append(-g)
            elif other == left:
                rhs[k] += g
    if active:
        operator = coo_matrix((values, (rows, cols)), shape=(len(active), len(active))).tocsc()
        if solver == "direct":
            potentials_array = spsolve(operator, rhs)
        elif solver == "cg":
            diagonal = operator.diagonal()
            from scipy.sparse import diags
            preconditioner = diags(1 / diagonal)
            potentials_array, info = cg(operator, rhs, rtol=rtol, atol=0,
                                        M=preconditioner, maxiter=20 * len(active))
            if info != 0:
                raise RuntimeError(f"cg failed: {info}")
        else:
            raise ValueError(solver)
        linear_residual = float(np.max(np.abs(operator @ potentials_array - rhs)))
        potentials = dict(zip(active, potentials_array))
    else:
        linear_residual = 0.0
        potentials = {}
    potential = lambda i: 1.0 if i == left else 0.0 if i == right else potentials.get(i, 0.0)
    flux_left = sum(g * (1 - potential(j if i == left else i))
                    for i, j, g in weighted if i == left or j == left)
    flux_right = sum(g * potential(j if i == right else i)
                     for i, j, g in weighted if i == right or j == right)
    scale = max(abs(flux_left), abs(flux_right), 1e-300)
    assert abs(flux_left - flux_right) / scale < 1e-7
    node_balance = {i: 0.0 for i in active}
    for i, j, g in weighted:
        if i in node_balance:
            node_balance[i] += g * (potential(i) - potential(j))
        if j in node_balance:
            node_balance[j] += g * (potential(j) - potential(i))
    conservation = max((abs(v) for v in node_balance.values()), default=0.0) / scale
    return {"G_m3s": flux_left, "P_m_s": flux_left / area,
            "potentials": potentials, "residual_relative": linear_residual / scale,
            "conservation_relative": conservation, "spanning_nodes": len(active)}


def baths(p, alpha, tau, times_s, phi=0.05, c0=600.0):
    """Bath activities and mass error; material P is m/s and alpha is m^-1."""
    times = np.asarray(times_s, dtype=float)
    rate = alpha * p
    if rate == 0:
        return {"receiver_c_molm3": np.zeros_like(times), "mass_relative_error": 0.0,
                "minimum_u": 0.0, "states": np.stack((np.ones_like(times), np.zeros_like(times)))}
    if tau == 0:
        receiver = -np.expm1(-2 * rate * times) / 2
        donor = 1 - receiver
        states = np.stack((donor, receiver))
        mass_error = float(np.max(np.abs(donor + receiver - 1)))
    else:
        matrix = np.array([
            [-rate, rate, 0],
            [rate / phi, -rate / phi - (1 - phi) / tau, (1 - phi) / tau],
            [0, phi / tau, -phi / tau],
        ])
        states = np.column_stack([expm(matrix * t) @ np.array([1.0, 0.0, 0.0]) for t in times])
        mass = states[0] + phi * states[1] + (1 - phi) * states[2]
        mass_error = float(np.max(np.abs(mass - 1)))
        receiver = states[2]
    return {"receiver_c_molm3": c0 * receiver, "mass_relative_error": mass_error,
            "minimum_u": float(states.min()), "states": states}
