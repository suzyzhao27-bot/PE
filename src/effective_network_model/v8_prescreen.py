"""Frozen v8.0.1 necessary-bound screen; no network linear solves or fitting."""

import csv
import json
import math
import sys
from collections import deque
from pathlib import Path

import numpy as np

from v7_network_geometry import generate, morphology
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.data_ingest.v8_step1_units import conductivity, DB, C0

PACKAGE = ROOT / "原始资料/赵珈绮/16_物理诊断与接续任务_v8.0.1"
if not PACKAGE.exists():
    PACKAGE = ROOT / "taskbook"
OUT = ROOT / "WP2"
RS = 0.35e-9


def conductance(a, length, port, r_s=RS):
    if a <= r_s:
        return 0.0
    lam = r_s / a
    hindrance = 1 - 2.104 * lam + 2.09 * lam**3 - 0.95 * lam**5
    k = (1 - lam) ** 2 * hindrance
    if not (0 < k <= 1):
        raise AssertionError((a, k))
    throat_resistance = length / (DB * math.pi * a * a * k)
    if port:
        return 1 / (throat_resistance + 1 / (4 * DB * (a - r_s)))
    return 1 / throat_resistance


def path_lower(node_count, live_edges):
    """One explicit series path is a valid lower bound; BFS avoids optimization."""
    source, target = node_count, node_count + 1
    neighbors = [[] for _ in range(node_count + 2)]
    for i, j, g in live_edges:
        neighbors[i].append((j, g))
        neighbors[j].append((i, g))
    parent = {source: (-1, 0.0)}
    queue = deque([source])
    while queue and target not in parent:
        i = queue.popleft()
        for j, g in neighbors[i]:
            if j not in parent:
                parent[j] = (i, g)
                queue.append(j)
    if target not in parent:
        return 0.0
    resistance = 0.0
    i = target
    while i != source:
        previous, g = parent[i]
        resistance += 1 / g
        i = previous
    return 1 / resistance


def screen_graph(graph, geom, q, r_s_nm):
    nodes, edges = graph["nodes"], graph["edges"]
    n = len(nodes)
    h_m = graph["h_nm"] * 1e-9
    s = geom["s_nm"] * 1e-9
    dz = h_m / geom["nz"]
    epsilon = geom["epsilon"]
    nxy = geom["nx"] * geom["ny"]
    cuts = np.zeros(geom["nz"] + 1)
    energy = 0.0
    live = []
    radii = []
    all_radii = []
    for e in edges:
        port = e["port_flag"]
        i, j = e["i"], e["j"]
        if port:
            node = j if i == "L" else i
            raw = q * nodes[node]["r_m"]
            face = s * s
        else:
            raw = q * min(nodes[i]["r_m"], nodes[j]["r_m"])
            zi, zj = i // nxy, j // nxy
            face = s * s if zi != zj else s * dz
        a = min(raw, math.sqrt(epsilon * face / math.pi))
        all_radii.append(a * 1e9)
        g = conductance(a, e["L_m"], port, r_s_nm * 1e-9)
        if g <= 0:
            continue
        radii.append(a * 1e9)
        ii = n if i == "L" else n + 1 if i == "R" else i
        jj = n if j == "L" else n + 1 if j == "R" else j
        live.append((ii, jj, g))
        zi = -0.5 if i == "L" else geom["nz"] - 0.5 if i == "R" else i // nxy
        zj = -0.5 if j == "L" else geom["nz"] - 0.5 if j == "R" else j // nxy
        delta = (zj - zi) / geom["nz"]
        energy += g * delta * delta
        if i == "L":
            cuts[0] += g
        elif j == "R":
            cuts[-1] += g
        elif zi != zj:
            cuts[min(zi, zj) + 1] += g
    cut = float(min(cuts))
    lower = path_lower(n, live)
    upper = 0.0 if lower == 0 else min(cut, energy)
    assert -1e-15 <= lower <= upper * (1 + 1e-9) + 1e-28
    area = (geom["nx"] * s) * (geom["ny"] * s)
    alpha = 100.0
    receiver_c = C0 * (-math.expm1(-2 * alpha * upper / area * 7200)) / 2
    return {
        "cut_upper_m3s": cut, "dirichlet_upper_m3s": energy,
        "path_lower_m3s": lower, "endpoint_upper_uScm": conductivity(receiver_c),
        "bound_method": "min(cut,dirichlet)", "chosen_upper_G_m3s": upper,
        "ion_disconnected": lower == 0,
        "ion_edges": len(live), "a_min_nm": min(radii) if radii else 0,
        "a_max_nm": max(radii) if radii else 0,
        "all_a_p05_nm": float(np.quantile(all_radii, 0.05)),
        "all_a_median_nm": float(np.median(all_radii)),
        "all_a_p95_nm": float(np.quantile(all_radii, 0.95)),
        "ion_edge_fraction": len(live) / len(edges),
    }


def main():
    spec = json.loads((PACKAGE / "configs/frozen_spec_v8.json").read_text(encoding="utf-8"))
    old = json.loads((PACKAGE / "inputs/legacy_v7_spec.json").read_text(encoding="utf-8"))
    frozen = spec["geometry"]
    assert frozen["nxy_main"] == 12 and frozen["p_z"] == 0.65 and frozen["sigma_ln"] == 0.35
    geometry_spec = old.copy()
    geometry_spec["geometry"] = dict(old["geometry"])
    geometry_spec["geometry"].update(nx=12, ny=12, p_z=0.65, sigma_ln=0.35, q=0.5)
    hs = frozen["B_h_nm"]
    seeds = frozen["main_seeds"]
    q_values = spec["models"]["Tcal"]["q_grid"]
    assert len(q_values) == 49 and len(hs) == 6 and len(seeds) == 8
    max_nz = max(math.ceil(h / sum(morphology(h, "S", geometry_spec["geometry"]["anchors"])[:2])) for h in hs)
    rows = []
    for h in hs:
        for seed in seeds:
            graph, meta = generate(h, "S", seed, geometry_spec, max_nz)
            info = {"s_nm": meta["s_nm"], "nz": meta["nz"], "epsilon": meta["epsilon"], "nx": 12, "ny": 12}
            for q in q_values:
                bounds = screen_graph(graph, info, q, 0.35)
                rows.append({"h_nm": h, "seed": seed, "q": q, "r_s_nm": 0.35,
                             **bounds, "geometry_hash": meta["geometry_hash"],
                             "geometry_hash_scope": "v7 q=0.5 base radius/topology; q separately recorded"})
    assert len(rows) == 2352
    groups = {}
    for row in rows:
        groups.setdefault((row["h_nm"], row["q"]), []).append(row)
    for group in groups.values():
        assert len(group) == 8
        mean = sum(r["endpoint_upper_uScm"] for r in group) / 8
        for row in group:
            row["seed_group_mean_uScm"] = mean
    threshold = 1363.25 / 3
    accepted = [q for q in q_values if groups[(2600, q)][0]["seed_group_mean_uScm"] >= threshold]
    accepted_set = set(accepted)
    for row in rows:
        row["excluded"] = row["q"] not in accepted_set
        row["reason"] = "G2: 2600 mean upper below observed/3" if row["excluded"] else "necessary upper not excluded; no proof of fit"
    fields = list(rows[0])
    with (OUT / "预筛全部记录.csv").open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    result = {
        "status": "NOT_STARTED_BY_GATE" if not accepted else "NECESSARY_BOUND_NOT_EXCLUDED",
        "candidate_q_count": len(q_values), "accepted_q_count": len(accepted),
        "accepted_q": accepted, "threshold_uScm": threshold,
        "accepted_q_min": min(accepted) if accepted else None,
        "accepted_q_max": max(accepted) if accepted else None,
        "accepted_all_edge_a_p05_range_nm": [min(r["all_a_p05_nm"] for r in rows if r["q"] in accepted_set),
                                              max(r["all_a_p05_nm"] for r in rows if r["q"] in accepted_set)] if accepted else None,
        "accepted_all_edge_a_median_range_nm": [min(r["all_a_median_nm"] for r in rows if r["q"] in accepted_set),
                                                 max(r["all_a_median_nm"] for r in rows if r["q"] in accepted_set)] if accepted else None,
        "observed_2600_uScm": 1363.25, "alpha_bound_m_inv": 100,
        "seed_count_each": 8, "row_count": len(rows),
        "max_2600_mean_upper_uScm": max(groups[(2600, q)][0]["seed_group_mean_uScm"] for q in q_values),
        "new_network_linear_solves": 0,
        "interpretation": "Upper bound is necessary only; no upper-bound ratio used as a P ratio. S anchors are not paired wet-state salt membrane geometry.",
    }
    (OUT / "预筛可行集.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"prescreen rows={len(rows)}; accepted q={len(accepted)}/{len(q_values)}; maximum 2600 mean upper={result['max_2600_mean_upper_uScm']:.6g} uS/cm")


if __name__ == "__main__":
    main()
