"""N0/NS salt transport on a saved v7 geometry; no fitting or ND branch."""

from __future__ import annotations

import gzip
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import spsolve


ROOT = Path(__file__).resolve().parents[2]
SPEC_PATH = ROOT / "原始资料/赵珈绮/14_论文补强学生任务书_v7.0.0/configs/frozen_spec.json"
GRAPHS = ROOT / "WP2/v7_network_tables.jsonl.gz"


def load_graph(series: str, h_nm: int, seed: int) -> dict:
    with gzip.open(GRAPHS, "rt", encoding="utf-8") as stream:
        for line in stream:
            graph = json.loads(line)
            if (graph["series"], graph["h_nm"], graph["seed"]) == (series, h_nm, seed):
                return graph
    raise KeyError((series, h_nm, seed))


def system(graph: dict, beta: float, alpha: float, spec: dict, access: bool = True):
    physics = spec["physics"]
    d_na, d_cl = physics["D_Na_m2_s"], physics["D_Cl_m2_s"]
    d = 2*d_na*d_cl/(d_na+d_cl) * physics.get("D_b_multiplier", 1)
    a_ref = physics["a_ref_m"]
    count = len(graph["nodes"])
    # The saved node x/y coordinates define the actual lateral cell size.
    xs = sorted({node["x_m"] for node in graph["nodes"]})
    ys = sorted({node["y_m"] for node in graph["nodes"]})
    spacing = xs[1]-xs[0]
    area = len(xs)*len(ys)*spacing*spacing
    capacities = np.empty(count+2)
    capacities[0] = physics["reservoir_ratio"] * area/alpha
    capacities[1] = area/alpha
    for node in graph["nodes"]:
        capacities[node["node_id"]+2] = node["V_m3"] * math.exp(-beta*a_ref/node["r_m"])
    edge_rows = []
    for edge in graph["edges"]:
        i = 0 if edge["i"] == "L" else 1 if edge["i"] == "R" else edge["i"]+2
        j = 0 if edge["j"] == "L" else 1 if edge["j"] == "R" else edge["j"]+2
        w = beta*a_ref/edge["a_m"]
        internal_resistance = edge["L_m"]*math.exp(w)/(d*edge["A_m2"])
        access_resistance = 1/(4*d*edge["a_m"]) if edge["port_flag"] and access else 0
        layer_resistance = (physics.get("external_layer_m", 0)/(d*edge["A_m2"])
                            if edge["port_flag"] else 0)
        g = 1/(internal_resistance+access_resistance+layer_resistance)
        edge_rows.append((i,j,g))
    return capacities, edge_rows, area


def generator(capacities: np.ndarray, edges: list[tuple[int,int,float]]):
    n = len(capacities)
    row, col, value = [], [], []
    for i,j,g in edges:
        row.extend([i,i,j,j])
        col.extend([i,j,i,j])
        value.extend([-g/capacities[i],g/capacities[i],g/capacities[j],-g/capacities[j]])
    return coo_matrix((value,(row,col)),shape=(n,n)).tocsc()


def closed_baths(capacities, edges, times_s: np.ndarray, initial: np.ndarray | None = None,
                 rtol: float = 1e-7, atol: float = 1e-11):
    operator = generator(capacities,edges)
    initial = np.array([1.0,0.0]+[0.0]*(len(capacities)-2)) if initial is None else initial
    result = solve_ivp(lambda _t,u: operator@u, (0,float(times_s[-1])), initial, t_eval=times_s,
                       method="BDF", jac=operator, rtol=rtol, atol=atol)
    if not result.success:
        raise RuntimeError(result.message)
    salt = capacities@result.y
    return {"u_R": result.y[1], "mass_relative_error": float(np.max(np.abs(salt-salt[0]))/salt[0]),
            "minimum_u": float(np.min(result.y)), "all_u": result.y}


def fixed_bath_p(capacities, edges, area: float) -> float:
    n = len(capacities)
    # Only membrane nodes in a component connected to both baths carry steady flow.
    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for i,j,_g in edges:
        a,b = find(i),find(j)
        if a != b:
            parent[b] = a
    if find(0) != find(1):
        return 0.0
    indices = [i for i in range(2,n) if find(i) == find(0)]
    if not indices:
        return sum(g for i,j,g in edges if {i,j} == {0,1})/area
    remap = {node:k for k,node in enumerate(indices)}
    row,col,val = [],[],[]
    rhs = np.zeros(len(indices))
    for i,j,g in edges:
        if i in remap:
            k = remap[i]
            row.append(k); col.append(k); val.append(g)
            if j in remap:
                row.append(k); col.append(remap[j]); val.append(-g)
            elif j == 0:
                rhs[k] += g
        if j in remap:
            k = remap[j]
            row.append(k); col.append(k); val.append(g)
            if i in remap:
                row.append(k); col.append(remap[i]); val.append(-g)
            elif i == 0:
                rhs[k] += g
    operator = coo_matrix((val,(row,col)),shape=(len(indices),len(indices))).tocsc()
    potentials = spsolve(operator,rhs)
    flux = sum(g*(1-potentials[remap[j]]) if i==0 and j in remap else
               g*(1-potentials[remap[i]]) if j==0 and i in remap else
               g if {i,j} == {0,1} else 0 for i,j,g in edges)
    return flux/area


def main() -> None:
    spec = json.loads(SPEC_PATH.read_text(encoding="utf-8-sig"))
    config_hash = hashlib.sha256(SPEC_PATH.read_bytes()).hexdigest()
    print(f"config sha256={config_hash}; probe H/200 nm/seed7001; alpha=1000 m^-1; beta=0,3")
    # Independent small-graph limits of the generic matrix and steady solver.
    exchange = closed_baths(np.array([1.0,1.0]),[(0,1,2.0)],np.array([0.0,0.5,1.0]))
    expected = 0.5*(1-np.exp(-4*np.array([0.0,0.5,1.0])))
    exchange_error = float(np.max(np.abs(exchange["u_R"]-expected)))
    series_p = fixed_bath_p(np.ones(3),[(0,2,2.0),(2,1,3.0)],1.0)
    parallel_p = fixed_bath_p(np.ones(2),[(0,1,2.0),(0,1,3.0)],1.0)
    blind_p = fixed_bath_p(np.ones(3),[(0,2,2.0)],1.0)
    blind = closed_baths(np.ones(3),[(0,2,2.0)],np.array([0.0,0.5,1.0]))
    disconnected_p = fixed_bath_p(np.ones(4),[(0,2,2.0),(1,3,3.0)],1.0)
    disconnected = closed_baths(np.ones(4),[(0,2,2.0),(1,3,3.0)],np.array([0.0,0.5,1.0]))
    equal = closed_baths(np.ones(3),[(0,2,2.0),(2,1,3.0)],np.array([0.0,0.5,1.0]),
                         initial=np.full(3,0.6))
    assert exchange_error < 1e-3
    assert abs(series_p-1.2) < 1e-12 and abs(parallel_p-5.0) < 1e-12
    assert blind_p == 0 and max(abs(blind["u_R"])) < 1e-12
    assert disconnected_p == 0 and max(abs(disconnected["u_R"])) < 1e-12
    assert float(np.max(np.abs(equal["all_u"]-0.6))) < 1e-12
    graph = load_graph("H",200,7001)
    rows=[]
    for beta in (0,3):
        capacities, edges, area = system(graph,beta,1e3,spec)
        full_window_minutes = np.array([0,0.001,0.01,0.1,1,5,10,20,40,60,90,120])
        result = closed_baths(capacities,edges,full_window_minutes*60.0)
        p = fixed_bath_p(capacities,edges,area)
        assert p >= 0
        assert result["mass_relative_error"] <= 1e-8
        initial_salt = float(capacities[0])
        negative_salt_fraction = float(max(0.0,-np.min(capacities[:,None]*result["all_u"]))/initial_salt)
        assert negative_salt_fraction <= 1e-12
        rows.append({"beta":beta,"P_m_s":p,"times_min":full_window_minutes.tolist(),
                     "u_R":result["u_R"].tolist(),
                     "max_mass_relative_error":result["mass_relative_error"],
                     "minimum_u":result["minimum_u"],
                     "max_negative_salt_fraction":negative_salt_fraction})
    evidence={"config_sha256":config_hash,"T2_two_node_exchange_max_abs_error":exchange_error,
              "T2_three_node_series_P":series_p,"T2_two_parallel_P":parallel_p,
              "T3_blind_fixed_P":blind_p,"T3_blind_right_max_u":float(np.max(np.abs(blind["u_R"]))),
              "T3_disconnected_fixed_P":disconnected_p,
              "T3_disconnected_right_max_u":float(np.max(np.abs(disconnected["u_R"]))),
              "T3_equal_activity_max_abs_change":float(np.max(np.abs(equal["all_u"]-0.6))),
              "network_probe":{"series":"H","h_nm":200,"seed":7001,"alpha_m_inv":1000,"times_min":[0,20,60,120],"models":rows}}
    (ROOT/"WP2/v7_network_transport_probe.json").write_text(json.dumps(evidence,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(f"T2/T3 focused limits PASS; exchange max abs error {exchange_error:.3g}; network beta0/3 mass errors "
          f"{rows[0]['max_mass_relative_error']:.3g}/{rows[1]['max_mass_relative_error']:.3g}")


if __name__ == "__main__":
    main()
