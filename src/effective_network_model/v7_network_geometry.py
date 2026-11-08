"""Frozen v7 H/S proxy-network geometry and T1 geometry-budget audit.

Output is geometry only. Passing T1 does not validate transport or fitting.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import ndtr, ndtri


ROOT = Path(__file__).resolve().parents[2]
SPEC_PATH = ROOT / "原始资料/赵珈绮/14_论文补强学生任务书_v7.0.0/configs/frozen_spec.json"
OUT = ROOT / "WP2"


def morphology(h_nm: float, series: str, anchors: list[dict],
               extrapolation: str = "nearest_endpoint") -> tuple[float, float, float]:
    if series == "H":
        anchor = next(a for a in anchors if a["draw"] == 357)
        return anchor["df_nm"], anchor["dp_nm"], anchor["phi_A"]
    assert series == "S"
    assert extrapolation in ("nearest_endpoint", "linear_in_log_h")
    ordered = sorted(anchors, key=lambda a: a["h_nm"])
    if h_nm <= ordered[0]["h_nm"]:
        if extrapolation == "linear_in_log_h":
            lo, hi = ordered[:2]
            w = (math.log(h_nm)-math.log(lo["h_nm"])) / (math.log(hi["h_nm"])-math.log(lo["h_nm"]))
            return tuple((1-w)*lo[key]+w*hi[key] for key in ("df_nm", "dp_nm", "phi_A"))
        a = ordered[0]
        return a["df_nm"], a["dp_nm"], a["phi_A"]
    if h_nm >= ordered[-1]["h_nm"]:
        if extrapolation == "linear_in_log_h":
            lo, hi = ordered[-2:]
            w = (math.log(h_nm)-math.log(lo["h_nm"])) / (math.log(hi["h_nm"])-math.log(lo["h_nm"]))
            return tuple((1-w)*lo[key]+w*hi[key] for key in ("df_nm", "dp_nm", "phi_A"))
        a = ordered[-1]
        return a["df_nm"], a["dp_nm"], a["phi_A"]
    for lo, hi in zip(ordered, ordered[1:]):
        if lo["h_nm"] <= h_nm <= hi["h_nm"]:
            w = (math.log(h_nm)-math.log(lo["h_nm"])) / (math.log(hi["h_nm"])-math.log(lo["h_nm"]))
            return tuple((1-w)*lo[key]+w*hi[key] for key in ("df_nm", "dp_nm", "phi_A"))
    raise AssertionError(h_nm)


def layer_count(h_nm: float, s_nm: float, rule: str) -> int:
    if rule == "max(1,ceil(h/s))":
        return max(1, math.ceil(h_nm/s_nm))
    if rule == "round_half_up":
        return max(1, math.floor(h_nm/s_nm + 0.5))
    raise ValueError(rule)


def truncated_lognormal(uniform: np.ndarray, mean: float, sigma: float, lo: float, hi: float) -> np.ndarray:
    if sigma == 0:
        return np.full_like(uniform, mean)

    def mean_at(mu: float) -> float:
        za, zb = (math.log(lo)-mu)/sigma, (math.log(hi)-mu)/sigma
        denominator = ndtr(zb)-ndtr(za)
        numerator = ndtr(zb-sigma)-ndtr(za-sigma)
        return math.exp(mu+sigma*sigma/2) * numerator/denominator

    assert lo < mean < hi
    center = math.log(mean)
    lower, upper = center-1, center+1
    while mean_at(lower) > mean:
        lower -= 1
    while mean_at(upper) < mean:
        upper += 1
    root = brentq(lambda mu: mean_at(mu)-mean, lower, upper, xtol=1e-13)
    cdf_lo = ndtr((math.log(lo)-root)/sigma)
    cdf_hi = ndtr((math.log(hi)-root)/sigma)
    return np.exp(root + sigma*ndtri(cdf_lo + uniform*(cdf_hi-cdf_lo)))


class UnionFind:
    def __init__(self, size: int):
        self.parent = list(range(size))

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra


def generate(h_nm: float, series: str, seed: int, spec: dict, max_nz: int) -> tuple[dict, dict]:
    geom = spec["geometry"]
    nx, ny = geom["nx"], geom["ny"]
    assert nx >= 3 and ny >= 3 and geom["periodic_xy"]
    df_nm, dp_nm, phi = morphology(h_nm, series, geom["anchors"],
                                  geom["S_extrapolation"])
    s = geom["spacing_factor"]*(df_nm+dp_nm)*1e-9
    h = h_nm*1e-9
    nz = layer_count(h_nm, s*1e9, geom["nz_rule"])
    dz = h/nz
    epsilon = geom["kappa_V"]*phi
    total_area = (nx*s)*(ny*s)
    radius_mean = dp_nm*0.5e-9
    lo, hi = 0.25*radius_mean, min(2*radius_mean, 0.49*s)
    radius_rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([seed, 0])))
    axial_rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([seed, 1])))
    radius_u = radius_rng.random((max_nz, ny, nx))[:nz]
    axial_u = axial_rng.random((max_nz-1, ny, nx))[:max(0, nz-1)]
    radius = truncated_lognormal(radius_u, radius_mean, geom["sigma_ln"], lo, hi)
    axial_kept = axial_u < geom["p_z"]
    node_volume = epsilon*s*s*dz
    nodes = []
    index = lambda z, y, x: (z*ny+y)*nx+x
    for z in range(nz):
        for y in range(ny):
            for x in range(nx):
                nodes.append({"node_id": index(z,y,x), "x_m": (x+0.5)*s, "y_m": (y+0.5)*s,
                              "z_m": (z+0.5)*dz, "V_m3": node_volume, "r_m": float(radius[z,y,x]),
                              "reservoir_role": "membrane"})
    edges = []
    seen = set()
    capped = 0

    def add_edge(i: int | str, j: int | str, length: float, face: float, raw_radius: float, port: bool) -> None:
        nonlocal capped
        key = (str(i), str(j))
        assert i != j and key not in seen
        seen.add(key)
        cap_radius = math.sqrt(epsilon*face/math.pi)
        a = float(min(raw_radius, cap_radius))
        was_capped = bool(raw_radius > cap_radius)
        capped += was_capped
        assert length > 0 and a > 0 and face > 0
        edges.append({"edge_id": len(edges), "i": i, "j": j, "L_m": length,
                      "a_m": a, "A_m2": math.pi*a*a, "port_flag": port,
                      "area_capped": was_capped})

    for z in range(nz):
        for y in range(ny):
            for x in range(nx):
                i = index(z,y,x)
                r = radius[z,y,x]
                for yy, xx in ((y,(x+1)%nx), ((y+1)%ny,x)):
                    j = index(z,yy,xx)
                    add_edge(i,j,s,s*dz,geom["q"]*min(r,radius[z,yy,xx]),False)
                if z+1 < nz and axial_kept[z,y,x]:
                    j = index(z+1,y,x)
                    add_edge(i,j,dz,s*s,geom["q"]*min(r,radius[z+1,y,x]),False)
                if z == 0:
                    add_edge("L",i,dz/2,s*s,geom["q"]*r,True)
                if z == nz-1:
                    add_edge(i,"R",dz/2,s*s,geom["q"]*r,True)

    union = UnionFind(len(nodes)+2)
    for edge in edges:
        i = len(nodes) if edge["i"] == "L" else len(nodes)+1 if edge["i"] == "R" else edge["i"]
        j = len(nodes) if edge["j"] == "L" else len(nodes)+1 if edge["j"] == "R" else edge["j"]
        union.union(i,j)
    left, right = union.find(len(nodes)), union.find(len(nodes)+1)
    by_component = {}
    for node in nodes:
        component = union.find(node["node_id"])
        node["component_id"] = component
        by_component[component] = by_component.get(component, 0)+1
    spanning_count = sum(v for component,v in by_component.items() if component == left == right)
    one_side_count = sum(v for component,v in by_component.items() if (component == left) != (component == right))
    isolated_count = len(nodes)-spanning_count-one_side_count
    left_area = sum(e["A_m2"] for e in edges if e["i"] == "L")
    right_area = sum(e["A_m2"] for e in edges if e["j"] == "R")
    volume_error = abs(sum(n["V_m3"] for n in nodes)-epsilon*total_area*h)/(epsilon*total_area*h)
    assert volume_error < 1e-12
    assert left_area <= epsilon*total_area*(1+1e-12)
    assert right_area <= epsilon*total_area*(1+1e-12)
    assert spanning_count+one_side_count+isolated_count == len(nodes)
    graph = {"series": series, "h_nm": h_nm, "seed": seed, "nodes": nodes, "edges": edges}
    encoded = json.dumps(graph, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")
    digest = hashlib.sha256(encoded).hexdigest()
    summary = {"series": series, "h_nm": h_nm, "seed": seed, "geometry_hash": digest,
               "df_nm": df_nm, "dp_nm": dp_nm, "phi_A": phi, "s_nm": s*1e9, "nz": nz,
               "node_count": len(nodes), "edge_count": len(edges), "sample_r_mean_nm": float(radius.mean()*1e9),
               "epsilon": epsilon, "volume_relative_error": volume_error,
               "left_port_area_fraction": left_area/(epsilon*total_area),
               "right_port_area_fraction": right_area/(epsilon*total_area),
               "area_capped_edges": capped, "spanning_nodes": spanning_count,
               "one_side_nodes": one_side_count, "isolated_nodes": isolated_count}
    return graph, summary


def main() -> None:
    spec = json.loads(SPEC_PATH.read_text(encoding="utf-8-sig"))
    geom = spec["geometry"]
    hs = spec["matrices"]["C1"]["h_nm"]
    series_list = spec["matrices"]["C1"]["series"]
    max_nz = max(math.ceil(h/((morphology(h, series, geom["anchors"])[0]+morphology(h, series, geom["anchors"])[1]))) for h in hs for series in series_list)
    config_hash = hashlib.sha256(SPEC_PATH.read_bytes()).hexdigest()
    print(f"config sha256={config_hash}; H/S h={hs} nm; nx=ny={geom['nx']}; seeds={geom['main_seeds']}; max_nz={max_nz}")
    summaries = []
    with gzip.open(OUT / "v7_network_tables.jsonl.gz", "wt", encoding="utf-8", newline="\n") as stream:
        for series in series_list:
            for h in hs:
                for seed in geom["main_seeds"]:
                    graph, summary = generate(h, series, seed, spec, max_nz)
                    stream.write(json.dumps(graph, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n")
                    summary["config_sha256"] = config_hash
                    summaries.append(summary)
    with (OUT / "v7_network_geometry_audit.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)
    assert len(summaries) == 6*2*8
    print(f"T1 geometry checks PASS: {len(summaries)} networks; max volume error {max(s['volume_relative_error'] for s in summaries):.3g}; "
          f"isolated-network nodes {sum(s['isolated_nodes'] for s in summaries)}")


if __name__ == "__main__":
    main()
