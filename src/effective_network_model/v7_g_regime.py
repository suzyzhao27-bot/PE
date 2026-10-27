"""Frozen G NS static-limit regime map: 2 series x 12 h x 9 beta x 8 seeds."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

from v7_f_s01 import reference_parameters
from v7_network_geometry import generate, layer_count, morphology
from v7_network_transport import ROOT, SPEC_PATH, system


AUDIT = ROOT / "WP5/v7_G_geometry_audit.csv"
ROWS = ROOT / "data/processed/figures/v7_G_static_regime.csv"
SUMMARY = ROOT / "WP5/v7_G_summary.json"
FIG = ROOT / "data/processed/figures"


def steady_metrics(graph, edges, area, spec):
    """P and the fraction of steady dissipation in both port access resistances."""
    count = len(graph["nodes"]) + 2
    parent = list(range(count))

    def find(k):
        while parent[k] != k:
            parent[k] = parent[parent[k]]
            k = parent[k]
        return k

    for i, j, _ in edges:
        a, b = find(i), find(j)
        if a != b:
            parent[b] = a
    if find(0) != find(1):
        return 0.0, None, 0.0
    indices = [i for i in range(2, count) if find(i) == find(0)]
    remap = {node: k for k, node in enumerate(indices)}
    row, col, value = [], [], []
    rhs = np.zeros(len(indices))
    for i, j, g in edges:
        if i in remap:
            k = remap[i]
            row.append(k); col.append(k); value.append(g)
            if j in remap:
                row.append(k); col.append(remap[j]); value.append(-g)
            elif j == 0:
                rhs[k] += g
        if j in remap:
            k = remap[j]
            row.append(k); col.append(k); value.append(g)
            if i in remap:
                row.append(k); col.append(remap[i]); value.append(-g)
            elif i == 0:
                rhs[k] += g
    matrix = coo_matrix((value, (row, col)), shape=(len(indices), len(indices))).tocsc()
    solved = spsolve(matrix, rhs)
    potentials = {0: 1.0, 1: 0.0}
    potentials.update((i, float(solved[k])) for i, k in remap.items())
    d_na, d_cl = spec["physics"]["D_Na_m2_s"], spec["physics"]["D_Cl_m2_s"]
    d = 2*d_na*d_cl/(d_na+d_cl)
    flux = total_power = access_power = 0.0
    for edge, (i, j, g) in zip(graph["edges"], edges):
        if i not in potentials or j not in potentials:
            continue
        drop = potentials[i]-potentials[j]
        current = g*drop
        total_power += current*drop
        if i == 0:
            flux += current
        elif j == 0:
            flux -= current
        if edge["port_flag"]:
            access_power += current*current/(4*d*edge["a_m"])
    if flux <= 0:
        raise ArithmeticError("Connected graph has nonpositive flux")
    energy_error = abs(total_power-flux)/flux
    if energy_error > 1e-8:
        raise ArithmeticError(f"Steady energy mismatch: {energy_error}")
    fraction = access_power/total_power
    if not 0 <= fraction <= 1+1e-10:
        raise ArithmeticError(f"Invalid access dissipation fraction: {fraction}")
    return flux/area, fraction, energy_error


def compute():
    if AUDIT.exists() or ROWS.exists():
        raise RuntimeError("G computation already recorded")
    spec = json.loads(SPEC_PATH.read_text(encoding="utf-8-sig"))
    frozen = spec["matrices"]["G"]
    assert frozen == {"series": ["H", "S"], "h_logspace_nm": [65, 2600, 12],
                      "beta_linspace": [0, 6, 9], "seeds": 8,
                      "steady_combinations": 1728, "fit_each_point": False,
                      "extra_dynamic_2D_sweep": False}
    geom = spec["geometry"]
    hs = np.geomspace(65, 2600, 12)
    betas = np.linspace(0, 6, 9)
    max_nz = max(layer_count(float(h), geom["spacing_factor"] *
                             sum(morphology(float(h), series, geom["anchors"],
                                            geom["S_extrapolation"])[:2]),
                             geom["nz_rule"])
                 for series in frozen["series"] for h in hs)
    alpha = 10 ** reference_parameters()[0]
    audits, rows = [], []
    for series in frozen["series"]:
        for h_index, h in enumerate(hs):
            for seed in range(7001, 7009):
                graph, audit = generate(float(h), series, seed, spec, max_nz)
                audits.append(audit)
                for beta_index, beta in enumerate(betas):
                    _, edges, area = system(graph, float(beta), alpha, spec)
                    p, access_fraction, energy_error = steady_metrics(graph, edges, area, spec)
                    rows.append({"series": series, "h_index": h_index,
                                 "h_nm": float(h), "beta_index": beta_index,
                                 "beta": float(beta), "seed": seed,
                                 "P_m_s": p,
                                 "spanning_node_fraction": audit["spanning_nodes"]/audit["node_count"],
                                 "spanning_nodes": audit["spanning_nodes"],
                                 "one_side_nodes": audit["one_side_nodes"],
                                 "isolated_nodes": audit["isolated_nodes"],
                                 "access_dissipation_fraction": access_fraction,
                                 "steady_energy_relative_error": energy_error})
            print(f"G {series} h {h_index+1}/12 complete", flush=True)
    assert len(audits) == 192 and len(rows) == 1728
    for path, items in ((AUDIT, audits), (ROWS, rows)):
        with path.open("w", encoding="utf-8-sig", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(items[0]))
            writer.writeheader(); writer.writerows(items)
    print(f"G computed {len(rows)} static combinations; max volume error "
          f"{max(a['volume_relative_error'] for a in audits):.3g}; "
          f"max energy error {max(r['steady_energy_relative_error'] for r in rows):.3g}")


def plot():
    with ROWS.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 1728
    os.environ.setdefault("MPLCONFIGDIR", str(FIG))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x = np.geomspace(65, 2600, 12)
    y = np.linspace(0, 6, 9)
    all_positive = [float(r["P_m_s"]) for r in rows if float(r["P_m_s"]) > 0]
    limits = (math.floor(math.log10(min(all_positive))),
              math.ceil(math.log10(max(all_positive))))
    summary = {"model": "NS_static_limit", "combination_count": len(rows),
               "geometry_count": 192, "P_positive_range_m_s": [min(all_positive), max(all_positive)],
               "log10_P_color_limits": limits, "series": {}}
    for series in ("H", "S"):
        part = [r for r in rows if r["series"] == series]
        p = np.empty((9, 12)); connected = np.empty((9, 12)); access = np.empty((9, 12))
        for yi in range(9):
            for xi in range(12):
                cell = [r for r in part if int(r["h_index"]) == xi and int(r["beta_index"]) == yi]
                assert len(cell) == 8
                p[yi, xi] = np.mean([float(r["P_m_s"]) for r in cell])
                connected[yi, xi] = np.mean([int(r["spanning_nodes"]) > 0 for r in cell])
                access_values = [float(r["access_dissipation_fraction"]) for r in cell
                                 if r["access_dissipation_fraction"]]
                access[yi, xi] = np.mean(access_values) if access_values else np.nan
        fig, ax = plt.subplots(figsize=(8.4, 5.2), constrained_layout=True)
        logp = np.where(p > 0, np.log10(p), np.nan)
        mesh = ax.pcolormesh(x, y, logp, shading="nearest", cmap="viridis",
                             vmin=limits[0], vmax=limits[1])
        fig.colorbar(mesh, ax=ax, label="log10 P (m/s)")
        valid = access[np.isfinite(access)]
        levels = [v for v in (0.01, 0.05, 0.1, 0.25, 0.5, 0.75)
                  if len(valid) and min(valid) < v < max(valid)]
        if levels:
            contours = ax.contour(x, y, access, levels=levels, colors="white", linewidths=0.9)
            ax.clabel(contours, inline=True, fontsize=8, fmt="access %.2g")
        incomplete = np.where(connected < 1)
        if len(incomplete[0]):
            ax.scatter(x[incomplete[1]], y[incomplete[0]], marker="x", color="red",
                       s=35, label="<8/8 spanning realizations")
            ax.legend(loc="lower left", fontsize=8)
        for anchor in (230, 370, 770, 1860):
            ax.axvline(anchor, color="white", alpha=0.35, linewidth=0.7, linestyle="--")
        ax.set(xscale="log", xlim=(65, 2600), ylim=(0, 6), xlabel="Thickness h (nm)",
               ylabel="Barrier beta", title=f"{series}: NS static-limit regime (8 seeds/cell)")
        connectivity_note = ("all 8/8 span" if not len(incomplete[0])
                             else "red x: fewer than 8/8 span")
        ax.text(0.01, 0.99,
                f"Dashed: structure anchors; contours: access dissipation; {connectivity_note}",
                transform=ax.transAxes, va="top", fontsize=8, color="white",
                bbox={"facecolor": "black", "alpha": 0.35, "edgecolor": "none"})
        stem = FIG / f"v7_G_{series}_static_regime"
        fig.savefig(stem.with_suffix(".svg"))
        fig.savefig(stem.with_suffix(".png"), dpi=180)
        plt.close(fig)
        summary["series"][series] = {
            "mean_P_range_m_s": [float(np.min(p)), float(np.max(p))],
            "spanning_realization_fraction_range": [float(np.min(connected)), float(np.max(connected))],
            "access_dissipation_fraction_range": [float(np.min(valid)), float(np.max(valid))]
            if len(valid) else None,
            "access_contour_levels": levels}
    SUMMARY.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("G plotted H/S regime maps")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("compute", "plot"))
    args = parser.parse_args()
    compute() if args.action == "compute" else plot()


if __name__ == "__main__":
    main()
