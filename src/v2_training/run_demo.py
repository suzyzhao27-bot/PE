#!/usr/bin/env python3
"""Generate labelled synthetic examples; no fitted or fabricated experiments."""
import argparse
import csv
import json
import math
import platform
import time
from dataclasses import asdict, replace
from pathlib import Path
from xml.sax.saxutils import escape
from models import (Geometry, ambipolar_diffusivity, cylindrical_pb, debye_length,
                    donnan_partition, neutral_curve, neutral_resistances,
                    reduced_curve, reduced_donnan_flux)


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def line_svg(path, title, xlabel, ylabel, series, xlimits, ylimits, footer):
    """Small standard-library SVG export; no plot dependencies needed."""
    w, h = 920, 580
    left, top, pw, ph = 95, 70, 600, 390
    xmin, xmax = xlimits
    ymin, ymax = ylimits
    colors = ["#1261a0", "#d75442", "#16877d", "#9962a7", "#d18b17", "#666666"]
    def xp(x):
        return left + (x - xmin) / (xmax - xmin) * pw
    def yp(y):
        return top + ph - (y - ymin) / (ymax - ymin) * ph
    items = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
             '<rect width="100%" height="100%" fill="white"/>',
             '<g font-family="Arial, sans-serif" fill="#222">',
             f'<text x="{left}" y="28" font-size="19">{escape(title)}</text>',
             f'<text x="{left}" y="48" font-size="12">SYNTHETIC SCENARIOS — NO EXPERIMENTAL FIT</text>']
    for i in range(6):
        x = xmin + (xmax - xmin) * i / 5
        y = ymin + (ymax - ymin) * i / 5
        items.extend([f'<line x1="{xp(x):.2f}" y1="{top}" x2="{xp(x):.2f}" y2="{top+ph}" stroke="#e2e2e2"/>',
                      f'<line x1="{left}" y1="{yp(y):.2f}" x2="{left+pw}" y2="{yp(y):.2f}" stroke="#e2e2e2"/>',
                      f'<text x="{xp(x):.2f}" y="{top+ph+23}" font-size="12" text-anchor="middle">{x:.3g}</text>',
                      f'<text x="{left-12}" y="{yp(y)+4:.2f}" font-size="12" text-anchor="end">{y:.3g}</text>'])
    items.extend([f'<rect x="{left}" y="{top}" width="{pw}" height="{ph}" fill="none" stroke="#333"/>',
                  f'<text x="{left+pw/2}" y="{top+ph+48}" font-size="14" text-anchor="middle">{escape(xlabel)}</text>',
                  f'<text transform="translate(24 {top+ph/2}) rotate(-90)" font-size="14" text-anchor="middle">{escape(ylabel)}</text>'])
    for i, (label, points) in enumerate(series):
        color = colors[i % len(colors)]
        coords = " ".join(f"{xp(x):.2f},{yp(y):.2f}" for x, y in points)
        items.append(f'<polyline points="{coords}" fill="none" stroke="{color}" stroke-width="2"/>')
        items.extend([f'<line x1="715" y1="{90+28*i}" x2="740" y2="{90+28*i}" stroke="{color}" stroke-width="2"/>',
                      f'<text x="746" y="{94+28*i}" font-size="12">{escape(label)}</text>'])
    for i, text in enumerate(footer):
        items.append(f'<text x="{left}" y="{535+i*17}" font-size="11">{escape(text)}</text>')
    items.append('</g></svg>')
    path.write_text("\n".join(items), encoding="utf-8")


def legacy_scale_audit():
    """Derived from AM2022 reported J; assumptions are NOT a fitted geometry."""
    j, delta_c, length = 3.23e-6, 270.0, 350e-9
    p_target = j / delta_c
    d_app = p_target * length
    d, radius, phi = 1.6e-9, 25e-9, .5
    g = Geometry(radius_m=radius, length_m=length, porosity=phi,
                 layer_left_m=0, layer_right_m=0)
    p_neutral = neutral_resistances(g, d)["permeability_m_s"]
    fw_required = p_target / p_neutral
    delta_total = d / p_target - (length + math.pi * radius / 2) / phi
    return dict(source="Sun et al. Advanced Materials 2022, DOI 10.1002/adma.202107941; legacy 04 audit",
                scope="derived scale check only; no digitized curve, no fit, no measured fw/delta",
                reported_salt_flux_mol_m2_s=j, stated_concentration_difference_mol_m3=delta_c,
                stated_length_m=length, derived_permeability_m_s=p_target,
                derived_apparent_D_m2_s=d_app, hypothetical_radius_m=radius,
                hypothetical_porosity=phi, hypothetical_free_D_m2_s=d,
                hypothetical_fw=1.0, assumed_external_layers_m=0.0,
                hypothetical_neutral_permeability_m_s=p_neutral,
                neutral_to_derived_permeability_ratio=p_neutral / p_target,
                hypothetical_fw_for_same_P_without_external_layers=fw_required,
                hypothetical_total_planar_external_layer_m_for_same_P_at_fw1=delta_total,
                caveat="Independent circular pores, known projected area, no convection; open-aperture controls must constrain outer resistance. This scale check cannot identify the mechanism.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path(__file__).with_name("demo_output"))
    parser.add_argument("--cells", type=int, default=512)
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    concentrations = [.1, 1.0, 10.0, 100.0, 270.0, 600.0]  # mol/m3 = mM numerically
    sigmas = [0.0, -.001, -.01, -.05]
    debye = [dict(c_bulk_mol_m3=c, c_bulk_mol_l=c / 1000,
                  debye_nm=debye_length(c) * 1e9, radius_nm=25.0,
                  kappa_a=25e-9 / debye_length(c), scenario_status="synthetic ideal 1:1 salt")
             for c in concentrations]
    write_csv(out / "debye_scale.csv", debye)
    sweep, profiles = [], []
    selected_profiles = []
    for c in concentrations:
        for sigma in sigmas:
            pb = cylindrical_pb(c, sigma, cells=args.cells)
            don = donnan_partition(c, sigma, 25e-9)
            sweep.append(dict(c_bulk_mol_m3=c, sigma_c_m2=sigma,
                              pb_coion_partition=pb["coion_partition"],
                              donnan_overlap_coion_partition=don["coion_partition"],
                              pb_center_mv=pb["center_potential_mv"],
                              pb_wall_extrapolated_mv=pb["wall_potential_extrapolated_mv"],
                              debye_nm=pb["debye_nm"], kappa_a=pb["kappa_a"],
                              a_over_gouy_chapman=pb["a_over_gouy_chapman"],
                              residual_scaled=pb["residual_scaled"], gauss_error_scaled=pb["gauss_error_scaled"],
                              cells=args.cells, model_level="radial equilibrium PB; NOT axial flux"))
            for row in pb["profile"]:
                profiles.append(dict(c_bulk_mol_m3=c, sigma_c_m2=sigma, **row))
            if sigma == -.01 and c in (.1, 1, 10, 270, 600):
                selected_profiles.append((f"{c:g} mM", [(row["r_over_a"], row["potential_mv"]) for row in pb["profile"]]))
    write_csv(out / "pb_vs_donnan.csv", sweep)
    write_csv(out / "pb_profiles.csv", profiles)
    line_svg(out / "pb_profiles.svg", "Radial nonlinear PB: a = 25 nm, sigma = -0.01 C/m2",
             "r/a (center 0; wall 1)", "Potential relative to virtual bulk (mV)", selected_profiles,
             (0, 1), (-155, 0),
             ["Level 2: equilibrium cross-section. Ideal fixed charge; epsilon_r = 78.3, T = 298.15 K.",
              "Wall potential is not zeta/KPFM. These profiles do not calculate membrane salt flux."])

    g = Geometry(layer_left_m=10e-6, layer_right_m=10e-6)
    dp, dm = 1.33e-9, 2.03e-9
    d = ambipolar_diffusivity(dp, dm)
    fluxes = []
    for c in concentrations:
        for fw in (1.0, .01, .0001):
            for sigma in (0.0, -.01, -.05):
                gg = replace(g, wetted_fraction=fw)
                neutral = neutral_resistances(gg, d)["permeability_m_s"] * c
                reduced = reduced_donnan_flux(gg, c, 0, sigma, dp, dm)
                fluxes.append(dict(c_left_mol_m3=c, c_right_mol_m3=0,
                                    fw_assumed=fw, sigma_assumed_c_m2=sigma,
                                    neutral_flux_mol_m2_s=neutral,
                                    reduced_donnan_flux_mol_m2_s=reduced,
                                    reduced_to_neutral_flux_ratio=reduced / neutral,
                                    external_layer_each_m=10e-6,
                                    applicability="neutral exact ideal geometry" if sigma == 0 else "Donnan homogeneous overlap approximation; high salt/nonuniform potential invalid for quantitative PE prediction",
                                    status="synthetic reduced model; NOT finite-length PNP or fit"))
    write_csv(out / "salt_wetting_charge_flux.csv", fluxes)
    outer_sweep = []
    for outer_each in (0, 1e-6, 10e-6, 100e-6):
        for fw in (1.0, .01, .0001):
            gg = replace(g, layer_left_m=outer_each, layer_right_m=outer_each, wetted_fraction=fw)
            outer_sweep.append(dict(external_layer_each_m=outer_each, fw_assumed=fw,
                                    **neutral_resistances(gg, d), status="synthetic geometry"))
    write_csv(out / "external_wetting_resistances.csv", outer_sweep)

    times = [float(t) for t in range(0, 7201, 120)]
    curves, series = [], []
    for label, fw, sigma, reduced in (
        ("neutral fw=1", 1., 0., False),
        ("reduced -0.01", 1., -.01, True),
        ("reduced -0.05", 1., -.05, True),
        ("neutral fw=0.01", .01, 0., False),
        ("neutral fw=0.0001", .0001, 0., False)):
        gg = replace(g, wetted_fraction=fw)
        for c0 in (1., 270.):
            rows = reduced_curve(gg, c0, 0, times, sigma, dp, dm, max_step_s=5) if reduced else neutral_curve(gg, c0, 0, times, d)
            for row in rows:
                curves.append(dict(scenario=label, c_left0_mol_m3=c0, sigma_c_m2=sigma, fw=fw,
                                   model="Donnan EN zero-current reduced" if reduced else "neutral exact", **row))
            if c0 == 270:
                series.append((label, [(row["time_s"] / 60, row["c_right_mol_m3"]) for row in rows]))
    write_csv(out / "finite_reservoir_curves.csv", curves)
    line_svg(out / "finite_reservoir_curves.svg", "Finite reservoirs: 270 mM to DI, 10 mL on each side",
             "Time (min)", "Receiving salt concentration (mol/m3 = mM)", series,
             (0, 120), (0, 140),
             ["Level 0 neutral exact / Level 1 Donnan electroneutral reduced; area = 1 cm2; a = 25 nm; L = 350 nm.",
              "phi = 0.4; outer layer = 10 um per side for all curves. Unknown parameters are assumptions, not data."])
    legacy = legacy_scale_audit()
    (out / "legacy_scale_audit.json").write_text(json.dumps(legacy, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    write_csv(out / "legacy_scale_audit.csv", [legacy])
    elapsed = time.perf_counter() - started
    manifest = dict(status="synthetic runnable demo; numerical checks separate from experiment validation",
                    python=platform.python_version(), runtime_seconds=elapsed, pb_cells=args.cells,
                    pb_scenarios=len(sweep), flux_scenarios=len(fluxes), finite_reservoir_rows=len(curves),
                    default_geometry=asdict(g), d_plus_m2_s=dp, d_minus_m2_s=dm, d_ambipolar_m2_s=d,
                    provenance="legacy 04_界面电动与盐输运核查.txt; all curves synthetic; source PDFs unchanged",
                    warning="Neither PB nor the reduced Donnan model is a finite-length PNP implementation; no charge origin inferred.")
    (out / "run_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(dict(output=str(out), runtime_seconds=round(elapsed, 3),
                          pb_scenarios=len(sweep), flux_scenarios=len(fluxes),
                          legacy_fw_scale=legacy["hypothetical_fw_for_same_P_without_external_layers"]), ensure_ascii=False))


if __name__ == "__main__":
    main()
