"""Meaningful numerical/physical checks, runnable with Python stdlib."""
import json
import math
import unittest
from dataclasses import replace
from pathlib import Path
from models import (F, R, EPS0, Geometry, ambipolar_diffusivity, cylindrical_pb,
                    debye_length, donnan_partition, donnan_transport_integral,
                    mg_l_to_si, neutral_curve, neutral_resistances,
                    reduced_curve, reduced_donnan_flux, si_to_mg_l)

METRICS = {}


def bessel_i0(x):
    term = total = 1.0
    for k in range(1, 200):
        term *= x * x / (4 * k * k)
        total += term
        if abs(term) < 1e-16 * total:
            break
    return total


def bessel_i1(x):
    term = total = x / 2
    for k in range(1, 200):
        term *= x * x / (4 * k * (k + 1))
        total += term
        if abs(term) < 1e-16 * abs(total):
            break
    return total


class ModelTests(unittest.TestCase):
    def test_units_debye(self):
        self.assertAlmostEqual(si_to_mg_l(1000), 58442.77)
        self.assertAlmostEqual(mg_l_to_si(58442.77), 1000)
        lam = debye_length(1000) * 1e9
        self.assertTrue(0.302 < lam < 0.306)
        self.assertAlmostEqual(debye_length(270) / debye_length(600), math.sqrt(600 / 270))
        METRICS["debye_at_1M_nm"] = lam

    def test_neutral_resistance_limits(self):
        d = 1.6e-9
        g = replace(Geometry(), layer_left_m=0, layer_right_m=0)
        r = neutral_resistances(g, d)
        density = g.porosity / (math.pi * g.radius_m ** 2)
        independent_per_pore = g.length_m / (d * math.pi * g.radius_m ** 2) + 1 / (2 * d * g.radius_m)
        self.assertAlmostEqual(r["total_s_m"], independent_per_pore / density)
        for length, expected_large in ((g.radius_m / 1000, False), (g.radius_m * 1000, True)):
            rr = neutral_resistances(replace(g, length_m=length), d)
            self.assertEqual(rr["internal_s_m"] > rr["access_s_m"], expected_large)
        r_half = neutral_resistances(replace(g, wetted_fraction=.5), d)
        self.assertAlmostEqual(r_half["permeability_m_s"], r["permeability_m_s"] / 2)
        self.assertEqual(neutral_resistances(replace(g, wetted_fraction=0), d)["permeability_m_s"], 0)
        # At fixed porosity larger circular pores increase access resistance per area.
        # At fixed individual pore density that conclusion would be different.
        r_two = neutral_resistances(replace(g, radius_m=2 * g.radius_m), d)
        self.assertAlmostEqual(r_two["access_s_m"], 2 * r["access_s_m"])

    def test_finite_reservoir_analytical_balance(self):
        g = replace(Geometry(), volume_left_m3=8e-6, volume_right_m3=13e-6)
        d = ambipolar_diffusivity()
        rows = neutral_curve(g, 270, 2, [0, 60, 1200, 7200, 1e9], d)
        mass = 270 * g.volume_left_m3 + 2 * g.volume_right_m3
        for row in rows:
            self.assertAlmostEqual(row["total_salt_mol"], mass, delta=mass * 1e-13)
        self.assertAlmostEqual(rows[-1]["c_left_mol_m3"], mass / (g.volume_left_m3 + g.volume_right_m3))
        self.assertAlmostEqual(rows[-1]["c_right_mol_m3"], rows[-1]["c_left_mol_m3"])
        p = neutral_resistances(g, d)["permeability_m_s"]
        initial_derivative = g.area_m2 * p * 268 / g.volume_right_m3
        tiny = neutral_curve(g, 270, 2, [1e-4], d)[0]
        self.assertAlmostEqual((tiny["c_right_mol_m3"] - 2) / 1e-4, initial_derivative,
                               delta=initial_derivative * 1e-6)
        METRICS["finite_reservoir_salt_mass_mol"] = mass

    def test_donnan_charge_and_sign(self):
        a = 25e-9
        negative = donnan_partition(1, -.01, a)
        positive = donnan_partition(1, .01, a)
        self.assertAlmostEqual(negative["c_plus_mol_m3"] * negative["c_minus_mol_m3"], 1)
        self.assertAlmostEqual(negative["c_plus_mol_m3"] - negative["c_minus_mol_m3"],
                               -negative["fixed_charge_signed_mol_m3"])
        self.assertAlmostEqual(negative["c_plus_mol_m3"], positive["c_minus_mol_m3"])
        self.assertAlmostEqual(negative["potential_v"], -positive["potential_v"])
        self.assertEqual(donnan_partition(0, -.01, a)["c_coion_mol_m3"], 0)

    def test_pb_neutral_and_nonlinear_gauss(self):
        zero = cylindrical_pb(1, 0, cells=64)
        self.assertEqual(zero["coion_partition"], 1)
        self.assertEqual(zero["wall_potential_extrapolated_mv"], 0)
        worst_residual = worst_gauss = 0
        for c in (.1, 1, 10, 270, 600):
            pb = cylindrical_pb(c, -.01, cells=384)
            self.assertLess(pb["residual_scaled"], 1e-10)
            self.assertLess(pb["gauss_error_scaled"], 1e-9)
            net = pb["mobile_charge_mean_c_m3"] + pb["wall_charge_per_pore_volume_c_m3"]
            self.assertLess(abs(net) / abs(pb["wall_charge_per_pore_volume_c_m3"]), 1e-8)
            self.assertTrue(0 < pb["coion_partition"] <= 1)
            worst_residual = max(worst_residual, pb["residual_scaled"])
            worst_gauss = max(worst_gauss, pb["gauss_error_scaled"])
        METRICS["pb_worst_scaled_residual"] = worst_residual
        METRICS["pb_worst_scaled_gauss_error"] = worst_gauss

    def test_pb_sign_exchange(self):
        minus = cylindrical_pb(1, -.02, cells=128)
        plus = cylindrical_pb(1, .02, cells=128)
        self.assertAlmostEqual(minus["c_plus_partition"], plus["c_minus_partition"], places=10)
        self.assertAlmostEqual(minus["wall_potential_extrapolated_mv"], -plus["wall_potential_extrapolated_mv"], places=10)
        self.assertAlmostEqual(minus["coion_partition"], plus["coion_partition"], places=11)

    def test_pb_grid_convergence(self):
        rows = []
        for c, sigma in ((1, -.01), (600, -.05)):
            values = [cylindrical_pb(c, sigma, cells=n) for n in (128, 256, 512)]
            differences = [abs(values[i + 1]["coion_partition"] - values[i]["coion_partition"]) for i in (0, 1)]
            self.assertLess(differences[1], differences[0] * .4)
            self.assertLess(differences[1], 5e-5)
            walls = [v["wall_potential_extrapolated_mv"] for v in values]
            self.assertLess(abs(walls[2] - walls[1]), abs(walls[1] - walls[0]) * .6)
            rows.append(dict(c_mol_m3=c, sigma_c_m2=sigma, cells=[128, 256, 512],
                             coion_partition=[v["coion_partition"] for v in values],
                             delta_256_512=differences[1], wall_mv=walls))
        METRICS["pb_grid_convergence"] = rows

    def test_pb_linear_debye_huckel_benchmark(self):
        c, sigma, a, n = 1.0, -1e-7, 25e-9, 256
        pb = cylindrical_pb(c, sigma, a, n)
        ka = a / debye_length(c)
        b = sigma * a / (EPS0 * 78.3 * R * 298.15 / F)
        factor = b / (ka * bessel_i1(ka))
        vt = R * 298.15 / F
        errors = [abs(row["potential_mv"] / (1000 * vt) - factor * bessel_i0(ka * row["r_over_a"]))
                  for row in pb["profile"]]
        relative = max(errors) / abs(factor * bessel_i0(ka))
        self.assertLess(relative, 2e-5)
        METRICS["pb_linear_bessel_max_relative_error"] = relative

    def test_reduced_flux_neutral_sign_and_quadrature(self):
        g = replace(Geometry(), layer_left_m=30e-6, layer_right_m=90e-6)
        dp, dm = 1.33e-9, 2.03e-9
        p = neutral_resistances(g, ambipolar_diffusivity(dp, dm))["permeability_m_s"]
        for cl, cr in ((270, 0), (0, 270), (1, .1), (.1, 1)):
            self.assertAlmostEqual(reduced_donnan_flux(g, cl, cr, 0, dp, dm), p * (cl - cr), delta=abs(p * (cl - cr)) * 1e-11)
        swapped = replace(g, layer_left_m=g.layer_right_m, layer_right_m=g.layer_left_m)
        self.assertAlmostEqual(reduced_donnan_flux(g, 10, 1, -.01),
                               -reduced_donnan_flux(swapped, 1, 10, -.01), places=12)
        self.assertAlmostEqual(reduced_donnan_flux(g, 10, 1, -.01, dp, dm),
                               reduced_donnan_flux(g, 10, 1, .01, dm, dp), places=12)
        # Independently integrate D_EN(q) by Simpson rule to test analytical primitive.
        charge = 2 * .01 / (F * g.radius_m)
        ql = donnan_partition(10, -.01, g.radius_m)["c_coion_mol_m3"]
        qr = donnan_partition(1, -.01, g.radius_m)["c_coion_mol_m3"]
        n, h = 2048, (ql - qr) / 2048
        def integrand(q):
            return dp * dm * (2 * q + charge) / (dp * (q + charge) + dm * q)
        weighted = integrand(qr) + integrand(ql)
        weighted += sum((4 if i % 2 else 2) * integrand(qr + i * h) for i in range(1, n))
        simpson = weighted * h / 3
        exact = donnan_transport_integral(10, 1, -.01, g.radius_m, dp, dm)
        self.assertAlmostEqual(simpson, exact, delta=abs(exact) * 1e-10)
        METRICS["donnan_transport_primitive_relative_error"] = abs(simpson - exact) / abs(exact)

    def test_reduced_transient_neutral_and_timestep(self):
        g = replace(Geometry(), volume_left_m3=8e-6, volume_right_m3=13e-6)
        times = [0, 1200, 2400, 3600, 7200]
        exact = neutral_curve(g, 270, 0, times, ambipolar_diffusivity())
        numeric = reduced_curve(g, 270, 0, times, 0, max_step_s=30)
        error = max(abs(a["c_right_mol_m3"] - b["c_right_mol_m3"]) for a, b in zip(exact, numeric))
        self.assertLess(error, 1e-7)
        a = reduced_curve(g, 1, 0, times, -.01, max_step_s=60)
        b = reduced_curve(g, 1, 0, times, -.01, max_step_s=30)
        delta = max(abs(x["c_right_mol_m3"] - y["c_right_mol_m3"]) for x, y in zip(a, b))
        self.assertLess(delta, 1e-7)
        mass = g.volume_left_m3
        for row in b:
            self.assertAlmostEqual(row["total_salt_mol"], mass, delta=mass * 1e-12)
        METRICS["neutral_RK4_vs_exact_max_error_mol_m3"] = error
        METRICS["charged_RK4_60_30s_max_difference_mol_m3"] = delta


if __name__ == "__main__":
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(ModelTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    report = dict(status="PASS" if result.wasSuccessful() else "FAIL", tests_run=result.testsRun,
                  failures=len(result.failures), errors=len(result.errors), metrics=METRICS,
                  scope="numerical model checks; not experiment validation")
    destination = Path(__file__).with_name("test_report.json")
    destination.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    raise SystemExit(0 if result.wasSuccessful() else 1)
