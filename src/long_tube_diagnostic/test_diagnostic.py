import hashlib
import json
import math
import time
import unittest
from pathlib import Path
from long_tube_np import local_pb_transport, local_donnan_transport, integrate_tube
from pb_v1_readonly_copy import (R, F, ambipolar_diffusivity,
                                 donnan_transport_integral)

ROOT = Path(__file__).resolve().parent
METRICS = {}


class DiagnosticTests(unittest.TestCase):
    def test_zero_sigma_ambipolar_and_henderson(self):
        dp, dm = 1.33e-9, 2.03e-9
        h = (dm - dp) / (dm + dp)
        for cl in (30, 270):
            result = integrate_tube(cl, 1., 0, cells=64, intervals=256)
            exact_flux = ambipolar_diffusivity(dp, dm) * (cl - 1) / 350e-9
            exact_delta_mv = R * 298.15 / F * h * math.log(1 / cl) * 1000
            self.assertAlmostEqual(result["pore_area_salt_flux_mol_m2_s"], exact_flux, delta=abs(exact_flux) * 2e-8)
            self.assertAlmostEqual(result["virtual_delta_phi_right_minus_left_mv"], exact_delta_mv, delta=1e-10)
        METRICS["neutral_virtual_delta_30_to_1_mv"] = exact_delta_mv if cl == 30 else R * 298.15 / F * h * math.log(1 / 30) * 1000

    def test_zero_current_local_and_pb_conservation(self):
        worst_i = worst_gauss = 0.
        for c in (.1, 1, 10, 100, 270):
            for sigma in (-.01, .01):
                value = local_pb_transport(c, sigma, 25e-9, 512)
                gp, gm = 1.33e-9 * c * value["k_plus"], 2.03e-9 * c * value["k_minus"]
                h = value["d_reduced_virtual_potential_d_log_c"]
                # Arbitrary axial dln(c)/dx=-1 m^-1; use NP independently.
                jp, jm = gp * (1 + h), gm * (1 - h)
                relative_i = abs(jp - jm) / max(jp, jm)
                self.assertLess(relative_i, 1e-12)
                self.assertAlmostEqual(jp / c, value["effective_salt_D_m2_s"], delta=jp / c * 1e-12)
                self.assertLess(value["pb_residual_scaled"], 1e-10)
                self.assertLess(value["pb_gauss_error_scaled"], 1e-8)
                worst_i = max(worst_i, relative_i)
                worst_gauss = max(worst_gauss, value["pb_gauss_error_scaled"])
        METRICS["worst_relative_zero_current_error"] = worst_i
        METRICS["worst_scaled_Gauss_error_test_points"] = worst_gauss

    def test_charge_ion_exchange_and_flow_direction(self):
        neg = integrate_tube(30., 1., -.01, cells=128, intervals=128)
        swapped = integrate_tube(30., 1., .01, cells=128, intervals=128,
                                 d_plus=2.03e-9, d_minus=1.33e-9)
        reverse = integrate_tube(1., 30., -.01, cells=128, intervals=128)
        self.assertAlmostEqual(neg["pore_area_salt_flux_mol_m2_s"], swapped["pore_area_salt_flux_mol_m2_s"], places=12)
        self.assertAlmostEqual(neg["virtual_delta_phi_right_minus_left_mv"], -swapped["virtual_delta_phi_right_minus_left_mv"], places=10)
        self.assertAlmostEqual(neg["pore_area_salt_flux_mol_m2_s"], -reverse["pore_area_salt_flux_mol_m2_s"], places=12)
        self.assertAlmostEqual(neg["virtual_delta_phi_right_minus_left_mv"], -reverse["virtual_delta_phi_right_minus_left_mv"], places=10)

    def test_donnan_uniform_overlap_limit(self):
        # Small radius and weak radial variation; not the 25nm sigma=.01 scenarios.
        radius, sigma = 1e-9, -1e-4
        max_relative = 0.
        for c in (.1, .3, 1.0):
            pb = local_pb_transport(c, sigma, radius, 256)
            don = local_donnan_transport(c, sigma, radius)
            rel = abs(pb["effective_salt_D_m2_s"] - don["effective_salt_D_m2_s"]) / don["effective_salt_D_m2_s"]
            max_relative = max(max_relative, rel)
            self.assertLess(rel, 1e-4)
        integral = integrate_tube(1., .1, sigma, radius=radius, cells=128, intervals=256, approximation="donnan")
        exact = donnan_transport_integral(1, .1, sigma, radius, 1.33e-9, 2.03e-9) / 350e-9
        self.assertAlmostEqual(integral["pore_area_salt_flux_mol_m2_s"], exact, delta=abs(exact) * 1e-8)
        METRICS["weak_radial_overlap_PB_Donnan_max_relative_D_difference"] = max_relative
        METRICS["Donnan_virtual_c_integral_vs_q_integral_relative_error"] = abs(integral["pore_area_salt_flux_mol_m2_s"] - exact) / abs(exact)

    def test_radial_grid_convergence(self):
        values = [integrate_tube(270., 1., -.01, cells=n, intervals=256) for n in (256, 512, 1024)]
        ratios = [v["charged_to_neutral_salt_flux_ratio"] for v in values]
        voltages = [v["virtual_delta_phi_right_minus_left_mv"] for v in values]
        differences = [abs(ratios[i + 1] - ratios[i]) for i in (0, 1)]
        self.assertLess(differences[1], differences[0] * .4)
        self.assertLess(differences[1], 1e-5)
        self.assertLess(abs(voltages[2] - voltages[1]), .01)
        METRICS["radial_grid_convergence"] = dict(cells=[256, 512, 1024], ratios=ratios,
                                                  virtual_delta_mv=voltages, ratio_delta_512_1024=differences[1])

    def test_axial_quadrature_convergence(self):
        values = [integrate_tube(270., 1., -.01, cells=512, intervals=n) for n in (64, 128, 256)]
        ratios = [v["charged_to_neutral_salt_flux_ratio"] for v in values]
        voltages = [v["virtual_delta_phi_right_minus_left_mv"] for v in values]
        differences = [abs(ratios[i + 1] - ratios[i]) for i in (0, 1)]
        self.assertLess(differences[1], differences[0] * .15)
        self.assertLess(differences[1], 1e-7)
        self.assertLess(abs(voltages[2] - voltages[1]), .001)
        METRICS["log_c_quadrature_convergence"] = dict(intervals=[64, 128, 256], ratios=ratios,
                                                      virtual_delta_mv=voltages, ratio_delta_128_256=differences[1])

    def test_registered_copy_hash_identity(self):
        source_identity = json.loads((ROOT / "pb_source_identity.json").read_text(encoding="utf-8"))
        actual = hashlib.sha256((ROOT / "pb_v1_readonly_copy.py").read_bytes()).hexdigest()
        self.assertEqual(actual, source_identity["registered_source_sha256"])
        METRICS["v1_PB_copy_sha256"] = actual
        METRICS["original_absolute_file_is_not_a_runtime_dependency"] = True


if __name__ == "__main__":
    started = time.perf_counter()
    tests = unittest.defaultTestLoader.loadTestsFromTestCase(DiagnosticTests)
    result = unittest.TextTestRunner(verbosity=2).run(tests)
    report = dict(status="PASS" if result.wasSuccessful() else "FAIL", tests_run=result.testsRun,
                  failures=len(result.failures), errors=len(result.errors),
                  runtime_seconds=time.perf_counter() - started, metrics=METRICS,
                  scope="restricted long-tube diagnostic tests; no finite-length/entrance PNP or PE validation")
    (ROOT / "test_report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    raise SystemExit(0 if result.wasSuccessful() else 1)
