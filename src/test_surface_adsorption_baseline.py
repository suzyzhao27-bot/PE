import csv
import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "src" / "mechanism_surface_adsorption_baseline.py"
spec = importlib.util.spec_from_file_location("mechanism_surface_adsorption_baseline", MODULE_PATH)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


class SurfaceAdsorptionBaselineTests(unittest.TestCase):
    def test_structure_constraints_are_three_group_contract(self):
        rows = module.structure_rows()
        self.assertEqual([row["draw"] for row in rows], ["182", "357", "604"])
        self.assertEqual([row["thickness_nm"] for row in rows], [770.0, 370.0, 230.0])
        self.assertEqual([row["second_draw_ratio"] for row in rows], [1.5, 2.1, 2.73])

    def test_ratio_baseline_matches_declared_charge_reference(self):
        summary = module.compute_baseline()
        self.assertAlmostEqual(summary["thin_over_thick"], 1.8373448688786695, places=12)
        self.assertAlmostEqual(summary["middle_over_thick"], 1.1786830627378397, places=12)
        self.assertAlmostEqual(summary["thin_over_middle"], 1.558811632204923, places=12)
        self.assertAlmostEqual(summary["noncrystalline_bound"], 1.2860247772223432, places=12)

    def test_audit_output_preserves_unknown_simulation_status(self):
        audit = module.build_audit_rows()
        statuses = {row["item"]: row["status"] for row in audit}
        self.assertEqual(statuses["dry_neutral_octanoic_acid_w_of_z"], "BLOCKED_MISSING_MD_INPUTS")
        self.assertEqual(statuses["wet_HA_and_A_minus"], "NOT_STARTED_BY_DESIGN")
        self.assertEqual(statuses["zeta_site_density_reference"], "REPRODUCED_FROM_SUPPLIED_CSV")


if __name__ == "__main__":
    unittest.main()
