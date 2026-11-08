import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "src" / "surface_geometry_proxy.py"
spec = importlib.util.spec_from_file_location("surface_geometry_proxy", MODULE_PATH)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


class SurfaceGeometryProxyTests(unittest.TestCase):
    def test_summary_contains_three_existing_family_geometries(self):
        rows = module.build_proxy_rows()
        self.assertEqual([row["draw"] for row in rows], ["182", "357", "604"])
        self.assertEqual([row["N"] for row in rows], [11, 6, 4])
        self.assertTrue(all(row["evidence_role"] == "legacy_geometry_output" for row in rows))

    def test_one_nm_outer_shell_is_geometry_proxy_only(self):
        rows = module.build_proxy_rows()
        shell = {row["draw"]: row["one_nm_shell_fraction"] for row in rows}
        self.assertGreater(shell["357"], shell["182"])
        self.assertGreater(shell["182"], shell["604"])
        self.assertEqual({row["model_status"] for row in rows}, {"GEOMETRY_PROXY_NOT_ATOMIC_SURFACE"})

    def test_proxy_does_not_emit_adsorption_quantities(self):
        rows = module.build_proxy_rows()
        self.assertTrue(all("W_z" not in row and "K_dry" not in row and "R_dry" not in row for row in rows))


if __name__ == "__main__":
    unittest.main()
