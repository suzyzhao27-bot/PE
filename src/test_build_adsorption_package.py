import csv
import importlib.util
import math
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "src" / "build_adsorption_package.py"
spec = importlib.util.spec_from_file_location("build_adsorption_package", MODULE_PATH)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


class AdsorptionPackageContractTests(unittest.TestCase):
    def test_baseline_rows_reproduce_declared_ratios_and_free_energy_scales(self):
        rows = module.build_baseline_rows()
        by_name = {row["comparison"]: row for row in rows}
        self.assertAlmostEqual(float(by_name["thin_over_thick"]["ratio"]), 1.837, places=3)
        self.assertAlmostEqual(float(by_name["middle_over_thick"]["ratio"]), 1.179, places=3)
        self.assertAlmostEqual(float(by_name["thin_over_middle"]["ratio"]), 1.559, places=3)
        self.assertAlmostEqual(float(by_name["thin_over_thick"]["delta_delta_g_kj_mol"]), -1.51, places=2)
        self.assertAlmostEqual(float(by_name["noncrystalline_bound"]["ratio"]), 1.29, places=2)

    def test_structure_table_has_three_processing_groups_and_source_labels(self):
        rows = module.structure_rows()
        self.assertEqual([row["draw"] for row in rows], ["182", "357", "604"])
        self.assertEqual([row["second_draw_ratio"] for row in rows], [1.5, 2.1, 2.73])
        self.assertTrue(all(row["evidence_role"] == "task_brief_constraint" for row in rows))

    def test_package_manifest_is_self_consistent_after_build(self):
        with tempfile.TemporaryDirectory() as temp:
            package = module.build_package(Path(temp))
            manifest = package / "manifest.csv"
            self.assertTrue(manifest.is_file())
            with manifest.open(encoding="utf-8-sig", newline="") as stream:
                rows = list(csv.DictReader(stream))
            self.assertGreaterEqual(len(rows), 12)
            listed = (package / "03_文件清单.txt").read_text(encoding="utf-8").splitlines()
            self.assertTrue(all(row["relative_path"] in listed for row in rows))
            for row in rows:
                self.assertTrue((package / row["relative_path"]).is_file(), row["relative_path"])
            self.assertTrue((package / "00_导师先读.txt").is_file())
            self.assertTrue((package / "01_科研推进与复算记录.md").is_file())
            self.assertTrue((package / "02_计算边界与输入缺口.txt").is_file())
            self.assertTrue((package / "03_文件清单.txt").is_file())
            self.assertTrue((package / "figures/图1_结构约束与模拟边界.png").is_file())
            self.assertTrue((package / "figures/图2_量级标尺与交付状态.png").is_file())
            self.assertTrue((package / "figures/图3_几何外层代理.png").is_file())
            self.assertTrue((package / "figure_data/20261109_几何外层代理_自由膜.csv").is_file())
            self.assertTrue((package / "figure_data/20261109_几何代理合同.json").is_file())
            self.assertEqual(package.name, "32_拉伸结构与小分子吸附_导师侧计算审阅包_20261109_v1.0.1")


if __name__ == "__main__":
    unittest.main()
