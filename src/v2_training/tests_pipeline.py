import unittest, tempfile, json, csv, math
from pathlib import Path
import pipeline as p

class FrozenPipeline(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.out=Path(self.tmp.name)
        self.result=p.demo(self.out)
    def tearDown(self): self.tmp.cleanup()
    def test_known_synthetic_scale_and_whole_runs(self):
        self.assertAlmostEqual(self.result['fit']['permeability_m_s']/3e-7,1.,places=3)
        self.assertEqual(self.result['validation']['whole_run_n'],4)
        self.assertFalse(self.result['validation']['experimental_mechanism_validated'])
    def test_training_mother_cannot_be_validation(self):
        with self.assertRaisesRegex(ValueError,'mother IDs'):
            p.evaluate(self.out/'frozen.json',self.out/'train.csv',self.out/'bad.csv')
    def test_code_change_after_freeze_rejected(self):
        q=json.loads((self.out/'frozen.json').read_text())
        q['model_source_sha256']['models.py']='changed'
        p.write_json(self.out/'tampered.json',q)
        with self.assertRaisesRegex(ValueError,'source changed'):
            p.evaluate(self.out/'tampered.json',self.out/'validation.csv',self.out/'bad.csv')
    def test_synthetic_to_real_cannot_validate_material(self):
        text=(self.out/'validation.csv').read_text().replace('synthetic,','real,')
        (self.out/'real.csv').write_text(text)
        with self.assertRaisesRegex(ValueError,'synthetic fit'):
            p.evaluate(self.out/'frozen.json',self.out/'real.csv',self.out/'bad.csv')
    def test_constant_volume_model_rejects_sampling_trajectory(self):
        with (self.out/'validation.csv').open() as f: rows=list(csv.DictReader(f))
        rows[2]['volume_right_m3']='9e-6'
        with (self.out/'sampling.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=p.FIELDS);w.writeheader();w.writerows(rows)
        with self.assertRaisesRegex(ValueError,'Constant-volume'):
            p.read_data(self.out/'sampling.csv')
    def test_unequal_volume_prediction_conserves_salt(self):
        row={'volume_left_m3':7e-6,'volume_right_m3':11e-6,
             'c_left0_mol_m3':30.,'c_right0_mol_m3':2.,'area_m2':2e-4,'time_s':3000}
        cl,cr=p.predict(row,1e-6)
        self.assertAlmostEqual(cl*7e-6+cr*11e-6,30*7e-6+2*11e-6,places=15)
        self.assertGreaterEqual(cr,2.);self.assertLessEqual(cl,30.)
if __name__=='__main__':unittest.main(verbosity=2)
