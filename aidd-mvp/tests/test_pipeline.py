"""在临时目录真实执行清洗、切分和 RF，避免覆盖历史全量结果。"""
import json
import os
import copy
import runpy
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
import joblib

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import _common as common


class PipelineTest(unittest.TestCase):
    def test_curation_split_training(self):
        old_cfg = copy.deepcopy(common.cfg)
        old_paths = common.P.copy()
        old_trees = common.cfg['model']['rf_n_estimators']
        try:
            with tempfile.TemporaryDirectory() as temp:
                for key in common.P:
                    common.P[key] = str(Path(temp) / key)
                common.cfg['model']['rf_n_estimators'] = 8
                common.ensure_dirs()
                source = pd.read_csv(ROOT / 'data/raw/raw_activity.csv').head(400)
                source.to_csv(common.raw('raw_activity.csv'), index=False)
                for stage in ('05_standardize.py', '06_scaffold_split.py', '07_rf_baseline.py'):
                    runpy.run_path(str(ROOT / 'scripts' / stage), run_name='__main__')
                df = pd.read_csv(common.processed('dataset_with_scaffold.csv'))
                x = np.load(common.processed('X_ecfp.npy'))
                splits = json.load(open(common.processed('split_indices.json')))
                self.assertEqual(len(x), len(df))
                self.assertTrue(df.smiles.is_unique)
                tr, te = splits['scaffold']['train'], splits['scaffold']['test']
                self.assertFalse(set(df.iloc[tr].scaffold.fillna('')) & set(df.iloc[te].scaffold.fillna('')))
                model = joblib.load(Path(common.P['models']) / 'rf_scaffold.joblib')
                predictions = pd.read_csv(common.processed('predictions_scaffold.csv'))
                np.testing.assert_allclose(model.predict(x[te]), predictions.pred_pIC50)
                self.assertTrue(np.isfinite(predictions.pred_pIC50).all())
                if os.environ.get('AIDD_TEST_FULL') == '1':
                    common.cfg['model']['chemprop_epochs'] = 1
                    common.cfg['docking'].update(n_screen=2, exhaustiveness=1, n_poses=2, cpu=1)
                    common.cfg['closed_loop'].update(n_rounds=2, batch_size=5, initial_labeled=20)
                    common.P['raw'] = str(ROOT / 'data/raw')
                    for stage in ('01_target_rank.py', '02_download_structure.py',
                                  '03_extract_pocket.py', '08_chemprop.py', '09_docking.py',
                                  '10_admet.py', '11_closed_loop.py', 'make_report.py', 'check_artifacts.py'):
                        runpy.run_path(str(ROOT / 'scripts' / stage), run_name='__main__')
                    docking = json.load(open(common.processed('09_docking.json')))
                    self.assertNotIn('error', docking['redock_check'])
                    self.assertEqual(docking['n_docked'], 2)
                    self.assertTrue((Path(common.P['models']) / 'chemprop.pt').exists())
        finally:
            common.P.update(old_paths)
            common.cfg.clear()
            common.cfg.update(old_cfg)
            common.P = common.cfg['paths']


if __name__ == '__main__':
    unittest.main()
