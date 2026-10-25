"""只读检查现有产物；不把历史产物检查等同于重新训练。"""
import json
import numpy as np
import pandas as pd
from _common import processed


def check():
    df = pd.read_csv(processed('dataset_with_scaffold.csv'))
    x = np.load(processed('X_ecfp.npy'), mmap_mode='r')
    assert len(df) == len(x), '特征和标签行数不一致'
    assert df.smiles.is_unique and df.pIC50_median.notna().all()
    splits = json.load(open(processed('split_indices.json')))
    for name, split in splits.items():
        tr, te = split['train'], split['test']
        assert not set(tr) & set(te)
        assert set(tr) | set(te) == set(range(len(df)))
        if name == 'scaffold':
            assert not set(df.iloc[tr].scaffold.fillna('')) & set(df.iloc[te].scaffold.fillna(''))
        pred = pd.read_csv(processed(f'predictions_{name}.csv'))
        assert len(pred) == len(te)
        assert pred.smiles.tolist() == df.iloc[te].smiles.tolist()
        assert np.allclose(pred.true_pIC50, df.iloc[te].pIC50_median)
    dock = pd.read_csv(processed('docking.csv'))
    assert len(dock) > 0 and np.isfinite(dock[['docking_score', 'consensus']]).all().all()
    assert dock.consensus.is_monotonic_decreasing
    hist = pd.read_csv(processed('iteration_history.csv'))
    for _, group in hist.groupby('strategy'):
        assert group.n_labeled.is_monotonic_increasing
        assert np.isfinite(group.MAE).all()
    print(f'PASS: {len(df)} molecules; split/prediction/docking/history integrity')


if __name__ == '__main__':
    check()
