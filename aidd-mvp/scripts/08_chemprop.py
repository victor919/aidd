"""08 · Advanced：Chemprop（D-MPNN）对照

输入：data/processed/dataset_with_scaffold.csv、split_indices.json
输出：models/chemprop.pt、data/processed/08_chemprop.json

铁律：Advanced model 必须有 baseline 对照——
    同一批分子、同一个 scaffold split、同一套 metric。
"""
import os
import json
import numpy as np
import pandas as pd
import torch
from rdkit import Chem, RDLogger

from _common import (ROOT, cfg, P, ensure_dirs, processed, path,
                     log, banner, save_json, exists)

RDLogger.DisableLog("rdApp.*")
ensure_dirs()
banner("08 · Advanced：Chemprop D-MPNN（与 07 同一 scaffold split 对照）")

from chemprop.data import MoleculeDatapoint, MoleculeDataset, build_dataloader
from chemprop.featurizers import SimpleMoleculeMolGraphFeaturizer
from chemprop.nn import BondMessagePassing, MeanAggregation, RegressionFFN
from chemprop.models import MPNN
import lightning.pytorch as pl

df = pd.read_csv(processed("dataset_with_scaffold.csv"))
y_all = df["pIC50_median"].values
splits = json.load(open(processed("split_indices.json")))
tr = np.array(splits["scaffold"]["train"])
te = np.array(splits["scaffold"]["test"])
log(f"沿用 07 的 scaffold split：train={len(tr)} test={len(te)}")

m = cfg["model"]
d_h = int(m["chemprop_hidden"])
depth = int(m["chemprop_depth"])


def make_dataset(idx, ys_for_scale=None):
    mols, ys = [], []
    for i in idx:
        mol = Chem.MolFromSmiles(df.iloc[i]["smiles"])
        if mol is None:
            continue
        mols.append(mol); ys.append([float(y_all[i])])
    dps = [MoleculeDatapoint(mol, y=np.array(y)) for mol, y in zip(mols, ys)]
    return MoleculeDataset(dps, featurizer=SimpleMoleculeMolGraphFeaturizer())


train_dset = make_dataset(tr)
test_dset = make_dataset(te)

# chemprop 的 ScaleTransform 做的是「正向」标准化 forward: (X - mean) / scale，
# 不是很多人以为的「反变换」。所以正确用法是：把 y 标准化后再训练，评估时手动乘回scale。
# 若误把它挂在 FFN.output_transform 上，预测阶段会被二次标准化 —— 这就是
# 「GNN 惨败」最常见的假象来源（本轮实跑过：MAE 从 0.6 直接飙到 12）。
USE_Y_NORM = os.environ.get("AIDD_CP_YNORM", "1") == "1"
if USE_Y_NORM:
    scaler_std = train_dset.normalize_targets()
    test_dset.normalize_targets(scaler_std)
    _mean = float(np.asarray(scaler_std.mean_).ravel()[0])
    _scale = float(np.asarray(scaler_std.scale_).ravel()[0])
    output_transform = None
else:
    _mean, _scale = 0.0, 1.0
    output_transform = None

log(f"分子图构建完成：train={len(train_dset)} test={len(test_dset)}")
log(f"target 标准化: {USE_Y_NORM}  mean={_mean:.3f} scale={_scale:.3f}")

train_loader = build_dataloader(train_dset, batch_size=int(m["chemprop_batch_size"]),
                                shuffle=True, num_workers=0)
test_loader = build_dataloader(test_dset, batch_size=256, shuffle=False, num_workers=0)

torch.manual_seed(int(cfg["random_seed"]))
np.random.seed(int(cfg["random_seed"]))

mp = BondMessagePassing(d_h=d_h, depth=depth, dropout=0.05)
agg = MeanAggregation()
ffn = RegressionFFN(n_tasks=1, input_dim=d_h, hidden_dim=d_h,
                    n_layers=int(m["chemprop_ffn_layers"]), dropout=0.05)
mpnn = MPNN(message_passing=mp, agg=agg, predictor=ffn,
            warmup_epochs=2, init_lr=1e-4, max_lr=1e-3, final_lr=1e-4)

n_params = sum(p.numel() for p in mpnn.parameters() if p.requires_grad)
log(f"MPNN 可训练参数 {n_params:,}")

trainer = pl.Trainer(max_epochs=int(m["chemprop_epochs"]),
                     logger=False, enable_checkpointing=False,
                     enable_progress_bar=False, accelerator="auto",
                     devices=1, num_sanity_val_steps=0)
log(f"开始训练 {m['chemprop_epochs']} 个 epoch（M1 Pro / MPS 或 CPU）")
trainer.fit(mpnn, train_loader)

# ---------- 评估 ----------
preds = trainer.predict(mpnn, test_loader)
pred_norm = np.concatenate([np.asarray(p).reshape(-1) for p in preds])
true_norm = np.array([d.y[0] for d in test_dset])

# 若做了 y 标准化，这里统一乘回原始尺度，保证和 baseline 在同一量纲上比较
pred = pred_norm * _scale + _mean
true = true_norm * _scale + _mean

log(f"sanity check: true[mean/min/max]={true.mean():.3f}/{true.min():.3f}/{true.max():.3f}"
    f"  pred[mean/min/max]={pred.mean():.3f}/{pred.min():.3f}/{pred.max():.3f}")
if abs(pred.mean() - true.mean()) > 1.0:
    log("预测均值与真值偏差 > 1.0：多半是 target 标准化没还原，先修这个再谈模型", "WARN")

mae = float(np.abs(pred - true).mean())
rmse = float(np.sqrt(((pred - true) ** 2).mean()))
ss_res = ((true - pred) ** 2).sum(); ss_tot = ((true - true.mean()) ** 2).sum()
r2 = float(1 - ss_res / ss_tot)

# ---------- 与 baseline 对照 ----------
base = json.load(open(processed("07_baseline.json")))["scaffold"]
print(f"\n{'模型':<28}{'MAE':>9}{'RMSE':>9}{'R2':>9}")
print("-" * 55)
print(f"{'ECFP4 + RF (baseline)':<28}{base['MAE']:>9.4f}{base['RMSE']:>9.4f}{base['R2']:>9.4f}")
print(f"{'Chemprop D-MPNN':<28}{mae:>9.4f}{rmse:>9.4f}{r2:>9.4f}")
gap = mae - base["MAE"]
verdict = "Chemprop 更好" if gap < 0 else "Chemprop 这次没赢 baseline"
print(f"\n  ΔMAE = {gap:+.4f} -> {verdict}")
print("  重点不是谁赢，而是：把表示从「指纹」换成「分子图」，模型到底多拿到了什么信息？")

torch.save(mpnn.state_dict(), path(P["models"], "chemprop.pt"))
save_json({"target_normalized": bool(USE_Y_NORM), "target_mean": round(_mean,4), "target_scale": round(_scale,4),
           "n_train": int(len(train_dset)), "n_test": int(len(test_dset)),
           "MAE": round(mae, 4), "RMSE": round(rmse, 4), "R2": round(r2, 4),
           "n_params": int(n_params), "epochs": int(m["chemprop_epochs"]),
           "baseline_MAE": base["MAE"], "delta_MAE": round(gap, 4),
           "verdict": verdict}, processed("08_chemprop.json"))
banner(f"08 done · Chemprop MAE = {mae:.4f}（baseline {base['MAE']:.4f}）")
