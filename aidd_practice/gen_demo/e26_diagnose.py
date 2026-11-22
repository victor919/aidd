#!/usr/bin/env python
"""e26 · 双模型分歧诊断（有标签测试集可验证；生成候选只报分布）

- 测试集：用 ⑫ 已保存的 scaffold seed42 RF/GNN npz（1770 分子，有真值），检查
  分歧与误差的关系、一致但双错案例 —— 可验证。
- 生成候选：重建 ⑫ 的 GNN seed42 模型（state_dict），对 e24/e25 生成的候选预测；
  只报分歧分布，不声称"筛掉了错的"（候选无实验标签）。
输出 results/e26_disagreement.json
"""
import os, json, time, random
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from rdkit import Chem, RDLogger
from gen_common import World

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
RES = os.path.join(HERE, "results")
GN = os.path.join(ROOT, "aidd_practice", "gnn_demo", "results")
DATA = os.path.join(ROOT, "aidd-mvp", "data", "processed")
T0 = time.time()
def log(*a): print(f"[{time.time()-T0:6.1f}s]", *a, flush=True)

out = {}

# ---------- 1. 有标签测试集：分歧能否预警误差 ----------
r = np.load(os.path.join(GN, "pred_scaffold_rf_seed42.npz"))
g = np.load(os.path.join(GN, "pred_scaffold_gnn_seed42.npz"))
y, pr, pg = r["true"], r["pred"], g["pred"]
er, eg = np.abs(pr - y), np.abs(pg - y)
d = np.abs(pr - pg)
out["test_set"] = dict(
    n=int(len(y)),
    RF_MAE=round(float(er.mean()), 4), GNN_MAE=round(float(eg.mean()), 4),
    gnn_better_frac=round(float((eg < er).mean()), 4),
    spearman_disagree_vs_rf_err=round(float(spearmanr(d, er).statistic), 4),
    spearman_disagree_vs_max_err=round(float(spearmanr(d, np.maximum(er, eg)).statistic), 4),
    spearman_disagree_vs_min_err=round(float(spearmanr(d, np.minimum(er, eg)).statistic), 4),
)
# 分箱：分歧越小，最大误差是否越小？（可验证的分箱统计）
qs = np.quantile(d, [.25, .5, .75])
bins = [("Q1 一致", d <= qs[0]), ("Q2", (d > qs[0]) & (d <= qs[1])),
        ("Q3", (d > qs[1]) & (d <= qs[2])), ("Q4 分歧", d > qs[2])]
out["test_set"]["bins"] = [dict(bin=name, n=int(m.sum()),
                                max_err_med=round(float(np.median(np.maximum(er, eg)[m])), 4),
                                frac_max_err_gt1=round(float((np.maximum(er, eg)[m] > 1).mean()), 4))
                           for name, m in bins]
# 一致但双错案例
mask = (d < 0.5) & (er > 2.0) & (eg > 2.0)
ii = np.where(mask)[0]
out["test_set"]["n_agree_but_both_wrong"] = int(mask.sum())
out["test_set"]["agree_but_both_wrong_top5"] = [
    dict(true=round(float(y[i]), 3), RF=round(float(pr[i]), 3), GNN=round(float(pg[i]), 3))
    for i in ii[np.argsort(-y[ii])][:5]]
log("test set done", out["test_set"]["n"], "n_agree_but_both_wrong", int(mask.sum()))

# ---------- 2. 生成候选：重建 GNN 模型 ----------
import torch
from chemprop.models import MPNN
from chemprop.nn import BondMessagePassing, MeanAggregation, RegressionFFN
from chemprop.data import MoleculeDatapoint, MoleculeDataset, build_dataloader
from chemprop.featurizers import SimpleMoleculeMolGraphFeaturizer
import lightning.pytorch as pl

sd = torch.load(os.path.join(GN, "model_scaffold_gnn_seed42.pt"), map_location="cpu", weights_only=False)
model = MPNN(BondMessagePassing(d_h=300, depth=3, dropout=0.05), MeanAggregation(),
             RegressionFFN(n_tasks=1, input_dim=300, hidden_dim=300, n_layers=2, dropout=0.05),
             warmup_epochs=2, init_lr=1e-4, max_lr=1e-3, final_lr=1e-4)
model.load_state_dict(sd); model.eval()
df = pd.read_csv(os.path.join(DATA, "dataset_with_scaffold.csv"))
yy = df["pIC50_median"].values
sp = json.load(open(os.path.join(DATA, "split_indices.json")))
tr_idx = np.array(sp["scaffold"]["train"])
mean, scale = float(yy[tr_idx].mean()), float(yy[tr_idx].std())

w = World()
def gnn_predict(smis):
    f = SimpleMoleculeMolGraphFeaturizer()
    dps = [MoleculeDatapoint(Chem.MolFromSmiles(s), y=np.array([0.])) for s in smis]
    dl = build_dataloader(MoleculeDataset(dps, featurizer=f), batch_size=256,
                          shuffle=False, num_workers=0)
    tr = pl.Trainer(logger=False, enable_checkpointing=False, enable_progress_bar=False,
                    accelerator="cpu", devices=1, num_sanity_val_steps=0)
    pn = np.concatenate([np.asarray(z).reshape(-1) for z in tr.predict(model, dl)])
    return pn * scale + mean

# 收集生成候选（e25 的保留记录 + e24 的候选池摘要）
rec = json.load(open(os.path.join(RES, "e25_closedloop.json")))["all_records"]  # 全量 27 次运行（原 audit_records 只存最后一次，已废弃）
smis = list(dict.fromkeys([x["smiles"] for x in rec]))
log("unique generated candidates", len(smis))
if len(smis) > 1200:
    smis = random.Random(42).sample(smis, 1200)   # 控制推理耗时，固定抽样
predg = gnn_predict(smis)
w.oracle_fit()
rows = []
for s, pgv in zip(smis, predg):
    p = w.pool.get(s) or w.evaluate(s)
    rows.append(dict(smiles=s, RF=p["pred"], GNN=round(float(pgv), 4),
                     disagree=round(abs(p["pred"] - float(pgv)), 4), nn=p["nn"],
                     SA=p["SA"], QED=p["QED"], MW=p["MW"]))
dis = np.array([x["disagree"] for x in rows])
nn_ = np.array([x["nn"] for x in rows])
out["generated"] = dict(
    n=len(rows), n_sampled=len(smis),
    disagree_median=round(float(np.median(dis)), 4),
    disagree_q90=round(float(np.quantile(dis, .9)), 4),
    frac_disagree_gt1=round(float((dis > 1).mean()), 4),
    frac_nn_lt_05=round(float((nn_ < .5).mean()), 4),
    spearman_disagree_vs_nn=round(float(spearmanr(dis, nn_).statistic), 4),
    note="生成候选无实验标签：这里只报告分歧分布，不声称筛掉错误候选",
    top_disagree=sorted(rows, key=lambda x: -x["disagree"])[:5])
# 固定规则选例（四条规则各选一个；同一结构可占多个位置，不事后替换）
import pandas as pd
_df = pd.read_csv(os.path.join(HERE, "..", "..", "aidd-mvp", "data", "processed", "dataset_with_scaffold.csv"))
_fit = set(_df.smiles.values[np.array(json.load(open(os.path.join(RES, "e23_indices.json")))["fit"])])
_rules = {
    "RF 预测最高": max(rows, key=lambda x: x["RF"]),
    "最近邻最低": min(rows, key=lambda x: x["nn"]),
    "分歧最大": max(rows, key=lambda x: x["disagree"]),
    "SA 最高": max(rows, key=lambda x: x["SA"]),
}
out["generated"]["cards"] = [
    dict(rule=k, smiles=v["smiles"], RF=v["RF"], GNN=v["GNN"], disagree=v["disagree"],
         nn=v["nn"], SA=v["SA"], QED=v["QED"], MW=v["MW"],
         in_fit=bool(v["smiles"] in _fit))
    for k, v in _rules.items()]
json.dump(out, open(os.path.join(RES, "e26_disagreement.json"), "w", encoding="utf-8"),
          ensure_ascii=False, indent=2)
log("wrote e26_disagreement.json")
