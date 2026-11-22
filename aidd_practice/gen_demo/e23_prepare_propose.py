#!/usr/bin/env python
"""e23 · ⑬ 阶段1：准备拟合池、oracle、BRICS 片段池与"可回查"的局部提议函数

输出（固定版本）：
  results/e23_pool.json          片段池与数据集口径
  results/e23_fragments.json     片段列表（供 e24/e25 复用）
  results/e23_oracle_meta.json   oracle 协议与拟合池索引摘要

纪律：
  - oracle 与最近邻参照库都只用拟合池 7148；测试集 1770 完全不进生成循环。
  - 提议函数是"局部片段替换"：同一函数供 A/B/C 三臂共用，保证差异只来自保留规则。
"""
import os, json, random, time
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger, DataStructs
from rdkit.Chem import BRICS, Descriptors, QED
from rdkit import RDConfig
import sys
sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PROC = os.path.join(ROOT, "aidd-mvp", "data", "processed")
RES = os.path.join(HERE, "results")
os.makedirs(RES, exist_ok=True)
T0 = time.time()
def log(*a): print(f"[{time.time()-T0:6.1f}s]", *a, flush=True)
def dump(name, obj):
    p = os.path.join(RES, name)
    json.dump(obj, open(p, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    log("wrote", name)

random.seed(42); np.random.seed(42)

df = pd.read_csv(os.path.join(PROC, "dataset_with_scaffold.csv"))
X = np.load(os.path.join(PROC, "X_ecfp.npy"))
y = df["pIC50_median"].values.astype(float)
sp = json.load(open(os.path.join(PROC, "split_indices.json")))
tr_all = np.array(sp["scaffold"]["train"]); te = np.array(sp["scaffold"]["test"])
log("scaffold train/test", len(tr_all), len(te))
# 拟合池/验证池：与 e12 同协议（7148 / 1550），seed 42
rng = np.random.RandomState(42); order = rng.permutation(len(tr_all))
n_val = 1550
val = tr_all[order[:n_val]]; fit = tr_all[order[n_val:]]
log("fit/val", len(fit), len(val))

# 片段池：只从拟合池提取（与 ⑥下 的全量 3310 口径不同，正文按本口径引用）
pieces = set()
n_ok = 0
for i in fit:
    m = Chem.MolFromSmiles(df.iloc[i]["smiles"])
    if m is None: continue
    try:
        for f in BRICS.BRICSDecompose(m):
            if isinstance(f, str) and "*" in f:
                pieces.add(f)
        n_ok += 1
    except Exception:
        pass
pieces = sorted(pieces)
log("BRICS 片段池", len(pieces), "来自", n_ok, "个拟合池分子")

# 按连接点数分组，供局部替换时保持"接口数"一致
by_dummy = {}
for p in pieces:
    k = p.count("*")
    by_dummy.setdefault(k, []).append(p)
log("片段连接点分布", {k: len(v) for k, v in sorted(by_dummy.items())})

# 起点：拟合池中 pIC50 高/中/低各一个（同一批起点供三臂复用）
ordr = np.argsort(y[fit])
starts = [int(fit[ordr[0]]), int(fit[ordr[len(fit)//2]]), int(fit[ordr[-1]])]
log("起点 pIC50", [round(float(y[i]), 3) for i in starts])

dump("e23_pool.json", dict(
    n=len(df), n_fit=len(fit), n_val=len(val), n_test=len(te),
    n_pieces=len(pieces), n_decomposed=n_ok,
    piece_dummy_counts={str(k): len(v) for k, v in sorted(by_dummy.items())},
    start_indices=starts, start_pic50=[round(float(y[i]), 4) for i in starts],
    oracle="ECFP4(2048,r=2)+RF(500 trees, seed42)，拟合池 7148；最近邻参照库同为拟合池",
    note="测试集 1770 不进生成循环，仅在 e26 分歧诊断中使用已保存的 npz 预测"))
dump("e23_fragments.json", dict(pieces=pieces, by_dummy={str(k): v for k, v in by_dummy.items()}))
json.dump(dict(fit=fit.tolist(), val=val.tolist(), test=te.tolist()),
          open(os.path.join(RES, "e23_indices.json"), "w"), indent=1)
log("done")
