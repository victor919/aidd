# 文章关联：AI制药实践篇B_实操篇_公众号文章.html
# 脚本职责：早期实验 4：性质面板、Pareto 与主动学习
# 使用注解：依赖 step3 的数据表及指纹；性质计算与离线主动学习不等于实验验证。
# 完整对应关系与运行顺序：README.md
"""EGFR MVP Step4: ADMET panel -> Pareto front -> active learning loop（本机实跑）"""
import os, json, time, warnings
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen, QED, rdMolDescriptors
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, r2_score

warnings.filterwarnings("ignore")
RDLogger.DisableLog("rdApp.*")

BASE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(BASE, "data", "egfr_mvp")
agg = pd.read_csv(os.path.join(DATA, "dataset_with_scaffold.csv"))
agg["scaffold"] = agg["scaffold"].fillna("")   # 无环分子 Murcko 骨架为空串
X = np.load(os.path.join(DATA, "X_ecfp.npy"))
y = agg["pIC50_median"].values
print(f"loaded: {len(agg)} molecules, X {X.shape}")

print("\n" + "=" * 72)
print("[0] 数据体检：先看看有没有物理上不可能的数字")
print("=" * 72)
print("pIC50 quantiles:", np.round(np.percentile(y, [0, 0.1, 1, 25, 50, 75, 99, 99.9, 100]), 2).tolist())
print("对应的 IC50 (nM):",
      [f"{10**(9-v):.3g}" for v in np.percentile(y, [0, 0.1, 1, 25, 50, 75, 99, 99.9, 100])])
bad_hi = y > 11.0     # IC50 < 0.01 nM
bad_lo = y < 3.0      # IC50 > 1 mM
print(f"\npIC50 > 11 (IC50 < 0.01 nM) : {bad_hi.sum()} 个分子")
print(f"pIC50 < 3  (IC50 > 1 mM)    : {bad_lo.sum()} 个分子")
print("\n-- 最极端的 5 个 pIC50（这些多半是单位/录入错误，不是超强分子）--")
ext = pd.DataFrame({"smiles": agg["std_smiles"], "pIC50": y,
                    "IC50_nM": 10 ** (9 - y), "n_meas": agg["n_meas"]})
print(ext.sort_values("pIC50", ascending=False).head(5).to_string(index=False,
      float_format=lambda v: f"{v:.4g}"))
keep_qc = ~(bad_hi | bad_lo)
print(f"\nQC 后保留 {keep_qc.sum()} / {len(agg)} 个分子 ({keep_qc.mean()*100:.1f}%)")

print("\n" + "=" * 72)
print("[1] ADMET / developability descriptor panel")
print("=" * 72)

def panel(smiles):
    m = Chem.MolFromSmiles(smiles)
    return dict(
        MW=Descriptors.MolWt(m),
        logP=Crippen.MolLogP(m),
        TPSA=Descriptors.TPSA(m),
        HBD=rdMolDescriptors.CalcNumHBD(m),
        HBA=rdMolDescriptors.CalcNumHBA(m),
        RotB=rdMolDescriptors.CalcNumRotatableBonds(m),
        AromaticRings=rdMolDescriptors.CalcNumAromaticRings(m),
        QED=QED.qed(m),
        FormalCharge=Chem.GetFormalCharge(m),
    )

t0 = time.time()
props = pd.DataFrame([panel(s) for s in agg["std_smiles"]])
print(f"computed {len(props)} molecules in {time.time()-t0:.1f}s")
full = pd.concat([agg[["std_smiles", "pIC50_median", "n_meas", "scaffold"]], props], axis=1)

print("\n-- distribution --")
print(props[["MW", "logP", "TPSA", "HBD", "HBA", "RotB", "QED"]].describe().round(2).to_string())

print("\n-- Lipinski & Veber rule pass rates (on the whole EGFR set) --")
lip = ((full.MW <= 500) & (full.logP <= 5) & (full.HBD <= 5) & (full.HBA <= 10))
veber = ((full.RotB <= 10) & (full.TPSA <= 140))
print(f"  Lipinski  (MW<=500, logP<=5, HBD<=5, HBA<=10) : {lip.mean()*100:.1f}%  ({lip.sum()}/{len(full)})")
print(f"  Veber     (RotB<=10, TPSA<=140)               : {veber.mean()*100:.1f}%  ({veber.sum()}/{len(full)})")
print(f"  both                                          : {(lip & veber).mean()*100:.1f}%")

print("\n-- 只看最强的一批分子（pIC50 >= 8，即 IC50 <= 10 nM）--")
strong = full[full.pIC50_median >= 8]
print(f"  n = {len(strong)}")
if len(strong) > 0:
    print("  Lipinski pass: %.1f%%" % (((strong.MW<=500)&(strong.logP<=5)&(strong.HBD<=5)&(strong.HBA<=10)).mean()*100))
    print(strong[["MW","logP","TPSA","QED"]].describe().round(2).to_string())
    weak = full[full.pIC50_median < 7]
    print("\n-- 对比：pIC50 < 7 的弱活性分子 --")
    print(f"  n = {len(weak)}")
    print(weak[["MW","logP","TPSA","QED"]].describe().round(2).to_string())
    print("\n  mean MW    : strong %.1f vs weak %.1f" % (strong.MW.mean(), weak.MW.mean()))
    print("  mean logP  : strong %.2f vs weak %.2f" % (strong.logP.mean(), weak.logP.mean()))
    print("  mean QED   : strong %.3f vs weak %.3f" % (strong.QED.mean(), weak.QED.mean()))

full.to_csv(os.path.join(DATA, "dataset_with_props.csv"), index=False)

print("\n" + "=" * 72)
print("[2] 多目标：为什么不能 sort_values('score') 取 Top1")
print("=" * 72)
# 先看：如果不做 QC，Top1 会是谁？
bad = full[full.pIC50_median > 11].sort_values("pIC50_median", ascending=False)
print(f"-- 若直接按活性排序（未做 QC），Top1 是：--")
if len(bad):
    r0 = bad.iloc[0]
    print(f"  pIC50={r0.pIC50_median:.2f}  (IC50 = {10**(9-r0.pIC50_median):.3g} nM)  MW={r0.MW:.0f}")
    print(f"  SMILES: {r0.std_smiles}")
    print("  这个数字比任何真实药物的 IC50 低了约 8 个数量级 —— 它来自单位录入错误。")
# 教学用：活性 + 成药性 的双目标，找 Pareto front（只在 QC 通过的分子里做）
sub = full[(full.pIC50_median >= 7) & (full.pIC50_median <= 11)].copy()
print(f"\ncandidates with 7 <= pIC50 <= 11 (QC 通过): {len(sub)}")
# 目标1：预测/实测活性越高越好；目标2：QED 越高越好；约束：MW<=500, logP<=5
sub = sub.reset_index(drop=True)
a = sub.pIC50_median.values
b = sub.QED.values

def pareto_mask(pts):
    """pts: (n,2) 越大越好，返回帕累托前沿布尔掩码"""
    n = len(pts)
    keep = np.ones(n, dtype=bool)
    for i in range(n):
        dom = np.all(pts >= pts[i], axis=1) & np.any(pts > pts[i], axis=1)
        if dom.any():
            keep[i] = False
    return keep

pts = np.stack([a, b], axis=1)
mask = pareto_mask(pts)
print(f"Pareto front size: {mask.sum()} / {len(sub)}  ({mask.mean()*100:.1f}%)")
pf = sub[mask].sort_values("pIC50_median", ascending=False)
print("\n-- Pareto front (sorted by activity) --")
print(pf[["std_smiles","pIC50_median","QED","MW","logP","TPSA"]].head(10).round(3).to_string(index=False))

# 三种不同"最佳"选择
print("\n-- 三种不同排序规则，选出三个不同的 Top1 --")
top_act = sub.sort_values("pIC50_median", ascending=False).iloc[0]
top_qed = sub.sort_values("QED", ascending=False).iloc[0]
sub["naive_sum"] = sub.pIC50_median + sub.QED   # 错误示范：量纲不同直接相加
top_sum = sub.sort_values("naive_sum", ascending=False).iloc[0]
from sklearn.preprocessing import StandardScaler
z = StandardScaler().fit_transform(sub[["pIC50_median","QED"]].values)
sub["z_sum"] = 0.7*z[:,0] + 0.3*z[:,1]
top_z = sub.sort_values("z_sum", ascending=False).iloc[0]
for name, row in [("Top activity", top_act), ("Top QED", top_qed),
                  ("naive sum (量纲未处理)", top_sum), ("z-score 加权", top_z)]:
    print(f"  {name:<22s} pIC50={row.pIC50_median:.2f}  QED={row.QED:.3f}  MW={row.MW:.0f}  logP={row.logP:.2f}")
print("\n  -> naive sum 与 z-score 选出的分子是否同一个: ",
      top_sum.std_smiles == top_z.std_smiles)
print(f"  QED 与 pIC50 的相关性 (Spearman): ", end="")
from scipy.stats import spearmanr
print(f"{spearmanr(sub.pIC50_median, sub.QED).statistic:.3f}")

print("\n" + "=" * 72)
print("[3] 脏数据值多少钱：同一个测试集，训练集带不带异常值")
print("=" * 72)
from sklearn.model_selection import GroupShuffleSplit
qc_pos = np.where(keep_qc)[0]          # QC 通过的行号
dirty_pos = np.where(~keep_qc)[0]
splitter = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=42)
tr_s, te_s = next(splitter.split(agg.iloc[qc_pos], groups=agg.iloc[qc_pos]["scaffold"]))
te = qc_pos[te_s]                       # 测试集：干净分子
tr = qc_pos[tr_s]
dirty_train = np.setdiff1d(dirty_pos, te)
print(f"test(all clean)={len(te)}  train(clean)={len(tr)}  dirty extra={len(dirty_train)}")

def rf_mae(train_idx):
    m = RandomForestRegressor(n_estimators=300, random_state=42, n_jobs=-1)
    m.fit(X[train_idx], y[train_idx])
    return mean_absolute_error(y[te], m.predict(X[te]))

mae_clean = rf_mae(tr)
mae_dirty = rf_mae(np.concatenate([tr, dirty_train]))
print(f"  train clean only          : MAE={mae_clean:.3f}")
print(f"  train clean + 脏数据({len(dirty_train)}) : MAE={mae_dirty:.3f}")
print(f"  -> 脏数据占比 {len(dirty_train)/(len(tr)+len(dirty_train))*100:.2f}%，MAE 变化 {mae_dirty-mae_clean:+.3f}")

print("\n" + "=" * 72)
print("[4] Active learning: uncertainty sampling vs random vs greedy")
print("=" * 72)
rng = np.random.RandomState(0)
pool0 = rng.permutation(tr)
seed_n = 300
pool = pool0[seed_n:]
seed = pool0[:seed_n]
print(f"test={len(te)}  initial labeled={len(seed)}  pool={len(pool)}")

def fit_eval(train_idx, label=None):
    """返回 (test MAE, test R2, model)；不确定性用森林内部各棵树的预测方差，一次训练即可"""
    m = RandomForestRegressor(n_estimators=200, random_state=42, n_jobs=-1)
    m.fit(X[train_idx], y[train_idx])
    p = m.predict(X[te])
    return mean_absolute_error(y[te], p), r2_score(y[te], p), m

def ensemble_uncertainty(model, pool_idx):
    P = np.stack([t.predict(X[pool_idx]) for t in model.estimators_])  # (n_trees, n_pool)
    return P.mean(axis=0), P.std(axis=0)

history = {"uncertainty": [], "random": [], "greedy": []}
BATCH = 300
N_ROUNDS = 6
T0 = time.time()

for strategy in ["uncertainty", "random", "greedy"]:
    train_idx = list(seed)
    pool_idx = list(pool)
    print(f"\n  -- {strategy} --")
    for r in range(N_ROUNDS):
        mae, r2, model = fit_eval(train_idx)
        history[strategy].append({"round": r, "n_labeled": len(train_idx),
                                  "mae": round(mae, 3), "r2": round(r2, 3)})
        print(f"    round {r}: n={len(train_idx):<5d} MAE={mae:.3f} R2={r2:.3f}   ({time.time()-T0:.0f}s)", flush=True)
        if len(pool_idx) <= BATCH:
            break
        mean, std = ensemble_uncertainty(model, pool_idx)
        if strategy == "random":
            pick = rng.choice(len(pool_idx), BATCH, replace=False)
        elif strategy == "uncertainty":
            pick = np.argsort(-std)[:BATCH]
        else:  # greedy: 只挑预测活性最高的（模型最自信的地方）
            pick = np.argsort(-mean)[:BATCH]
        pick = np.asarray(pick)
        train_idx = train_idx + [pool_idx[i] for i in pick]
        keep = np.ones(len(pool_idx), dtype=bool); keep[pick] = False
        pool_idx = [pool_idx[i] for i in np.where(keep)[0]]

print("\n-- summary: MAE by labeled-set size --")
n_grid = sorted({h["n_labeled"] for s in history for h in history[s]})
rows = []
for n in n_grid:
    row = {"n_labeled": n}
    for s in history:
        match = [h for h in history[s] if h["n_labeled"] == n]
        if match:
            row[s] = match[0]["mae"]
    rows.append(row)
print(pd.DataFrame(rows).to_string(index=False))

# 用 random 需要多少数据才能达到 uncertainty 在 1800 时的 MAE
unc = history["uncertainty"]; rnd = history["random"]
target = [h for h in unc if h["n_labeled"] == 1800]
if target:
    t = target[0]["mae"]
    hit = [h for h in rnd if h["mae"] <= t]
    print(f"\nuncertainty @1800 达到 MAE={t:.3f}；random 达到同等水平需要 "
          f"{hit[0]['n_labeled'] if hit else '>本次实验上限'} 个样本")

json.dump(history, open(os.path.join(DATA, "al_history.json"), "w"), indent=2)
print("\nsaved -> al_history.json / dataset_with_props.csv")
