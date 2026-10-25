"""11 · 干湿闭环：实验结果回来以后怎么办

输入：data/processed/dataset_with_scaffold.csv、split_indices.json、04 的 pIC50 真值
输出：data/processed/iteration_history.csv、11_closed_loop.json

把「已发表的 ChEMBL 活性数据」当作**真实湿实验结果**来用：
    每轮由策略挑一批分子 → 查它们真实的 pIC50 当作回来了的实验值 → 加入训练集重训，
    看哪个策略用更少的实验预算把模型推得更好。

三种策略：
    uncertainty  挑模型「最没把握」的样本
    greedy      只挑模型预测活性最高的样本
    random      随机采样（对照基线）

跑完一次大约需要 15~20 分钟。若只想基于已有的 iteration_history.csv 重算结论，
加环境变量 AIDD_CL_REPORT_ONLY=1，可跳过全部重训练。
"""
import os
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, r2_score

from _common import cfg, P, ensure_dirs, processed, path, log, banner, save_json

ensure_dirs()
banner("11 · 干湿闭环：同一种预算，挑谁去测，结果差很多")

C = cfg["closed_loop"]
SEED = int(cfg["random_seed"])
N_ROUNDS = int(C["n_rounds"])
BATCH = int(C["batch_size"])
INIT = int(C["initial_labeled"])
STRATS = C["strategies"]

M = cfg["model"]
FP_SIZE = int(M["fingerprint_size"])
FP_RADIUS = int(M["fingerprint_radius"])

# ---------------- 1. 特征与标签 ----------------
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit import DataStructs

RDLogger.DisableLog("rdApp.*")

full = pd.read_csv(processed("dataset_with_scaffold.csv"))
splits = __import__("json").load(open(processed("split_indices.json")))
te = np.array(splits["scaffold"]["test"])
tr = np.array(splits["scaffold"]["train"])
log(f"沿用 07 的 scaffold split：train={len(tr)} test={len(te)}")

gen = rdFingerprintGenerator.GetMorganGenerator(radius=FP_RADIUS, fpSize=FP_SIZE)


def featurize(smiles):
    mol = Chem.MolFromSmiles(smiles)
    arr = np.zeros((FP_SIZE,), dtype=np.uint8)
    if mol is None:
        return arr
    DataStructs.ConvertToNumpyArray(gen.GetFingerprint(mol), arr)
    return arr


log("生成 ECFP4 特征…")
X = np.stack([featurize(s) for s in full["smiles"]])
y = full["pIC50_median"].values
log(f"X shape = {X.shape}，y 均值 = {y.mean():.3f}")


# ---------------- 2. 一次训练同时得到预测与不确定性 ----------------
def rf_uncertainty(train_idx, pool_idx, test_idx, n_models=3):
    """用 n_models 个不同随机种子的 RF 组成 ensemble，std 即不确定性"""
    preds, tests = [], []
    for s in range(n_models):
        m = RandomForestRegressor(n_estimators=int(M["rf_n_estimators"]) // 2,
                                  random_state=SEED + s, n_jobs=1)
        m.fit(X[train_idx], y[train_idx])
        preds.append(m.predict(X[pool_idx]))
        tests.append(m.predict(X[test_idx]))
    Pm = np.stack(preds)
    Tm = np.stack(tests)
    return Pm.mean(0), Pm.std(0), Tm.mean(0)


# ---------------- 3. 三策略跑同一段轮次 ----------------
if os.environ.get("AIDD_CL_REPORT_ONLY") == "1":
    # 只重算结论：跳过全部重训练（跑一轮要 15~20 分钟）
    banner("跳过重训练，直接基于已有 history 重算结论")
    hist_df = pd.read_csv(processed("iteration_history.csv"))
    strategy_order = {s: i for i, s in enumerate(STRATS)}
    history = {}
    for s, g in hist_df.groupby("strategy"):
        g = g.sort_values("round")
        history[s] = g[["round", "n_labeled", "MAE", "R2"]].to_dict("records")
    hist_df["_o"] = hist_df["strategy"].map(strategy_order)
    hist_df = hist_df.sort_values(["_o", "round"]).drop(columns="_o")
else:
    history = {s: [] for s in STRATS}
rng = np.random.RandomState(SEED)
pool0 = rng.permutation(tr)          # 只在训练集内模拟「还有哪些分子没测」

if os.environ.get("AIDD_CL_REPORT_ONLY") != "1":
    for strategy in STRATS:
        banner(f"策略：{strategy}")
        labeled = list(pool0[:INIT])
        pool = list(pool0[INIT:])
        for r in range(N_ROUNDS):
            _, _, pred_test = rf_uncertainty(np.array(labeled), np.array(labeled[:1]), te)
            mae = mean_absolute_error(y[te], pred_test)
            r2 = r2_score(y[te], pred_test)
            log(f"  round {r}: n={len(labeled):<5d} MAE={mae:.4f} R2={r2:.4f}")
            history[strategy].append({"round": r, "n_labeled": len(labeled),
                                      "MAE": round(float(mae), 4), "R2": round(float(r2), 4)})
            if r == N_ROUNDS - 1 or len(pool) < BATCH:
                break
            if strategy == "random":
                pick_pos = rng.choice(len(pool), BATCH, replace=False)
            elif strategy == "uncertainty":
                _, std, _ = rf_uncertainty(np.array(labeled), np.array(pool), te)
                pick_pos = np.argsort(-std)[:BATCH]
            else:  # greedy
                mean, _, _ = rf_uncertainty(np.array(labeled), np.array(pool), te, n_models=1)
                pick_pos = np.argsort(-mean)[:BATCH]
            picked = [pool[i] for i in np.asarray(pick_pos)]
            labeled.extend(picked)
            remain = set(range(len(pool))) - set(np.asarray(pick_pos))
            pool = [pool[i] for i in sorted(remain)]

    hist_df = pd.DataFrame([dict(rows, strategy=s) for s, rs in history.items() for rows in rs])
    hist_df.to_csv(processed("iteration_history.csv"), index=False)
    log("轮次历史已写入 data/processed/iteration_history.csv")


# ---------------- 4. 结论：不要只看终点，要看整条曲线 ----------------
banner("对照结论")
pivot = hist_df.pivot_table(index="n_labeled", columns="strategy", values="MAE")
pivot = pivot[[c for c in STRATS if c in pivot.columns]]
print("\n同预算下的测试集 MAE（越低越好）：")
print(f"\n{'样本数':>8}" + "".join(f"{c:>14}" for c in pivot.columns))
print("-" * (8 + 14 * len(pivot.columns)))
for n, row in pivot.iterrows():
    print(f"{int(n):>8}" + "".join(f"{row[c]:>14.4f}" for c in pivot.columns))

summary = {}
for s in pivot.columns:
    h = hist_df[hist_df.strategy == s].sort_values("round")
    summary[s] = {"start_MAE": float(h.MAE.iloc[0]), "final_MAE": float(h.MAE.iloc[-1]),
                  "final_R2": float(h.R2.iloc[-1]), "final_n": int(h.n_labeled.iloc[-1])}
best = min(summary, key=lambda s: summary[s]["final_MAE"])
worst = max(summary, key=lambda s: summary[s]["final_MAE"])
print(f"\n终点最好的策略：{best}（MAE {summary[best]['final_MAE']:.4f}）"
      f"　最差：{worst}（MAE {summary[worst]['final_MAE']:.4f}）")

# 关键：uncertainty 什么时候才真正超过 random？（如果全程没超过，就要实话实说）
cross = None
if "uncertainty" in pivot.columns and "random" in pivot.columns:
    for n, row in pivot.iterrows():
        if row["uncertainty"] < row["random"]:
            cross = int(n)
            break
    if cross is not None:
        log(f"uncertainty 从第 {cross} 个样本起才反超 random —— "
            f"在更早的阶段它反而落后")
    else:
        log("在本轮预算内，uncertainty 全程没有反超 random", "WARN")

# greedy 为什么会「越挑越差」：它只挑模型已经很自信且高预测值的区域
if "greedy" in summary:
    g0, g1 = summary["greedy"]["start_MAE"], summary["greedy"]["final_MAE"]
    log(f"greedy：从 {g0:.4f} 到 {g1:.4f}（{(g1-g0)/g0*100:+.1f}%）。"
        f"它只挑模型最自信的高分区，等于反复告诉模型同一件事 —— "
        f"这就是为什么「顺着模型的偏好一路挑下去」往往会越挑越窄。")

save_json({"config": {"rounds": N_ROUNDS, "batch": BATCH, "initial": INIT},
           "summary": summary,
           "cross_over_n": cross,
           "best_final": best, "worst_final": worst,
           "mae_by_budget": pivot.round(4).to_dict(),
           "history": history}, processed("11_closed_loop.json"))
banner("11 done · 闭环的价值在于：实验结果不是终点，而是下一轮的起点")
