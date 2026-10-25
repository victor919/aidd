# -*- coding: utf-8 -*-
# 文章关联：AI制药实践篇_公众号文章.html
# 脚本职责：实验：第七站，主动学习与随机选样对照
# 使用注解：TDC Caco-2 真实标签用于离线模拟回传；没有执行湿实验。
# 完整对应关系与运行顺序：README.md
"""干湿闭环迷你演示：主动学习 vs 随机选样（真实数据 TDC Caco-2 / 真实模型 RandomForest）"""
import warnings, numpy as np, pandas as pd
warnings.filterwarnings("ignore")
from rdkit import Chem, RDLogger, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score, mean_squared_error
from tdc.single_pred import ADME
RDLogger.DisableLog('rdApp.*')

data = ADME(name="Caco2_Wang").get_data()
print("数据集 Caco2_Wang：%d 个分子，标签 = Caco-2 渗透性 logPapp（cm/s）" % len(data))
print("标签范围：%.2f ~ %.2f，均值 %.2f" % (data.Y.min(), data.Y.max(), data.Y.mean()))

gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
X, y, ok = [], [], []
for smi, val in zip(data.Drug, data.Y):
    m = Chem.MolFromSmiles(smi)
    if m is None: continue
    X.append(np.array(gen.GetFingerprint(m), dtype=np.float64)); y.append(val)
X = np.array(X); y = np.array(y)
print("可解析分子：%d 个，指纹维度：%d" % (X.shape[0], X.shape[1]))

def run(seed, strategy, n_init=100, n_batch=50, n_round=6):
    rng = np.random.default_rng(seed)
    idx = rng.permutation(len(y))
    test_idx = idx[:180]              # 固定测试集
    pool_all = idx[180:]
    start = rng.permutation(pool_all)[:n_init]
    pool = np.array([i for i in pool_all if i not in start])
    labeled = list(start)
    hist = []
    for r in range(n_round + 1):
        model = RandomForestRegressor(n_estimators=300, random_state=0, n_jobs=-1)
        model.fit(X[labeled], y[labeled])
        pred = model.predict(X[test_idx])
        hist.append((len(labeled), r2_score(y[test_idx], pred),
                     float(np.sqrt(mean_squared_error(y[test_idx], pred)))))
        if r == n_round: break
        if strategy == "random":
            pick = rng.permutation(pool)[:n_batch]
        else:  # uncertainty：用森林各树预测的标准差挑选最不确定的
            preds = np.array([t.predict(X[pool]) for t in model.estimators_])
            unc = preds.std(axis=0)
            pick = pool[np.argsort(-unc)[:n_batch]]
        labeled.extend(pick.tolist())
        pool = np.array([i for i in pool if i not in set(pick)])
    return hist

results = {"random": [], "uncertainty": []}
for strat in results:
    for seed in [0, 1, 2]:
        results[strat].append(run(seed, strat))

print()
print("=" * 88)
print("表 7｜干湿闭环迷你实验：同样做 400 次“实验”，选样策略不同，模型差多少")
print("=" * 88)
print("%-14s %14s %14s %14s" % ("已标注样本数", "随机选样 R²", "主动学习 R²", "R² 差值"))
for r in range(7):
    n = results["random"][0][r][0]
    rr = np.mean([h[r][1] for h in results["random"]])
    ur = np.mean([h[r][1] for h in results["uncertainty"]])
    print("%-14d %14.3f %14.3f %14.3f" % (n, rr, ur, ur - rr))

print()
print("RMSE（三次随机种子平均）")
print("%-14s %14s %14s" % ("已标注样本数", "随机选样", "主动学习"))
for r in range(7):
    n = results["random"][0][r][0]
    rr = np.mean([h[r][2] for h in results["random"]])
    ur = np.mean([h[r][2] for h in results["uncertainty"]])
    print("%-14d %14.3f %14.3f" % (n, rr, ur))

# 达到相同精度，主动学习省了多少次实验
target = 0.60
def first_hit(strat, col=1, thr=target):
    for r in range(7):
        v = np.mean([h[r][col] for h in results[strat]])
        if v >= thr: return results[strat][0][r][0], v
    return None, None
n_r, v_r = first_hit("random"); n_u, v_u = first_hit("uncertainty")
print()
if n_r and n_u:
    print("测试集 R² 达到 %.2f：随机选样需要 %d 个样本，主动学习只需 %d 个 → 少做 %d 次“实验”（省 %.0f%%）"
          % (target, n_r, n_u, n_r - n_u, 100 * (n_r - n_u) / n_r))
else:
    print("三次种子平均后，随机选样与主动学习在 400 样本内达到 R²=%.2f 的情况：" % target, n_r, n_u)
