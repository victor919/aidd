# -*- coding: utf-8 -*-
"""
AI 制药工具篇 · 第五节配套脚本：ECFP / Morgan 分子指纹
------------------------------------------------------
环境：conda activate aidd
      Python 3.11.14 / RDKit 2026.03.5 / scikit-learn 1.9.0
运行：python aidd_ecfp_demo.py
产出：molecules.csv（描述符表）、ecfp_X.npy（指纹矩阵）、
      ecfp_similarity_matrix.png（相似度热力图）

脚本按 ①→⑧ 顺序执行，每一步都可以单独拿出来塞进 Notebook 里跑。
"""
import time
import numpy as np
import pandas as pd

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem import Descriptors, Crippen, QED

RDLogger.DisableLog("rdApp.*")   # 关掉 RDKit 的啰嗦日志

# 10 个熟悉的小分子药 + 1 个用于演示相似性的水杨酸
DRUGS = [
    ("阿司匹林",     "CC(=O)Oc1ccccc1C(=O)O"),
    ("布洛芬",       "CC(C)Cc1ccc(C(C)C(=O)O)cc1"),
    ("对乙酰氨基酚", "CC(=O)Nc1ccc(O)cc1"),
    ("咖啡因",       "Cn1c(=O)c2c(ncn2C)n(C)c1=O"),
    ("二甲双胍",     "CN(C)C(=N)NC(=N)N"),
    ("尼古丁",       "CN1CCC[C@H]1c1cccnc1"),
    ("维生素C",      "O=C1O[C@H]([C@H](O)CO)C(O)=C1O"),
    ("地西泮",       "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21"),
    ("西地那非",     "CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)nc12"),
    ("阿莫西林",     "CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccc(O)cc3)C(=O)N2C1C(=O)O"),
    ("水杨酸",       "O=C(O)c1ccccc1O"),
]


def line(title):
    print("\n" + "=" * 58)
    print(title)
    print("=" * 58)


# ---------------------------------------------------------------
# ① SMILES → Mol（指纹之前必须先有 Mol 对象）
# ---------------------------------------------------------------
line("① SMILES → Mol")
mols, names = {}, []
for name, smi in DRUGS:
    m = Chem.MolFromSmiles(smi)
    if m is None:                      # 永远先判 None，别让脏数据溜到后面
        print(f"  [跳过] {name} 解析失败：{smi}")
        continue
    mols[name] = m
    names.append(name)
print(f"  成功解析 {len(mols)} 个分子：{' / '.join(names)}")

mol = mols["阿司匹林"]
print(f"  阿司匹林：{mol.GetNumAtoms()} 个重原子，"
      f"canonical SMILES = {Chem.MolToSmiles(mol)}")

# ---------------------------------------------------------------
# ② 生成 ECFP4：现在推荐的写法是 rdFingerprintGenerator
# ---------------------------------------------------------------
line("② 生成 ECFP4 指纹（radius=2, 2048 位）")
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
fp = gen.GetFingerprint(mol)

print(f"  type(fp)          = {type(fp).__name__}")
print(f"  len(fp)           = {fp.GetNumBits()} 位")
print(f"  fp.GetNumOnBits() = {fp.GetNumOnBits()} 位被点亮")

bits = list(fp.GetOnBits())            # 被点亮的位下标
print(f"  点亮的位下标        = {bits}")
print(f"  前 40 位长什么样    = {''.join('1' if fp.GetBit(i) else '0' for i in range(40))}")

# 说明：1 万个分子里绝大多数位永远是 0，所以存成 0/1 向量非常稀疏

# 旧写法（大量老教程里还是这个，功能等价，见到要认得）
fp_old = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
print(f"  [旧 API] 点亮位数    = {fp_old.GetNumOnBits()}")
print(f"  新旧写法是否一致      = {list(fp.GetOnBits()) == list(fp_old.GetOnBits())}")

# ---------------------------------------------------------------
# ③ radius 实验：一圈一圈往外看，看看到底点亮多少位
# ---------------------------------------------------------------
line("③ radius 实验（阿司匹林，2048 位）")
for r in range(4):
    g = rdFingerprintGenerator.GetMorganGenerator(radius=r, fpSize=2048)
    f = g.GetFingerprint(mol)
    tag = "   ← ECFP4，工业界默认" if r == 2 else ""
    print(f"  radius={r}  点亮 {f.GetNumOnBits():>3} 位  "
          f"稀疏度 {f.GetNumOnBits() / 2048 * 100:5.2f}%{tag}")

# ---------------------------------------------------------------
# ④ 位宽实验：2048 位够不够？位宽变小会发生"撞车"
# ---------------------------------------------------------------
line("④ 位宽实验（阿司匹林，radius=2）")
for n in (128, 512, 1024, 2048, 4096):
    g = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=n)
    f = g.GetFingerprint(mol)
    print(f"  nBits={n:>5}  点亮 {f.GetNumOnBits():>3} 位  "
          f"稀疏度 {f.GetNumOnBits() / n * 100:5.2f}%")

# ---------------------------------------------------------------
# ⑤ 速度：化学信息学对硬件有多友好
# ---------------------------------------------------------------
line("⑤ 速度测试")
N = 10000
smi_list = [s for _, s in DRUGS] * (N // len(DRUGS) + 1)
smi_list = smi_list[:N]
m_list = [Chem.MolFromSmiles(s) for s in smi_list]

t0 = time.perf_counter()
_ = [gen.GetFingerprint(m) for m in m_list]
dt = time.perf_counter() - t0
print(f"  {N} 个分子生成 ECFP4(2048bit)：{dt:.3f} 秒 "
      f"（约 {N / dt:.0f} 个/秒）")

# ---------------------------------------------------------------
# ⑥ Tanimoto 相似度：两个分子到底有多像
# ---------------------------------------------------------------
line("⑥ Tanimoto 相似度（以阿司匹林为查询分子）")
fps = {n: gen.GetFingerprint(m) for n, m in mols.items()}
q = fps["阿司匹林"]

rank = sorted(fps.items(), key=lambda kv: -DataStructs.TanimotoSimilarity(q, kv[1]))
for n, f in rank:
    print(f"  {n:<8} {DataStructs.TanimotoSimilarity(q, f):.3f}")

# 批量版：一次算 1 个 vs N 个，比循环快得多
sims = DataStructs.BulkTanimotoSimilarity(q, list(fps.values()))
print(f"\n  BulkTanimotoSimilarity 一次算完 {len(sims)} 个："
      f"{[round(s, 3) for s in sims]}")

# ---------------------------------------------------------------
# ⑦ 相似度矩阵 + 热力图（这就是最朴素的"虚拟筛选"）
# ---------------------------------------------------------------
line("⑦ 相似度矩阵（11 × 11）")
n = len(names)
sim = np.zeros((n, n))
for i in range(n):
    # 对角线以下复用对称结果，省一半计算
    for j in range(i, n):
        s = DataStructs.TanimotoSimilarity(fps[names[i]], fps[names[j]])
        sim[i, j] = sim[j, i] = s

sim_df = pd.DataFrame(np.round(sim, 3), index=names, columns=names)
print(sim_df.to_string())

# 找出除自己以外最相似的一对（先看上三角，排除对角线）
iu = np.triu_indices(n, k=1)
idx = np.argmax(sim[iu])
i, j = iu[0][idx], iu[1][idx]
print(f"\n  最相似的一对：{names[i]} ↔ {names[j]}  Tanimoto = {sim[i, j]:.3f}")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["font.sans-serif"] = ["PingFang SC"]
    plt.rcParams["axes.unicode_minus"] = False

    fig, ax = plt.subplots(figsize=(6.6, 5.6), dpi=150)
    im = ax.imshow(sim, cmap="Greens", vmin=0, vmax=1)
    ax.set_xticks(range(n)); ax.set_xticklabels(names, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(n)); ax.set_yticklabels(names, fontsize=9)
    for a in range(n):
        for b in range(n):
            ax.text(b, a, f"{sim[a, b]:.2f}", ha="center", va="center",
                    fontsize=7.5, color="white" if sim[a, b] > 0.5 else "#333")
    ax.set_title("Tanimoto 相似度矩阵（ECFP4, 2048bit）", fontsize=12, pad=12)
    fig.colorbar(im, shrink=0.8)
    fig.tight_layout()
    fig.savefig("ecfp_similarity_matrix.png")
    plt.close(fig)
    print("  已保存热力图：ecfp_similarity_matrix.png")
except Exception as e:
    print(f"  [热力图跳过] {e}")

# ---------------------------------------------------------------
# ⑧ 指纹 → 机器学习特征矩阵 X（ECFP + Random Forest 的入口）
# ---------------------------------------------------------------
line("⑧ 指纹 → 特征矩阵 X → sklearn")
X = np.zeros((n, 2048), dtype=np.int8)
for k, nm in enumerate(names):
    DataStructs.ConvertToNumpyArray(fps[nm], X[k])
print(f"  X.shape = {X.shape}   dtype = {X.dtype}")
print(f"  非零元素占比 = {(X != 0).mean() * 100:.2f}%   （非常稀疏，树模型最爱）")
print(f"  第 0 行前 30 维 = {X[0][:30]}")

# 顺手把描述符表也存下来，做课后项目用
rows = []
for nm in names:
    m = mols[nm]
    rows.append((nm, Chem.MolToSmiles(m),
                 round(Descriptors.MolWt(m), 1),
                 round(Crippen.MolLogP(m), 2),
                 round(Descriptors.TPSA(m), 1),
                 round(QED.qed(m), 2)))
df = pd.DataFrame(rows, columns=["name", "smiles", "MW", "logP", "TPSA", "QED"])
df.to_csv("molecules.csv", index=False)
np.save("ecfp_X.npy", X)
print("\n  已保存：molecules.csv（描述符表）、ecfp_X.npy（指纹矩阵）")

# ---------------------------------------------------------------
# ⑨ 选做：在真实数据集上跑一遍 ECFP + Random Forest baseline
#    （需要 deepchem，会自动下载 Delaney/ESOL 数据集）
# ---------------------------------------------------------------
line("⑨ 选做：Delaney 溶解度数据集上的 ECFP + RF baseline")
try:
    import deepchem as dc
    from deepchem.molnet import load_delaney
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.metrics import r2_score, mean_squared_error

    tasks, datasets, _ = load_delaney(featurizer="ECFP", splitter="scaffold",
                                      transformer=[], reload=False)
    train, valid, test = datasets
    rf = RandomForestRegressor(n_estimators=500, random_state=0, n_jobs=-1)
    rf.fit(train.X, train.y.ravel())
    for split_name, d in (("Train", train), ("Test", test)):
        pred = rf.predict(d.X)
        print(f"  {split_name:<5} R²   = {r2_score(d.y.ravel(), pred):.3f}")
        print(f"  {split_name:<5} RMSE = {mean_squared_error(d.y.ravel(), pred) ** 0.5:.3f}")
    print("  （对照：换成 splitter='random'，Test R² 会明显更高——见第九节）")
except Exception as e:
    print(f"  [跳过] {type(e).__name__}: {e}")

print("\n全部跑完 ✅")
