# 文章关联：AI制药实践篇B_实操篇_公众号文章.html
# 脚本职责：反例实验与配图：硬编码盒子、演示库、SA 与 docking
# 使用注解：读取正式 aidd-mvp 的配置和产物，写 mvp_trap_results.json 及 t1–t3。
# 完整对应关系与运行顺序：README.md
"""验证一份流传的 AIDD MVP 伪代码：三处硬伤 + 两条真实数据的反面证据。

素材来自 @long-text 里的 AIDD_MVP_Pipeline 伪代码：
    ① pocket_center 硬编码 {"x":12.5,"y":-8.3,"z":24.1}
    ② demo_library = [Aspirin, Caffeine, Ibuprofen] 去跑 EGFR pipeline
    ③ 用 MolWt<500 / QED>0.5 充当活性筛选
再用真实数据看两条避坑指南：
    ④ SAScore 是否应该在 Day-1 就进筛选
    ⑤ docking score 与真实活性到底有没有关系

所有数字本机实跑，数据源：aidd-mvp 工程（config.yaml / seed=42）。
"""
import os
import json
import time
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Descriptors, QED, rdFingerprintGenerator
from rdkit import RDConfig
from scipy.stats import spearmanr, pearsonr
from sklearn.ensemble import RandomForestRegressor

RDLogger.DisableLog("rdApp.*")

ROOT = os.path.expanduser("~/WorkBuddy/Work/aidd-mvp")
PROC = os.path.join(ROOT, "data", "processed")
INTER = os.path.join(ROOT, "data", "interim")
RAW = os.path.join(ROOT, "data", "raw")
OUT_DIR = os.path.expanduser("~/WorkBuddy/Work/aidd_practice")
FIG_DIR = os.path.expanduser("~/WorkBuddy/Work/AI制药实操篇_配图")
os.makedirs(OUT_DIR, exist_ok=True)

import yaml
CFG = yaml.safe_load(open(os.path.join(ROOT, "config.yaml"), encoding="utf-8"))
M = CFG["model"]
FP_R, FP_SIZE = int(M["fingerprint_radius"]), int(M["fingerprint_size"])
SEED = int(CFG["random_seed"])
R = {}
T0 = time.time()


def log(msg):
    print(f"[{time.time()-T0:7.1f}s] {msg}", flush=True)


def banner(t):
    print("\n" + "=" * 66 + f"\n{t}\n" + "=" * 66, flush=True)


# ============================================================ Part A 盒子
banner("A · 硬编码的口袋坐标，落在蛋白的什么地方？")

# 手工解析 1M17（不依赖 Biopython）
prot_xyz, aq4_xyz = [], []
pocket_res = set()
for line in open(os.path.join(RAW, "1M17.pdb")):
    if line.startswith(("ATOM", "HETATM")):
        xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        resname = line[17:20].strip()
        if resname == "AQ4" and line.startswith("HETATM"):
            aq4_xyz.append(xyz)
        elif line.startswith("ATOM"):
            prot_xyz.append(xyz)
prot_xyz = np.array(prot_xyz)
aq4_xyz = np.array(aq4_xyz)
log(f"1M17 蛋白原子 {len(prot_xyz)} 个，配体 AQ4 原子 {len(aq4_xyz)} 个")

# pocket.json 里的 20 个口袋残基
pkj = json.load(open(os.path.join(PROC, "pocket.json")))
pocket_resnumbers = {r["pdb_num"] for r in pkj["pocket_residues"]}

true_center = aq4_xyz.mean(axis=0)
hard_center = np.array([12.5, -8.3, 24.1])      # 伪代码里写死的那个
BOX = 20.0
half = BOX / 2.0


def in_box(xyz, center, half):
    return (np.abs(xyz - center) <= half).all(axis=1)


d_center = float(np.linalg.norm(true_center - hard_center))
lig_in_hard = int(in_box(aq4_xyz, hard_center, half).sum())
lig_in_true = int(in_box(aq4_xyz, true_center, half).sum())

# 口袋残基 / 蛋白原子在两个盒子里的覆盖情况
resatoms = {}
for line in open(os.path.join(RAW, "1M17.pdb")):
    if line.startswith("ATOM"):
        xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        resname = line[17:20].strip()
        try:
            resnum = int(line[22:26])
        except ValueError:
            continue
        if resnum in pocket_resnumbers:
            resatoms.setdefault(resnum, []).append(xyz)

cov = {"hard": 0, "true": 0}
for resnum, atoms in resatoms.items():
    arr = np.array(atoms)
    cov["hard"] += int(in_box(arr, hard_center, half).any())
    cov["true"] += int(in_box(arr, true_center, half).any())

prot_in_hard = int(in_box(prot_xyz, hard_center, half).sum())
prot_in_true = int(in_box(prot_xyz, true_center, half).sum())

R["box"] = {
    "true_center": [round(float(x), 3) for x in true_center],
    "hard_coded_center": [12.5, -8.3, 24.1],
    "center_distance": round(d_center, 2),
    "ligand_atoms": int(len(aq4_xyz)),
    "ligand_atoms_in_hard_box": lig_in_hard,
    "ligand_atoms_in_true_box": lig_in_true,
    "pocket_residues_total": len(resatoms),
    "pocket_residues_in_hard_box": cov["hard"],
    "pocket_residues_in_true_box": cov["true"],
    "protein_atoms_in_hard_box": prot_in_hard,
    "protein_atoms_in_true_box": prot_in_true,
}
print(f"  真实口袋中心 {np.round(true_center,2)} vs 硬编码 {hard_center} → 相距 {d_center:.1f} 埃")
print(f"  配体 29 个原子落进盒子：硬编码 {lig_in_hard} 个 / 正确 {lig_in_true} 个")
print(f"  20 个口袋残基覆盖：硬编码 {cov['hard']} 个 / 正确 {cov['true']} 个")
print(f"  盒内蛋白原子数：硬编码 {prot_in_hard} / 正确 {prot_in_true}")

# ============================================================ Part B demo 分子库
banner("B · 用阿司匹林/咖啡因/布洛芬去跑 EGFR pipeline，会输出什么？")

df = pd.read_csv(os.path.join(PROC, "dataset.csv"))
X = np.load(os.path.join(PROC, "X_ecfp.npy"))
y = df["pIC50_median"].values
splits = json.load(open(os.path.join(PROC, "split_indices.json")))
tr = np.array(splits["scaffold"]["train"])
log(f"数据 {X.shape}，标签 [{y.min():.2f}, {y.max():.2f}]，骨架切分训练 {len(tr)} 条")

rf = RandomForestRegressor(n_estimators=int(M["rf_n_estimators"]),
                           random_state=SEED, n_jobs=-1)
rf.fit(X[tr], y[tr])
log("RF（ECFP4 + 500 棵，骨架切分训练集）训练完成")

gen = rdFingerprintGenerator.GetMorganGenerator(radius=FP_R, fpSize=FP_SIZE)


def fp_of(smiles):
    mol = Chem.MolFromSmiles(smiles)
    arr = np.zeros((FP_SIZE,), dtype=np.uint8)
    if mol is None:
        return None, arr
    DataStructs.ConvertToNumpyArray(gen.GetFingerprint(mol), arr)
    return mol, arr


def props(mol):
    return dict(MW=round(Descriptors.MolWt(mol), 1),
                logP=round(Descriptors.MolLogP(mol), 2),
                QED=round(QED.qed(mol), 3))


demo = [
    ("Aspirin（伪代码 demo）", "CC(=O)Oc1ccccc1C(=O)O"),
    ("Caffeine（伪代码 demo）", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"),
    ("Ibuprofen（伪代码 demo）", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"),
]

# 对照组 1：1M17 共晶配体（从真实晶体文件还原的化学结构，不是手写 SMILES）
ctrl = []
try:
    from rdkit.Chem import rdMolDescriptors  # noqa
    pdbmol = Chem.MolFromPDBFile(os.path.join(INTER, "aq4_crystal.pdb"), removeHs=False)
    log(f"从 aq4_crystal.pdb 读入 {pdbmol.GetNumAtoms()} 个原子")
except Exception as e:  # pragma: no cover
    pdbmol = None
    log(f"读取晶体配体失败：{e}")

# 对照组 2：ChEMBL 里叫 Erlotinib 的那个分子（联网取官方 canonical SMILES）
erlo_smiles, erlo_source = None, None
try:
    import requests
    for attempt in range(3):
        try:
            r = requests.get(
                "https://www.ebi.ac.uk/chembl/api/data/molecule.json",
                params={"pref_name": "ERLOTINIB", "limit": 1}, timeout=30)
            if r.status_code == 200:
                j = r.json()["molecules"][0]["molecule_structures"]
                erlo_smiles = j["canonical_smiles"]
                erlo_source = "ChEMBL API (pref_name=ERLOTINIB)"
                break
        except Exception:
            time.sleep(3)
except Exception as e:  # pragma: no cover
    log(f"ChEMBL 查询失败：{e}")

rows = []
for name, smi in demo:
    mol, arr = fp_of(smi)
    pred = float(rf.predict(arr.reshape(1, -1))[0])
    rows.append(dict(name=name, group="demo", smiles=smi,
                     pred_pIC50=round(pred, 2), **props(mol)))

if erlo_smiles:
    mol, arr = fp_of(erlo_smiles)
    pred = float(rf.predict(arr.reshape(1, -1))[0])
    rows.append(dict(name=f"Erlotinib（{erlo_source}）", group="real",
                     smiles=erlo_smiles, pred_pIC50=round(pred, 2), **props(mol)))

# 对照组 3：数据集里真实高活性的分子（取 pIC50 最高的 3 个不同分子）
top_idx = np.argsort(-y)[:400]
seen, n_add = set(), 0
for i in top_idx:
    smi = df.iloc[i]["smiles"]
    key = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
    if key in seen:
        continue
    seen.add(key)
    mol, arr = fp_of(smi)
    rows.append(dict(name=f"数据集真实高分 #{n_add+1}", group="real", smiles=smi,
                     pred_pIC50=round(float(rf.predict(arr.reshape(1, -1))[0]), 2),
                     true_pIC50=round(float(y[i]), 2), **props(mol)))
    n_add += 1
    if n_add >= 3:
        break

demo_pred = [r["pred_pIC50"] for r in rows if r["group"] == "demo"]
R["demo_library"] = {
    "rows": rows,
    "erlotinib_source": erlo_source,
    "demo_pass_MW500": all(r["MW"] < 500 for r in rows if r["group"] == "demo"),
    "demo_pass_QED05": all(r["QED"] > 0.5 for r in rows if r["group"] == "demo"),
    "demo_pred_mean": round(float(np.mean(demo_pred)), 2),
    "demo_pred_max": round(float(max(demo_pred)), 2),
    "dataset_activity_threshold_pIC50_7_rate": round(float((y >= 7).mean()), 4),
}
print("\n  伪代码 demo 库三件套：")
for r in rows:
    if r["group"] == "demo":
        print(f"    {r['name']:<28} MW={r['MW']:<7} QED={r['QED']:<6} 预测 pIC50={r['pred_pIC50']}")
print("  对照组：")
for r in rows:
    if r["group"] == "real":
        extra = f" 实测={r['true_pIC50']}" if r.get("true_pIC50") else ""
        print(f"    {r['name']:<28} MW={r['MW']:<7} QED={r['QED']:<6} 预测 pIC50={r['pred_pIC50']}{extra}")

# ============================================================ Part C SAScore Day-1
banner("C · SAScore：Day-1 就进筛选，会砍掉什么？")

import sys
sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer  # noqa: E402

sa, ok = [], 0
for smi in df["smiles"]:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        sa.append(np.nan)
        continue
    sa.append(sascorer.calculateScore(mol))
    ok += 1
sa = np.array(sa)
log(f"SAScore 计算完成 {ok} 个分子")

valid = ~np.isnan(sa)
sa_v, y_v = sa[valid], y[valid]
rho_sa, _ = spearmanr(sa_v, y_v)
strong, weak = sa_v[y_v >= 8], sa_v[y_v < 7]

GEN_TH = 4.5   # 生成模型常用阈值
R["sascore"] = {
    "n": int(valid.sum()),
    "mean": round(float(sa_v.mean()), 2),
    "median": round(float(np.median(sa_v))),
    "max": round(float(sa_v.max()), 2),
    "min": round(float(sa_v.min()), 2),
    "pass_le_4_5": round(float((sa_v <= GEN_TH).mean()), 4),
    "pass_le_6_0": round(float((sa_v <= 6.0).mean()), 4),
    "spearman_sa_vs_pIC50": round(float(rho_sa), 4),
    "strong_mean_sa": round(float(strong.mean()), 2),
    "weak_mean_sa": round(float(weak.mean()), 2),
    "strong_n": int(len(strong)),
    "weak_n": int(len(weak)),
}
print(f"  SAScore 均值 {sa_v.mean():.2f}，中位数 {np.median(sa_v):.2f}，最难 {sa_v.max():.2f}")
print(f"  SA<=4.5 通过率 {(sa_v<=4.5).mean():.1%}；SA<=6.0 通过率 {(sa_v<=6.0).mean():.1%}")
print(f"  SA vs pIC50 Spearman = {rho_sa:.3f}")
print(f"  强活性组(n={len(strong)}) 平均 SA {strong.mean():.2f} vs 弱活性组(n={len(weak)}) {weak.mean():.2f}")

# 漏斗：加 SA 约束前后
act = y >= CFG["admet"]["activity_min"]
funnel_noSA = int(act.sum())
funnel_SA45 = int((act & (sa <= GEN_TH)).sum())
funnel_SA60 = int((act & (sa <= 6.0)).sum())
R["sascore"]["funnel"] = {
    "activity_min": CFG["admet"]["activity_min"],
    "n_active": funnel_noSA,
    "n_active_and_SA45": funnel_SA45,
    "n_active_and_SA60": funnel_SA60,
}
print(f"  活性达标(pIC50>={CFG['admet']['activity_min']}) {funnel_noSA} 个 "
      f"→ 加 SA<=4.5 后 {funnel_SA45} 个 → 加 SA<=6.0 后 {funnel_SA60} 个")

# ============================================================ Part D docking vs 真实活性
banner("D · Docking 分数，和真实活性到底有没有关系？")

dock = pd.read_csv(os.path.join(PROC, "docking.csv"))
truth = dict(zip(df["smiles"], y))
dock["true_pIC50"] = dock["smiles"].map(truth)
sub = dock.dropna(subset=["true_pIC50"])
r_dp, p_dp = pearsonr(sub["pred_pIC50"], sub["docking_score"])
r_dt, p_dt = pearsonr(sub["docking_score"], sub["true_pIC50"])
rho_dt, _ = spearmanr(sub["docking_score"], sub["true_pIC50"])
r_pt, _ = pearsonr(sub["pred_pIC50"], sub["true_pIC50"])
rho_pt, _ = spearmanr(sub["pred_pIC50"], sub["true_pIC50"])

best_by_dock = sub.loc[sub["docking_score"].idxmin()]
best_by_true = sub.loc[sub["true_pIC50"].idxmax()]
best_rank_of_trueTop = int((sub["docking_score"] < best_by_true["docking_score"]).sum()) + 1

R["docking"] = {
    "n": int(len(sub)),
    "pearson_pred_vs_docking": round(float(r_dp), 3),
    "pearson_docking_vs_true": round(float(r_dt), 3),
    "spearman_docking_vs_true": round(float(rho_dt), 3),
    "pearson_pred_vs_true": round(float(r_pt), 3),
    "spearman_pred_vs_true": round(float(rho_pt), 3),
    "docking_score_range": [round(float(sub["docking_score"].min()), 2),
                            round(float(sub["docking_score"].max()), 2)],
    "true_pIC50_range": [round(float(sub["true_pIC50"].min()), 2),
                         round(float(sub["true_pIC50"].max()), 2)],
    "best_by_dock": {"id": best_by_dock["molecule_id"],
                     "docking": round(float(best_by_dock["docking_score"]), 2),
                     "true_pIC50": round(float(best_by_dock["true_pIC50"]), 2)},
    "best_by_true": {"id": best_by_true["molecule_id"],
                     "docking": round(float(best_by_true["docking_score"]), 2),
                     "true_pIC50": round(float(best_by_true["true_pIC50"]), 2),
                     "docking_rank": best_rank_of_trueTop},
}
print(f"  {len(sub)} 个分子：")
print(f"    预测活性 vs 对接分数      Pearson r = {r_dp:+.3f}")
print(f"    对接分数 vs 真实活性      Pearson r = {r_dt:+.3f}   Spearman = {rho_dt:+.3f}")
print(f"    预测活性 vs 真实活性      Pearson r = {r_pt:+.3f}   Spearman = {rho_pt:+.3f}")
print(f"    按对接分排第 1：{best_by_dock['molecule_id']} score={best_by_dock['docking_score']:.2f} "
      f"真实 pIC50={best_by_dock['true_pIC50']:.2f}")
print(f"    真实活性第 1：{best_by_true['molecule_id']} 在对接榜排第 {best_rank_of_trueTop} 名")

with open(os.path.join(OUT_DIR, "mvp_trap_results.json"), "w", encoding="utf-8") as f:
    json.dump(R, f, indent=2, ensure_ascii=False)
log(f"结果写入 {OUT_DIR}/mvp_trap_results.json")

# ============================================================ 图
banner("E · 出图")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Hiragino Sans GB"
plt.rcParams["axes.unicode_minus"] = False
GREEN, RED, BLUE, GREY = "#0FA87E", "#C0392B", "#2B5F8F", "#9A9A9A"

# --- 图 A：两个盒子
fig, ax = plt.subplots(figsize=(9.2, 5.4))
ax.scatter(prot_xyz[:, 0], prot_xyz[:, 1], s=3, c="#dfe4e8", label="1M17 蛋白原子", rasterized=True)
pocket_pts = np.array([a for at in resatoms.values() for a in at])
ax.scatter(pocket_pts[:, 0], pocket_pts[:, 1], s=14, c=GREEN, label="真实口袋残基（20 个）")
ax.scatter(aq4_xyz[:, 0], aq4_xyz[:, 1], s=46, c=RED, marker="*", label="共晶配体 AQ4（厄洛替尼）")
for c, lab, col, ls in [(true_center, "按配体算出的盒子", GREEN, "-"),
                        (hard_center, "伪代码写死的盒子", RED, "--")]:
    ax.add_patch(plt.Rectangle((c[0] - half, c[1] - half), BOX, BOX,
                               fill=False, edgecolor=col, lw=2, ls=ls))
    ax.plot(c[0], c[1], "o" if ls == "-" else "s", ms=7, color=col)
ax.set_xlabel("x (埃)"); ax.set_ylabel("y (埃)")
ax.set_title(f"两个相差 {d_center:.1f} 埃的对接口：正确盒子包住 {cov['true']}/20 个口袋残基，写死的包住 {cov['hard']} 个",
             fontsize=12)
ax.legend(fontsize=10, loc="upper left")
ax.grid(alpha=.25, ls=":")
fig.tight_layout(); fig.savefig(f"{FIG_DIR}/t1_box.png", dpi=170); plt.close(fig)
log("t1_box.png")

# --- 图 B：demo 库 vs 真抑制剂
fig, axes = plt.subplots(2, 1, figsize=(9.2, 6.4))
labels = [r["name"].replace("（伪代码 demo）", "").replace("（ChEMBL API (pref_name=ERLOTINIB)）", "") for r in rows]
cols = [GREY if r["group"] == "demo" else GREEN for r in rows]
vals = [r["pred_pIC50"] for r in rows]
axes[0].barh(range(len(rows))[::-1], vals, color=cols, height=.62)
axes[0].set_yticks(range(len(rows))[::-1]); axes[0].set_yticklabels(labels, fontsize=10)
axes[0].axvline(7, color=BLUE, ls="--", lw=1.4)
axes[0].text(7.05, len(rows) - .4, "pIC50 = 7（IC50 = 100 nM）", fontsize=9, color=BLUE)
for i, v in enumerate(vals):
    axes[0].text(v + .06, len(rows) - 1 - i, f"{v}", va="center", fontsize=10, color="#333")
axes[0].set_xlim(0, max(vals) + 1.2)
axes[0].set_xlabel("模型预测的 pIC50"); axes[0].set_title(
    "同一条 pipeline 跑出来的「候选」：灰色三个是伪代码的 demo 库，绿色是真实抑制剂", fontsize=12)
axes[0].grid(axis="x", alpha=.25, ls=":")

# 下：MW/QED 说明了为什么它们能"通过过滤"
xs = np.arange(len(demo))
w = .38
demo_rows = [r for r in rows if r["group"] == "demo"]
axes[1].bar(xs - w / 2, [r["MW"] for r in demo_rows], w, color="#cfd8dc", label="分子量 MW（<500 ✓）")
axes[1].set_xticks(xs); axes[1].set_xticklabels([r["name"].replace("（伪代码 demo）", "") for r in demo_rows], fontsize=10)
axes[1].axhline(500, color=RED, ls="--", lw=1.2)
axes[1].text(-.4, 512, "MW=500 过滤线", fontsize=9, color=RED)
axes[1].set_ylabel("分子量")
ax2 = axes[1].twinx()
ax2.bar(xs + w / 2, [r["QED"] for r in demo_rows], w, color=BLUE, alpha=.75, label="QED（>0.5 ✓）")
ax2.set_ylabel("QED"); ax2.set_ylim(0, 1)
axes[1].set_title("它们全部通过 MolWt<500 与 QED>0.5 —— 于是被输出成「Lead」", fontsize=12)
fig.tight_layout(); fig.savefig(f"{FIG_DIR}/t2_demo.png", dpi=170); plt.close(fig)
log("t2_demo.png")

# --- 图 C：SAScore + docking
fig, axes = plt.subplots(2, 1, figsize=(9.2, 6.8))
axes[0].hist(sa_v, bins=60, color="#b8c4cc", edgecolor="white")
axes[0].axvline(4.5, color=RED, ls="--", lw=1.6)
axes[0].axvline(6.0, color=BLUE, ls="--", lw=1.6)
ymax = axes[0].get_ylim()[1]
axes[0].text(4.6, ymax * .88, f"生成模型常用 SA≤4.5\n通过率 {(sa_v<=4.5).mean():.1%}", fontsize=9, color=RED)
axes[0].text(6.1, ymax * .62, f"SA≤6.0\n通过率 {(sa_v<=6.0).mean():.1%}", fontsize=9, color=BLUE)
axes[0].set_xlabel("SAScore（越大越难合成）"); axes[0].set_ylabel("分子数")
axes[0].set_title(f"10468 个分子的可合成性分布：强活性组平均 SA {strong.mean():.2f}，弱活性组 {weak.mean():.2f}",
                  fontsize=12)

x = sub["docking_score"]; yy = sub["true_pIC50"]
sc = axes[1].scatter(x, yy, c=sub["pred_pIC50"], cmap="viridis", s=58, edgecolor="white")
cb = fig.colorbar(sc, ax=axes[1]); cb.set_label("模型预测 pIC50", fontsize=9)
axes[1].set_xlabel("Vina 对接分数 (kcal/mol，越负越好)"); axes[1].set_ylabel("真实 pIC50")
axes[1].set_title(f"30 个真实分子：对接分数 vs 真实活性 Pearson r = {r_dt:+.3f}（几乎正交）", fontsize=12)
axes[1].grid(alpha=.25, ls=":")
fig.tight_layout(); fig.savefig(f"{FIG_DIR}/t3_sa_docking.png", dpi=170); plt.close(fig)
log("t3_sa_docking.png")

banner(f"全部完成，用时 {time.time()-T0:.0f}s")
