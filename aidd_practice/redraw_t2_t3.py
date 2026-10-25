# 文章关联：AI制药实践篇B_实操篇_公众号文章.html
# 脚本职责：版式修订：t2 / t3
# 使用注解：依赖 mvp_trap_results.json 及正式工程数据；覆盖 mvp_trap_demo.py 生成的同名图。
# 完整对应关系与运行顺序：README.md
"""重绘 t2 / t3（修文字截顶问题）。数据来自 mvp_trap_results.json"""
import os, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Hiragino Sans GB"
plt.rcParams["axes.unicode_minus"] = False
GREEN, RED, BLUE, GREY = "#0FA87E", "#C0392B", "#2B5F8F", "#9A9A9A"

ROOT = os.path.expanduser("~/WorkBuddy/Work/aidd-mvp")
PROC = os.path.join(ROOT, "data", "processed")
FIG = os.path.expanduser("~/WorkBuddy/Work/AI制药实操篇_配图")
R = json.load(open(os.path.expanduser("~/WorkBuddy/Work/aidd_practice/mvp_trap_results.json")))
rows = R["demo_library"]["rows"]

# ---------- t2 ----------
fig, axes = plt.subplots(2, 1, figsize=(9.2, 6.4))
labels = [r["name"].split("（")[0] for r in rows]
cols = [GREY if r["group"] == "demo" else GREEN for r in rows]
vals = [r["pred_pIC50"] for r in rows]
axes[0].barh(range(len(rows))[::-1], vals, color=cols, height=.62)
axes[0].set_yticks(range(len(rows))[::-1]); axes[0].set_yticklabels(labels, fontsize=10)
axes[0].axvline(7, color=BLUE, ls="--", lw=1.4)
axes[0].text(7.08, 0.06, "pIC50 = 7（IC50 = 100 nM）", fontsize=9, color=BLUE,
             transform=matplotlib.transforms.blended_transform_factory(axes[0].transData, axes[0].transAxes))
for i, v in enumerate(vals):
    axes[0].text(v + .06, len(rows) - 1 - i, f"{v}", va="center", fontsize=10, color="#333")
axes[0].set_xlim(0, max(vals) + 1.2)
axes[0].set_xlabel("模型预测的 pIC50")
axes[0].set_title("同一条 pipeline 跑出来的「候选」：灰色三个是伪代码的 demo 库，绿色是真实抑制剂", fontsize=12)
axes[0].grid(axis="x", alpha=.25, ls=":")

demo_rows = [r for r in rows if r["group"] == "demo"]
xs = np.arange(len(demo_rows)); w = .38
axes[1].bar(xs - w / 2, [r["MW"] for r in demo_rows], w, color="#cfd8dc", label="分子量 MW（<500 ✓）")
axes[1].set_xticks(xs); axes[1].set_xticklabels([r["name"].split("（")[0] for r in demo_rows], fontsize=10)
axes[1].axhline(500, color=RED, ls="--", lw=1.2)
axes[1].text(-.42, 452, "MW = 500 过滤线", fontsize=9, color=RED)
axes[1].set_ylabel("分子量"); axes[1].set_ylim(0, 545)
ax2 = axes[1].twinx()
ax2.bar(xs + w / 2, [r["QED"] for r in demo_rows], w, color=BLUE, alpha=.75, label="QED（>0.5 ✓）")
ax2.set_ylabel("QED"); ax2.set_ylim(0, 1)
axes[1].set_title("它们全部通过 MolWt<500 与 QED>0.5 —— 于是被输出成「Lead」", fontsize=12)
fig.tight_layout(); fig.savefig(f"{FIG}/t2_demo.png", dpi=170); plt.close(fig)
print("t2 ok")

# ---------- t3 ----------
df = pd.read_csv(os.path.join(PROC, "dataset.csv"))
y = df["pIC50_median"].values
import sys
sys.path.append("/Users/dujiepeng/miniconda3/envs/aidd/lib/python3.11/site-packages")
from rdkit import RDConfig
sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer
from rdkit import Chem
sa = np.array([sascorer.calculateScore(Chem.MolFromSmiles(s)) if Chem.MolFromSmiles(s) else np.nan
               for s in df["smiles"]])
v = ~np.isnan(sa); sa_v, y_v = sa[v], y[v]
strong, weak = sa_v[y_v >= 8], sa_v[y_v < 7]

dock = pd.read_csv(os.path.join(PROC, "docking.csv"))
truth = dict(zip(df["smiles"], y))
dock["true_pIC50"] = dock["smiles"].map(truth)
sub = dock.dropna(subset=["true_pIC50"])
r_dt = R["docking"]["pearson_docking_vs_true"]

fig, axes = plt.subplots(2, 1, figsize=(9.2, 6.8))
axes[0].hist(sa_v, bins=60, color="#b8c4cc", edgecolor="white")
axes[0].axvline(4.5, color=RED, ls="--", lw=1.6)
axes[0].axvline(6.0, color=BLUE, ls="--", lw=1.6)
ymax = axes[0].get_ylim()[1]
axes[0].text(4.6, ymax * .76, f"生成模型常用 SA≤4.5\n通过率 {(sa_v<=4.5).mean():.1%}", fontsize=9, color=RED, va="top")
axes[0].text(6.1, ymax * .50, f"SA≤6.0\n通过率 {(sa_v<=6.0).mean():.1%}", fontsize=9, color=BLUE, va="top")
axes[0].set_xlabel("SAScore（越大越难合成）"); axes[0].set_ylabel("分子数")
axes[0].set_title(f"10468 个分子的可合成性分布：强活性组平均 SA {strong.mean():.2f}，弱活性组 {weak.mean():.2f}",
                  fontsize=12)
sc = axes[1].scatter(sub["docking_score"], sub["true_pIC50"], c=sub["pred_pIC50"],
                     cmap="viridis", s=58, edgecolor="white")
cb = fig.colorbar(sc, ax=axes[1]); cb.set_label("模型预测 pIC50", fontsize=9)
axes[1].set_xlabel("Vina 对接分数 (kcal/mol，越负越好)"); axes[1].set_ylabel("真实 pIC50")
axes[1].set_title(f"30 个真实分子：对接分数 vs 真实活性 Pearson r = {r_dt:+.3f}（几乎正交）", fontsize=12)
axes[1].grid(alpha=.25, ls=":")
fig.tight_layout(); fig.savefig(f"{FIG}/t3_sa_docking.png", dpi=170); plt.close(fig)
print("t3 ok")
