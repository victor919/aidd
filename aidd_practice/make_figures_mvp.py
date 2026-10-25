# 文章关联：AI制药实践篇B_实操篇_公众号文章.html
# 脚本职责：早期配图：m1–m6
# 使用注解：读取 data/egfr_mvp/；最终 m2 / m4 / m5 / m6 由 redraw_m_series.py 按正式工程口径重画。
# 完整对应关系与运行顺序：README.md
"""实操篇配图：全部基于本机实跑数据"""
import os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib import font_manager

for f in ["Hiragino Sans GB", "PingFang SC", "Heiti SC"]:
    try:
        font_manager.findfont(f, fallback_to_default=False)
        plt.rcParams["font.sans-serif"] = [f]
        print("font:", f)
        break
    except Exception:
        continue
plt.rcParams["axes.unicode_minus"] = False

BASE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(BASE, "data", "egfr_mvp")
OUT = "/Users/dujiepeng/WorkBuddy/Work/AI制药实操篇_配图"
os.makedirs(OUT, exist_ok=True)

GREEN, BLUE, RED, GREY, AMBER = "#0FA87E", "#2B5F8F", "#C0392B", "#9A9A9A", "#F5A623"
INK = "#1A1A1A"

# ============================================================
# m1: MVP 脚本地图
# ============================================================
fig, ax = plt.subplots(figsize=(9.4, 6.4))
ax.set_xlim(0, 11); ax.set_ylim(0, 7.4); ax.axis("off")

def box(x, y, w, h, text, fc, tc="#1A1A1A", fs=9.5, bold=False, lw=1.2):
    ax.add_patch(FancyBboxPatch((x - w/2, y - h/2), w, h, boxstyle="round,pad=0.05,rounding_size=0.12",
                                fc=fc, ec="#FFFFFF" if fc != "#FFFFFF" else GREY, lw=lw))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, color=tc,
            fontweight="bold" if bold else "normal", linespacing=1.5, zorder=3)

def arrow(x1, y1, x2, y2, color=GREY, style="-|>", lw=1.4, rad=0.0):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style, mutation_scale=13,
                                 color=color, lw=lw, shrinkA=2, shrinkB=2,
                                 connectionstyle=f"arc3,rad={rad}", zorder=2))

ax.text(5.5, 7.05, "EGFR MVP：一张图跑完的端到端流水线", ha="center", fontsize=14,
        fontweight="bold", color=INK)
ax.text(5.5, 6.72, "本机实跑：MacBook Pro M1 Pro / conda aidd / Python 3.11.14 · RDKit 2026.03.5",
        ha="center", fontsize=9, color=GREY)

steps = [
    (1.10, "① target_rank",   "target_rank.csv",   "机制与证据", AMBER),
    (2.82, "② download",      "1M17.pdb",          "实验结构",   BLUE),
    (4.54, "③ pocket",        "pocket.json",       "口袋残基",   BLUE),
    (6.26, "④ fetch ChEMBL",  "raw_activity.csv",  "真实活性",   GREEN),
    (7.98, "⑤ standardize",   "dataset.csv",       "标准化",     GREEN),
    (9.70, "⑥ scaffold split","train / test",      "切分",       GREEN),
]
for i, (x, script, out, tag, c) in enumerate(steps):
    y = 5.55
    ax.add_patch(FancyBboxPatch((x-0.80, y-0.31), 1.60, 0.62,
                 boxstyle="round,pad=0.05,rounding_size=0.12", fc="#FFFFFF", ec=c, lw=1.6, zorder=1))
    ax.text(x, y, script, ha="center", va="center", fontsize=9.2, color=INK, zorder=3)
    box(x, 4.62, 1.34, 0.42, out, "#F2F4F6", fs=8.2)
    ax.text(x, 4.28, tag, ha="center", fontsize=8, color=c, fontweight="bold")
    arrow(x, 5.20, x, 4.86, c)
    if i < len(steps) - 1:
        arrow(x + 0.82, 5.55, steps[i+1][0] - 0.82, 5.55, GREY, rad=0.0)

# 分叉：RF / GCN
box(7.98, 3.42, 1.50, 0.58, "⑦ rf_baseline\nECFP4 + RF", "#FFFFFF", fs=8.4)
ax.add_patch(FancyBboxPatch((7.23, 3.13), 1.50, 0.58, boxstyle="round,pad=0.05,rounding_size=0.12",
                            fc="#FFFFFF", ec=GREEN, lw=1.6))
box(9.70, 3.42, 1.50, 0.58, "⑧ gcn\n分子图卷积", "#FFFFFF", fs=8.4)
ax.add_patch(FancyBboxPatch((8.95, 3.13), 1.50, 0.58, boxstyle="round,pad=0.05,rounding_size=0.12",
                            fc="#FFFFFF", ec=GREEN, lw=1.6))
arrow(7.98, 4.26, 7.98, 3.75, GREEN); arrow(9.70, 4.26, 9.70, 3.75, GREEN)
box(8.84, 2.55, 2.60, 0.44, "model_prediction.csv（MAE 对照）", "#F2F4F6", fs=8.4)
arrow(7.98, 3.11, 8.60, 2.80, GREEN); arrow(9.70, 3.11, 9.10, 2.80, GREEN)

box(8.84, 1.72, 1.70, 0.44, "⑨ docking", "#FFFFFF", fs=8.6)
ax.add_patch(FancyBboxPatch((7.99, 1.50), 1.70, 0.44, boxstyle="round,pad=0.05,rounding_size=0.12",
                            fc="#FFFFFF", ec=GREY, lw=1.2, ls=(0, (3, 2))))
ax.text(8.84, 1.22, "需单独环境（Vina / Meeko）", ha="center", fontsize=7.5, color=GREY)
arrow(8.84, 2.31, 8.84, 1.97, GREY)

box(8.84, 0.72, 2.60, 0.50, "⑩ admet → candidate_ranking.csv", "#FFFFFF", fs=8.4)
ax.add_patch(FancyBboxPatch((7.54, 0.47), 2.60, 0.50, boxstyle="round,pad=0.05,rounding_size=0.12",
                            fc="#FFFFFF", ec=AMBER, lw=1.6))
arrow(8.84, 1.18, 8.84, 1.00, AMBER)

box(4.60, 0.72, 2.30, 0.50, "⑪ 实验接口 → assay_results.csv", "#FFFFFF", fs=8.4)
ax.add_patch(FancyBboxPatch((3.45, 0.47), 2.30, 0.50, boxstyle="round,pad=0.05,rounding_size=0.12",
                            fc="#FFFFFF", ec=RED, lw=1.6))
box(1.60, 0.72, 1.80, 0.50, "⑫ closed_loop\n→ Dataset v2", "#FFFFFF", fs=8.4)
ax.add_patch(FancyBboxPatch((0.70, 0.47), 1.80, 0.50, boxstyle="round,pad=0.05,rounding_size=0.12",
                            fc="#FFFFFF", ec=GREEN, lw=1.6))
arrow(3.42, 0.72, 2.53, 0.72, RED)
arrow(7.51, 0.72, 5.78, 0.72, AMBER)
ax.add_patch(FancyArrowPatch((1.60, 1.00), (1.60, 3.05), arrowstyle="-|>", mutation_scale=13,
                             color=GREEN, lw=1.6, ls=(0, (5, 3)),
                             connectionstyle="arc3,rad=-0.55", zorder=2))
ax.text(0.42, 2.05, "新标签\n回流重训", ha="center", fontsize=8.2, color=GREEN, fontweight="bold")

plt.tight_layout(rect=[0, 0, 1, 0.965])
plt.savefig(f"{OUT}/m1_mvp_map.png", dpi=170, bbox_inches="tight", facecolor="white")
plt.close(); print("m1 done")

# ============================================================
# m2: random vs scaffold split
# ============================================================
res = {
    "ECFP4+RF\nrandom split":  (0.524, 0.720, 0.716, 0.851),
    "ECFP4+RF\nscaffold split": (0.657, 0.883, 0.633, 0.797),
}
fig, axes = plt.subplots(3, 1, figsize=(8.6, 8.0))
labels = list(res.keys()); colors = [GREEN, RED]

axes[0].bar(labels, [res[k][0] for k in labels], color=colors, width=0.55)
for i, k in enumerate(labels):
    axes[0].text(i, res[k][0] + 0.012, f"{res[k][0]:.3f}", ha="center", fontsize=10.5, fontweight="bold", color=INK)
axes[0].set_title("MAE（越低越好）", fontsize=11, color=INK)
axes[0].set_ylim(0, 0.80); axes[0].grid(axis="y", ls=":", alpha=0.4); axes[0].set_axisbelow(True)
axes[0].tick_params(labelsize=8.5)
axes[0].annotate("", xy=(1, 0.657), xytext=(1, 0.524),
                 arrowprops=dict(arrowstyle="-|>", color=RED, lw=1.6))
axes[0].text(1.02, 0.59, "+25%", fontsize=9.5, color=RED, fontweight="bold")

axes[1].bar(labels, [res[k][2] for k in labels], color=colors, width=0.55)
for i, k in enumerate(labels):
    axes[1].text(i, res[k][2] + 0.015, f"{res[k][2]:.3f}", ha="center", fontsize=10.5, fontweight="bold", color=INK)
axes[1].set_title("R2（越高越好）", fontsize=11, color=INK)
axes[1].set_ylim(0, 0.85); axes[1].grid(axis="y", ls=":", alpha=0.4); axes[1].set_axisbelow(True)
axes[1].tick_params(labelsize=8.5)

# 骨架重叠示意
ov = {"random split": 652, "scaffold split": 0}
b = axes[2].bar(list(ov.keys()), list(ov.values()), color=[GREEN, RED], width=0.55)
for i, (k, v) in enumerate(ov.items()):
    axes[2].text(i, v + 18, str(v), ha="center", fontsize=10.5, fontweight="bold", color=INK)
axes[2].set_title("测试集骨架在训练集出现过的数量", fontsize=10.5, color=INK)
axes[2].set_ylim(0, 800); axes[2].grid(axis="y", ls=":", alpha=0.4); axes[2].set_axisbelow(True)
axes[2].tick_params(labelsize=8.5)

fig.suptitle("同一个模型、同一批分子，只换切分方式：指标掉了一大截（本机实跑，n=10 463）",
             fontsize=12.5, fontweight="bold", color=INK, y=0.985)
plt.tight_layout(rect=[0, 0, 1, 0.965])
plt.savefig(f"{OUT}/m2_split.png", dpi=170, bbox_inches="tight", facecolor="white")
plt.close(); print("m2 done")

# ============================================================
# m3: 1M17 口袋与编号偏移
# ============================================================
pocket = json.load(open(os.path.join(DATA, "pocket.json")))
res5 = pocket["residues"]
OFFSET = 24
fig, axes = plt.subplots(2, 1, figsize=(8.6, 7.0))

ax = axes[0]
# 重新计算 5埃 接触次数（与 step1 同逻辑）
import sys
sys.path.insert(0, os.path.join(BASE, "vendor"))
from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB.Polypeptide import is_aa
parser = PDBParser(QUIET=True)
st = parser.get_structure("E", os.path.join(DATA, "1M17.pdb"))
lig_atoms, prot_atoms = [], []
for ch in st[0]:
    for res in ch:
        if res.get_resname().strip() == "AQ4":
            lig_atoms.extend(list(res.get_atoms()))
        elif is_aa(res, standard=True):
            prot_atoms.extend(list(res.get_atoms()))
ns = NeighborSearch(prot_atoms)
cnt = {}
for a in lig_atoms:
    for r in ns.search(a.coord, 5.0, level="R"):
        key = (r.id[1], r.resname)
        cnt[key] = cnt.get(key, 0) + 1
items = sorted(cnt.items(), key=lambda kv: -kv[1])[:20][::-1]
labels3 = [f"{nm} {num}" for (num, nm), _ in items]
vals3 = [c for _, c in items]
ax.barh(np.arange(len(items)), vals3, color=GREEN, alpha=0.85, height=0.62)
ax.set_yticks(np.arange(len(items))); ax.set_yticklabels(labels3, fontsize=8.4, family="monospace")
ax.set_xlabel("与配体原子的近距离接触次数（5 埃 内）", fontsize=9.5, color="#3F3F3F")
ax.set_title("erlotinib 周围 5 埃 的 20 个口袋残基", fontsize=11.5, color=INK, fontweight="bold")
for i, v in enumerate(vals3):
    ax.text(v + 0.3, i, str(v), va="center", fontsize=8, color="#3F3F3F")
ax.grid(axis="x", ls=":", alpha=0.35); ax.set_axisbelow(True)
for s in ["top", "right"]:
    ax.spines[s].set_visible(False)

ax = axes[1]
pairs = [("L694", "L718"), ("G695", "G719"), ("K721", "K745 催化赖氨酸"), ("T766", "T790 gatekeeper"),
         ("M769", "M793 hinge"), ("C773", "C797 共价残基"),
         ("L820", "L844"), ("T830", "T854"), ("D831", "D855 DFG-Asp")]
y = np.arange(len(pairs))[::-1]
ax.barh(y, [1]*len(pairs), color="#FFFFFF", edgecolor="#E0E3E6")
ax.set_yticks(y); ax.set_yticklabels([p[0] for p in pairs], fontsize=9.5, color=BLUE, fontweight="bold")
ax.set_xticks([])
ax2 = ax.twinx()
ax2.set_ylim(ax.get_ylim()); ax2.set_yticks(y)
ax2.set_yticklabels([p[1] for p in pairs], fontsize=9.5, color=RED, fontweight="bold")
ax.set_xlim(0, 1)
for s in ["top", "right", "left", "bottom"]:
    ax.spines[s].set_visible(False); ax2.spines[s].set_visible(False)
ax.set_title("同一位点，两套编号：差整整 24", fontsize=11.5, color=INK, fontweight="bold")
ax.text(0.5, -0.9, "左：1M17 结构编号    →    右：UniProt P00533 / 临床通用编号",
        ha="center", fontsize=9.5, color=GREY)
ax.text(0.5, -1.5, "校验方式：1M17 链 A 的 312 个观测残基与 P00533 序列逐位比对，偏移 24 时 312/312 全中",
        ha="center", fontsize=8.6, color=GREY)
for i, (a_, b_) in enumerate(pairs):
    yy = y[i]
    ax.annotate("", xy=(0.52, yy), xytext=(0.44, yy),
                arrowprops=dict(arrowstyle="-|>", color=GREY, lw=1.2))
plt.tight_layout(rect=[0, 0, 1, 0.965])
plt.savefig(f"{OUT}/m3_pocket.png", dpi=170, bbox_inches="tight", facecolor="white")
plt.close(); print("m3 done")

# ============================================================
# m4: 活性-QED 张力 + 帕累托
# ============================================================
full = pd.read_csv(os.path.join(DATA, "dataset_with_props.csv"))
_agg = pd.read_csv(os.path.join(DATA, "dataset.csv"))[["std_smiles", "pIC50_min", "pIC50_max", "range"]]
full = full.merge(_agg, on="std_smiles", how="left")
fig, axes = plt.subplots(2, 1, figsize=(8.6, 7.2))

ax = axes[0]
sub = full[(full.pIC50_median >= 3) & (full.pIC50_median <= 12)]
ax.scatter(sub.pIC50_median, sub.QED, s=7, c=GREEN, alpha=0.28, edgecolors="none")
strong = sub[sub.pIC50_median >= 8]; weak = sub[sub.pIC50_median < 7]
ax.scatter([strong.pIC50_median.mean()], [strong.QED.mean()], s=150, marker="s", c=RED,
           edgecolors="white", linewidths=1.5, zorder=5, label=f"强活性组 n={len(strong)}")
ax.scatter([weak.pIC50_median.mean()], [weak.QED.mean()], s=150, marker="^", c=BLUE,
           edgecolors="white", linewidths=1.5, zorder=5, label=f"弱活性组 n={len(weak)}")
ax.set_xlabel("pIC50（越大越活）", fontsize=10.5, color="#3F3F3F")
ax.set_ylabel("QED（越接近 1 越像药）", fontsize=10.5, color="#3F3F3F")
ax.set_title("活性越强，分子反而越不像药", fontsize=11.5, color=INK, fontweight="bold")
ax.legend(fontsize=9, loc="lower left", framealpha=0.9)
ax.grid(ls=":", alpha=0.35); ax.set_axisbelow(True)
ax.annotate(f"QED 均值\n强 0.346 / 弱 0.441", xy=(9.6, 0.08), fontsize=9, color=RED,
            ha="center", bbox=dict(boxstyle="round,pad=0.35", fc="#FBF8F1", ec=RED, lw=1))
ax.annotate(f"MW 均值\n强 511.9 / 弱 455.2", xy=(4.35, 0.90), fontsize=9, color=BLUE,
            ha="center", bbox=dict(boxstyle="round,pad=0.35", fc="#F2F4F6", ec=BLUE, lw=1))

ax = axes[1]
s2 = full[(full.pIC50_median >= 7) & (full.pIC50_median <= 11)]
a = s2.pIC50_median.values; b = s2.QED.values
pts = np.stack([a, b], axis=1)
n = len(pts); keep = np.ones(n, bool)
for i in range(n):
    dom = np.all(pts >= pts[i], axis=1) & np.any(pts > pts[i], axis=1)
    if dom.any():
        keep[i] = False
ax.scatter(a[~keep], b[~keep], s=7, c="#D8DEE3", edgecolors="none", label=f"被支配 {int((~keep).sum())}")
ax.scatter(a[keep], b[keep], s=90, c=RED, edgecolors="white", linewidths=1.2, zorder=5,
           label=f"帕累托前沿 {int(keep.sum())}")
pf = pts[keep]; o = np.argsort(pf[:, 0])
ax.plot(pf[o, 0], pf[o, 1], c=RED, lw=1.6, ls="--", zorder=4)
ax.set_xlabel("pIC50", fontsize=10.5, color="#3F3F3F"); ax.set_ylabel("QED", fontsize=10.5, color="#3F3F3F")
ax.set_title("没有唯一的“最佳分子”，只有一条前沿", fontsize=11.5, color=INK, fontweight="bold")
ax.legend(fontsize=9.5, framealpha=0.9); ax.grid(ls=":", alpha=0.35); ax.set_axisbelow(True)
ax.annotate("Top 活性", xy=(11.0, 0.591), xytext=(9.9, 0.42), fontsize=9, color=RED,
            arrowprops=dict(arrowstyle="->", color=RED, lw=1.3))
ax.annotate("Top QED", xy=(7.04, 0.935), xytext=(7.5, 0.80), fontsize=9, color=BLUE,
            arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.3))
plt.tight_layout(rect=[0, 0, 1, 0.965])
plt.savefig(f"{OUT}/m4_qed_pareto.png", dpi=170, bbox_inches="tight", facecolor="white")
plt.close(); print("m4 done")

# ============================================================
# m5: 主动学习三策略
# ============================================================
hist = json.load(open(os.path.join(DATA, "al_history.json")))
fig, ax = plt.subplots(figsize=(8.6, 5.2))
style = {"uncertainty": (GREEN, "o", "挑模型最不确定的（uncertainty）"),
         "random": (GREY, "s", "随机挑（random）"),
         "greedy": (RED, "^", "只挑预测最强的（greedy）")}
for k, (c, m, lab) in style.items():
    xs = [h["n_labeled"] for h in hist[k]]; ys = [h["mae"] for h in hist[k]]
    ax.plot(xs, ys, marker=m, color=c, lw=2.0, ms=7, label=lab, zorder=3)
    ax.annotate(f"{ys[-1]:.3f}", (xs[-1], ys[-1]), textcoords="offset points", xytext=(8, -3),
                fontsize=9.5, color=c, fontweight="bold")
ax.axhline(0.592, color=BLUE, ls="--", lw=1.4, zorder=2)
ax.text(1950, 0.596, "全量数据训练的上限 MAE = 0.592", fontsize=9, color=BLUE)
ax.set_xlabel("已标记分子数（每轮新增 300）", fontsize=10.5, color="#3F3F3F")
ax.set_ylabel("测试集 MAE（越低越好）", fontsize=10.5, color="#3F3F3F")
ax.set_title("同样做 6 轮实验，挑谁去测，结果差很多（本机实跑）",
             fontsize=12, color=INK, fontweight="bold", pad=30)
ax.legend(fontsize=9.5, framealpha=0.95, loc=(0.35, 0.10))
ax.grid(ls=":", alpha=0.35); ax.set_axisbelow(True)
ax.annotate("越挑越差：300 → 600 时，\nMAE 反而从 0.858 升到 0.877",
            xy=(600, 0.877), xytext=(660, 0.945), fontsize=9, color=RED,
            arrowprops=dict(arrowstyle="->", color=RED, lw=1.3),
            bbox=dict(boxstyle="round,pad=0.35", fc="#FFF5F5", ec=RED, lw=1))
ax.set_ylim(0.56, 0.98)
ax.set_xlim(280, 2080)
plt.tight_layout(rect=[0, 0, 1, 0.965])
plt.savefig(f"{OUT}/m5_active_learning.png", dpi=170, bbox_inches="tight", facecolor="white")
plt.close(); print("m5 done")

# ============================================================
# m6: 数据体检
# ============================================================
fig, axes = plt.subplots(2, 1, figsize=(8.6, 6.4))
y_all = full.pIC50_median.values
ax = axes[0]
ax.hist(y_all, bins=60, color=GREEN, alpha=0.75, edgecolor="white", linewidth=0.4)
ax.set_xlabel("pIC50", fontsize=10.5, color="#3F3F3F"); ax.set_ylabel("分子数", fontsize=10.5, color="#3F3F3F")
ax.set_title("pIC50 分布：右边那条尾巴有问题", fontsize=11.5, color=INK, fontweight="bold")
ax.axvspan(11, 18, color=RED, alpha=0.16)
ax.annotate("pIC50 > 11\n= IC50 < 0.01 nM\n物理上不可能", xy=(13.5, 320), fontsize=9, color=RED,
            ha="center", bbox=dict(boxstyle="round,pad=0.35", fc="#FFF5F5", ec=RED, lw=1))
ax.annotate("", xy=(12.0, 8), xytext=(13.5, 300),
            arrowprops=dict(arrowstyle="->", color=RED, lw=1.2))
ax.text(14.6, 60, "24 条脏数据\n藏在这里（太矮看不见）\n但 Top1 排序会毁在这", fontsize=8.5,
        color=RED, ha="center", linespacing=1.5)
ax.grid(axis="y", ls=":", alpha=0.35); ax.set_axisbelow(True)

ax = axes[1]
# 同一分子多次测量的分歧
d = full[full.n_meas >= 2]
ax.hist(d["range"].clip(0, 4), bins=50, color=BLUE, alpha=0.75, edgecolor="white", linewidth=0.4)
ax.axvline(0.71, color=RED, lw=1.6, ls="--")
ax.text(0.78, ax.get_ylim()[1]*0.75, "中位数分歧 0.71\n个对数单位", fontsize=9, color=RED)
ax.set_xlabel("同一分子多次测量，pIC50 的最大-最小跨度", fontsize=10, color="#3F3F3F")
ax.set_ylabel("分子数", fontsize=10.5, color="#3F3F3F")
ax.set_title("35.5% 的分子被测过不止一次，而且测得不一样", fontsize=11.5, color=INK, fontweight="bold")
ax.grid(axis="y", ls=":", alpha=0.35); ax.set_axisbelow(True)
plt.tight_layout(rect=[0, 0, 1, 0.965])
plt.savefig(f"{OUT}/m6_data_qc.png", dpi=170, bbox_inches="tight", facecolor="white")
plt.close(); print("m6 done")
print("all done ->", OUT)
