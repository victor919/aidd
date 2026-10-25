# 文章关联：AI制药实践篇B_实操篇_公众号文章.html
# 脚本职责：配图：六坑地图 A / B
# 使用注解：读取正式工程数据及 mvp_trap_results.json；不重算实验结论。
# 完整对应关系与运行顺序：README.md
"""六坑地图：把实践篇里踩到的六个坑画成两张"导游图"。
所有数字来自 aidd-mvp 工程的实跑产物，不重算结论，只读数据画图。
"""
import os
import json
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

ROOT = os.path.expanduser("~/WorkBuddy/Work/aidd-mvp")
PROC = os.path.join(ROOT, "data", "processed")
RAW = os.path.join(ROOT, "data", "raw")
FIG = os.path.expanduser("~/WorkBuddy/Work/AI制药实操篇_配图")
os.makedirs(FIG, exist_ok=True)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

plt.rcParams["font.family"] = "Hiragino Sans GB"
plt.rcParams["axes.unicode_minus"] = False

GREEN, RED, BLUE = "#0FA87E", "#C0392B", "#2B5F8F"
GREY, AMBER, PURPLE = "#9A9A9A", "#E8A33D", "#7B5CA8"

TRAP = json.load(open(os.path.expanduser("~/WorkBuddy/Work/aidd_practice/mvp_trap_results.json")))
BOX = TRAP["box"]

# ============================================================ 数据来源准备
# --- 1M17 蛋白原子 + 共晶配体 AQ4（手工解析 PDB，不依赖 Biopython）
prot_xyz, aq4_xyz = [], []
with open(os.path.join(RAW, "1M17.pdb"), encoding="utf-8") as f:
    for line in f:
        if line.startswith("ATOM"):
            prot_xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        elif line.startswith("HETATM") and line[17:20].strip() == "AQ4":
            aq4_xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
prot_xyz = np.array(prot_xyz)
aq4_xyz = np.array(aq4_xyz)

true_c = np.array(BOX["true_center"])
hard_c = np.array(BOX["hard_coded_center"])
HALF = 10.0  # 20 埃盒子的一半

# 口袋残基原子（按 5 埃邻域，标准氨基酸残基）
resatoms = {}
AA = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())
with open(os.path.join(RAW, "1M17.pdb"), encoding="utf-8") as f:
    for line in f:
        if line.startswith(("ATOM", "HETATM")):
            res = line[17:20].strip()
            key = (line[21], line[22:27].strip())
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            if line.startswith("ATOM") and res in AA:
                d = np.linalg.norm(xyz - aq4_xyz, axis=1).min()
                if d <= 5.0:
                    resatoms.setdefault(key, []).append(xyz)
            elif res == "AQ4":
                pass
pocket_pts = np.array([a for v in resatoms.values() for a in v])
n_pocket_res = len({k for k in resatoms})

# --- ChEMBL：P00533 下的 16 个靶点条目（2026-10-17 实查 ChEMBL API）
chembl16 = [
    ("CHEMBL203", "SINGLE PROTEIN"),
    ("CHEMBL2111431", "PROTEIN FAMILY"),
    ("CHEMBL2363049", "PROTEIN FAMILY"),
    ("CHEMBL3137284", "CHIMERIC PROTEIN"),
    ("CHEMBL4523680", "PPI"), ("CHEMBL4523747", "PPI"), ("CHEMBL4523998", "PPI"),
    ("CHEMBL4802031", "PPI"), ("CHEMBL5465557", "PPI"), ("CHEMBL6066839", "PPI"),
    ("CHEMBL6066845", "PPI"), ("CHEMBL6193792", "PPI"), ("CHEMBL6193830", "PPI"),
    ("CHEMBL6193837", "PPI"), ("CHEMBL6193841", "PPI"), ("CHEMBL6195769", "PPI"),
]
type_color = {"SINGLE PROTEIN": GREEN, "PROTEIN FAMILY": AMBER,
              "CHIMERIC PROTEIN": BLUE, "PPI": "#c7cdd3"}
type_cn = {"SINGLE PROTEIN": "单蛋白（可用）", "PROTEIN FAMILY": "蛋白家族",
           "CHIMERIC PROTEIN": "嵌合蛋白", "PPI": "蛋白-蛋白相互作用"}

# --- 清洗前活性分布
raw = pd.read_csv(os.path.join(RAW, "raw_activity.csv"))
m = ((raw.standard_units == "nM") & (raw.standard_relation == "=")
     & raw.standard_value.notna() & (raw.standard_value > 0))
p_before = (9 - np.log10(raw.loc[m, "standard_value"])).values
p_before = p_before[np.isfinite(p_before)]

# --- docking vs 真实活性
df = pd.read_csv(os.path.join(PROC, "dataset.csv"))
truth = dict(zip(df["smiles"], df["pIC50_median"]))
dock = pd.read_csv(os.path.join(PROC, "docking.csv"))
dock["true_pIC50"] = dock["smiles"].map(truth)
sub = dock.dropna(subset=["true_pIC50"]).copy()
best_true = sub.loc[sub["true_pIC50"].idxmax()]
best_dock = sub.loc[sub["docking_score"].idxmin()]
rank_of_true_top = int((sub["docking_score"] < best_true["docking_score"]).sum()) + 1

# --- 主动学习曲线
hist = pd.read_csv(os.path.join(PROC, "iteration_history.csv"))

rows = TRAP["demo_library"]["rows"]


def panel(ax, tag, title, subtitle=None):
    """统一的小图抬头：左上角坑号 + 标题"""
    ax.set_title("", fontsize=0)
    ax.text(0.0, 1.14, tag, transform=ax.transAxes, fontsize=13, fontweight="bold",
            color="#fff", va="bottom", ha="left",
            bbox=dict(boxstyle="round,pad=0.32", fc=RED, ec="none"))
    t = ax.text(0.155, 1.145, title, transform=ax.transAxes, fontsize=12.5,
                fontweight="bold", color="#1a1a1a", va="bottom", ha="left")
    if subtitle:
        ax.text(0.0, 1.02, subtitle, transform=ax.transAxes, fontsize=10,
                color="#5a5a5a", va="bottom", ha="left")
    return t


# ==================================================================== 坑地图 A（坑 1-3）
figA, axes = plt.subplots(3, 1, figsize=(9.2, 11.6))

# ---------- 坑 1：写死的口袋坐标
ax = axes[0]
panel(ax, "坑 1", "写死的口袋坐标，落在蛋白外面 31.4 埃",
      "1M17 坐标系 XY 投影：绿框按共晶配体算出，红虚框是伪代码里的常量")
ax.scatter(prot_xyz[:, 0], prot_xyz[:, 1], s=4, c="#dde3e8", label="1M17 蛋白原子", rasterized=True)
ax.scatter(pocket_pts[:, 0], pocket_pts[:, 1], s=16, c=GREEN, label=f"真实口袋残基原子（{n_pocket_res} 个残基）")
ax.scatter(aq4_xyz[:, 0], aq4_xyz[:, 1], s=70, c=RED, marker="*",
           label="共晶配体 AQ4（厄洛替尼）", zorder=5)
for c, ls in [(true_c, "-"), (hard_c, "--")]:
    ax.plot([c[0] - HALF, c[0] + HALF, c[0] + HALF, c[0] - HALF, c[0] - HALF],
            [c[1] - HALF, c[1] - HALF, c[1] + HALF, c[1] + HALF, c[1] - HALF],
            color=GREEN if ls == "-" else RED, lw=2.1, ls=ls, zorder=4)
ax.annotate("", xy=hard_c[:2], xytext=true_c[:2],
            arrowprops=dict(arrowstyle="-|>", color="#333", lw=1.6, ls=":"))
mx, my = (true_c[:2] + hard_c[:2]) / 2
ax.text(mx - 1.5, my + 7.5, f"两个中心相距 {BOX['center_distance']} 埃", fontsize=11, color="#333",
        ha="center", bbox=dict(boxstyle="round,pad=0.25", fc="#fff", ec="#bbb"))
ax.text(hard_c[0], hard_c[1] - HALF - 2.0, "红框里：0 个蛋白原子、0 个口袋残基\nVina 照样会给你一组负分",
        fontsize=10, color=RED, ha="center", va="top")
ax.text(true_c[0], true_c[1] + HALF + 0.8, "绿框里：295 个蛋白原子、19/20 个残基",
        fontsize=10, color=GREEN, ha="center", va="bottom")
ax.set_xlabel("x (埃)", fontsize=10.5)
ax.set_ylabel("y (埃)", fontsize=10.5)
ax.legend(fontsize=9.5, loc="lower left", framealpha=.92)
ax.grid(alpha=.22, ls=":")
ax.set_xlim(-25, 52)
ax.set_ylim(-31, 32)

# ---------- 坑 2：demo 分子被收，真药被拒
ax = axes[1]
panel(ax, "坑 2", "三个解热镇痛药被输出成 Lead，真药厄洛替尼被 QED 拒之门外",
      "横轴 QED（药性），纵轴模型预测 pIC50（活性）——两件不同的事被塞进了同一道门槛")
demo = [r for r in rows if r["group"] == "demo"]
real = [r for r in rows if r["group"] == "real"]
ax.axvline(0.5, color="#333", lw=1.4, ls="--")
ax.axhline(7.0, color="#333", lw=1.4, ls="--")
ax.axhspan(7.0, 11.2, alpha=.05, color=GREEN)
ax.axvspan(0.5, 0.9, alpha=.05, color=GREEN)
ax.plot([r["QED"] for r in demo], [r["pred_pIC50"] for r in demo], "o", ms=13,
        color=GREY, mec="#666", label="伪代码 demo 分子（阿司匹林 / 咖啡因 / 布洛芬）")
ax.plot([real[0]["QED"]], [real[0]["pred_pIC50"]], "*", ms=23, color=RED,
        label="厄洛替尼（真正的 EGFR 药物）")
ax.plot([r["QED"] for r in real[1:]], [r["pred_pIC50"] for r in real[1:]], "D", ms=9,
        color=GREEN, label="数据集中真实的高活性分子")
for r, lab, tx, ty, ha in [(demo[0], "阿司匹林", 0.415, 4.78, "right"),
                           (demo[1], "咖啡因", 0.560, 5.50, "left"),
                           (demo[2], "布洛芬", 0.805, 4.62, "right"),
                           (real[0], "厄洛替尼\n被 QED>0.5 拒收", 0.438, 7.05, "left")]:
    ax.annotate(lab, (r["QED"], r["pred_pIC50"]), xytext=(tx, ty), fontsize=9.5, ha=ha,
                color=RED if lab.startswith("厄洛") else "#555",
                arrowprops=dict(arrowstyle="-", color="#aaa", lw=0.8))
ax.text(0.71, 5.75, "右下：三道「活性筛选」全过\n→ 被 pipeline 当 Lead 输出\n（实际 pIC50 只有 4.9~5.3，等于没活性）",
        fontsize=10, color="#444", ha="center", va="center",
        bbox=dict(boxstyle="round,pad=0.4", fc="#fff9e8", ec="#e0d5b8"))
ax.set_xlabel("QED（越大越像药）", fontsize=10.5)
ax.set_ylabel("模型预测 pIC50", fontsize=10.5)
ax.set_xlim(0.32, 0.9)
ax.set_ylim(4.3, 11.2)
ax.legend(fontsize=9.5, loc="upper left", framealpha=.92)
ax.grid(alpha=.22, ls=":")

# ---------- 坑 3：ChEMBL 16 选 1
ax = axes[2]
panel(ax, "坑 3", "ChEMBL 里搜「EGFR」返回 16 个条目，只有 1 个能直接用",
      "按 UniProt P00533 检索的实跑结果：同名的「EGFR」，语义完全不同")
ax.set_xlim(-0.15, 8.6); ax.set_ylim(-0.55, 4.75); ax.set_aspect("equal")
ax.axis("off")
from matplotlib.patches import Rectangle
counts = {}
for i, (tid, ttype) in enumerate(chembl16):
    r_, c_ = 3 - i // 4, i % 4
    counts[ttype] = counts.get(ttype, 0) + 1
    col = type_color[ttype]
    ax.add_patch(Rectangle((c_ + 0.05, r_ + 0.05), 0.9, 0.9, fc=col,
                           ec=GREEN if ttype == "SINGLE PROTEIN" else "#98a1a8",
                           lw=2.6 if ttype == "SINGLE PROTEIN" else 1.0, zorder=3))
    if ttype == "SINGLE PROTEIN":
        ax.text(c_ + 0.5, r_ + 0.5, "CHEMBL\n203", fontsize=8.0, color="#fff", fontweight="bold",
                ha="center", va="center", zorder=4, linespacing=1.1)
ax.annotate("唯一能用的这一个（人源、单蛋白、EGFR）", xy=(0.5, 3.98), xytext=(-0.05, 4.52),
            fontsize=10.5, color=GREEN, fontweight="bold", va="center", ha="left",
            arrowprops=dict(arrowstyle="-|>", color=GREEN, lw=1.8))
# 右侧图例
lx = 4.9
ax.text(lx, 3.95, "16 个「EGFR」，各自代表什么", fontsize=11.5, fontweight="bold", color="#333")
for j, (k, cn) in enumerate(type_cn.items()):
    y_ = 3.35 - j * 0.62
    ax.add_patch(Rectangle((lx, y_ - 0.19), 0.42, 0.38, fc=type_color[k], ec="#98a1a8", lw=1))
    ax.text(lx + 0.62, y_, f"{cn}：{counts.get(k, 0)} 个", fontsize=10.5, va="center", color="#333")
ax.text(lx, 0.28, "家族靶点的 IC50 混着好几个蛋白的测量；\nPPI 条目是分子胶 / PROTAC 筛选专用。",
        fontsize=10, color="#666", va="center")

figA.tight_layout(h_pad=3.4)
figA.savefig(f"{FIG}/trap_map_A.png", dpi=170, bbox_inches="tight", facecolor="#fff")
plt.close(figA)
print("trap_map_A.png")

# ==================================================================== 坑地图 B（坑 4-6）
figB, axes = plt.subplots(3, 1, figsize=(9.2, 11.6))

# ---------- 坑 4：pIC50 = 17.3
ax = axes[0]
panel(ax, "坑 4", "脏数据里藏着 pIC50 = 17.3 的「超级分子」",
      "清洗前 18381 条记录的 pIC50 分布：右边那个孤零零的点，换算成 IC50 是 5e-18 摩尔每升")
bins = np.arange(2, 18.6, 0.35)
ax.hist(p_before, bins=bins, color="#b9c4cc", edgecolor="none")
ax.set_yscale("log")
top = p_before.max()
ax.plot([top], [1.2], "v", ms=15, color=RED)
ax.annotate(f"pIC50 = {top:.2f}\n一升溶液里摊不到一个分子\n（大概率是单位录错）",
            xy=(top, 1.6), xytext=(10.2, 45), fontsize=10.5, color=RED, fontweight="bold",
            arrowprops=dict(arrowstyle="-|>", color=RED, lw=1.7))
ax.axvspan(11, 18.5, alpha=.10, color=RED)
ax.text(14.7, 500, "pIC50 > 11 的脏数据：16 条\n训练指标几乎不掉，\n但排序 Top 1 必是它",
        fontsize=10.5, color=RED, va="center", ha="center",
        bbox=dict(boxstyle="round,pad=0.35", fc="#fff", ec="#e8b4b4"))
ax.axvline(9 - np.log10(1e-9), color=GREEN, lw=1.6, ls="--")
ax.text(9.15, 20, "IC50 = 1 nM", fontsize=10, color=GREEN)
ax.set_xlabel("pIC50（清洗前，原始记录换算）", fontsize=10.5)
ax.set_ylabel("记录数（对数轴）", fontsize=10.5)
ax.grid(axis="y", alpha=.22, ls=":")
ax.set_axisbelow(True)

# ---------- 坑 5：docking 与真实活性正交
ax = axes[1]
panel(ax, "坑 5", "对接分数与真实活性几乎正交：30 个分子，Pearson r = +0.10",
      "横轴 AutoDock Vina 对接分（越负越好），纵轴真实测得的 pIC50")
ax.plot(sub["docking_score"], sub["true_pIC50"], "o", ms=9, color="#7f9fb5", mec="#4a7a99")
ax.plot([best_dock["docking_score"]], [best_dock["true_pIC50"]], "s", ms=15, color=BLUE,
        mec="#123", label="对接榜第 1 名（它的真实活性排第 18）")
ax.plot([best_true["docking_score"]], [best_true["true_pIC50"]], "*", ms=25, color=RED, mec="#600",
        label=f"真实活性第 1 名（对接榜只排第 {rank_of_true_top}）")
ax.annotate("这里的每一个点，都是几十万次构象搜索换来的",
            xy=(-9.1, 6.9), xytext=(-9.53, 6.15), fontsize=10, color=BLUE)
ax.annotate("活性最强的分子：对接第 19 名", xy=(best_true["docking_score"], best_true["true_pIC50"]),
            xytext=(-8.68, 10.85), fontsize=10.5, color=RED, fontweight="bold",
            arrowprops=dict(arrowstyle="-|>", color=RED, lw=1.6))
ax.set_xlabel("Vina 对接分数 (kcal/mol)", fontsize=10.5)
ax.set_ylabel("真实 pIC50（实验测得）", fontsize=10.5)
from scipy.stats import pearsonr
r, _ = pearsonr(sub["docking_score"], sub["true_pIC50"])
ax.text(0.985, 0.055, f"Pearson r = {r:+.2f}（若真相关应为明显的负斜率）",
        transform=ax.transAxes, fontsize=11, ha="right", color="#333", fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.35", fc="#fff9e8", ec="#e0d5b8"))
ax.legend(fontsize=10, loc="upper left", framealpha=.93)
ax.grid(alpha=.22, ls=":")
ax.set_ylim(5.3, 11.3)

# ---------- 坑 6：greedy 原地踏步
ax = axes[2]
panel(ax, "坑 6", "只挑模型预测最强的分子去做实验，六轮下来几乎原地踏步",
      "三条曲线是同一份预算（每轮 300 个分子、共 6 轮）下的测试集 MAE，越低越好")
style = {"uncertainty": (GREEN, "o", "-", "uncertainty：专挑模型最拿不准的"),
         "random": (BLUE, "s", "--", "random：随机挑"),
         "greedy": (RED, "^", "-", "greedy：只挑预测活性最高的")}
for k, (c, mk, ls, lab) in style.items():
    g = hist[hist.strategy == k]
    ax.plot(g.n_labeled, g.MAE, marker=mk, ls=ls, lw=2.3, ms=8, color=c, label=lab)
g = hist[hist.strategy == "greedy"]
ax.annotate(f"0.90 → {g.MAE.iloc[-1]:.2f}\n六轮预算，MAE 只降了 0.02",
            xy=(g.n_labeled.iloc[-1], g.MAE.iloc[-1]), xytext=(1030, 0.862),
            fontsize=10.5, color=RED, fontweight="bold",
            arrowprops=dict(arrowstyle="-|>", color=RED, lw=1.6))
gu = hist[hist.strategy == "uncertainty"]
ax.annotate("前四轮落后随机，\n1800 样本才反超，拿到最优 0.747",
            xy=(gu.n_labeled.iloc[-1], gu.MAE.iloc[-1]), xytext=(1080, 0.756),
            fontsize=10.5, color=GREEN,
            arrowprops=dict(arrowstyle="-|>", color=GREEN, lw=1.6))
ax.set_xlabel("累计已标记的分子数（实验预算）", fontsize=10.5)
ax.set_ylabel("测试集 MAE（对数单位）", fontsize=10.5)
ax.legend(fontsize=10, loc="lower left", framealpha=.93)
ax.grid(alpha=.22, ls=":")
ax.set_axisbelow(True)

figB.tight_layout(h_pad=3.4)
figB.savefig(f"{FIG}/trap_map_B.png", dpi=170, bbox_inches="tight", facecolor="#fff")
plt.close(figB)
print("trap_map_B.png")

print(f"坑1 距离 {BOX['center_distance']} 埃 / 盒内蛋白原子 {BOX['protein_atoms_in_hard_box']}")
print(f"坑2 demo 通过 {TRAP['demo_library']['demo_pass_QED05']} / 厄洛替尼 QED {real[0]['QED']}")
print(f"坑3 ChEMBL 条目 {len(chembl16)} 类型分布 {counts}")
print(f"坑4 清洗前最大 pIC50 {p_before.max():.2f} / >11 的 {int((p_before>11).sum())} 条")
print(f"坑5 n={len(sub)} r={r:+.3f} 真实活性第一在对接榜第 {rank_of_true_top} 名")
print("坑6 greedy MAE:", round(float(g.MAE.iloc[0]), 3), "->", round(float(g.MAE.iloc[-1]), 3))
