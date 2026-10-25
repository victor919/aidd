# 文章关联：AI制药实践篇B_实操篇_公众号文章.html
# 脚本职责：正式口径配图：m2 / m4 / m5 / m6
# 使用注解：读取 ../aidd-mvp/data/raw 与 processed；应在 make_figures_mvp.py 之后运行以保留最终口径。
# 完整对应关系与运行顺序：README.md
"""用 aidd-mvp 工程（run_all.sh 一键复现）的数据，重画 m2 / m4 / m5 / m6。

旧版四图基于早期试跑，数字口径与最终工程不一致；本文全部数字统一到工程口径。
"""
import os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Hiragino Sans GB"
plt.rcParams["axes.unicode_minus"] = False
GREEN, BLUE, RED, GREY, AMBER = "#0FA87E", "#2B5F8F", "#C0392B", "#9A9A9A", "#F5A623"
INK = "#1A1A1A"

ROOT = os.path.expanduser("~/WorkBuddy/Work/aidd-mvp")
PROC = os.path.join(ROOT, "data", "processed")
RAW = os.path.join(ROOT, "data", "raw")
OUT = os.path.expanduser("~/WorkBuddy/Work/AI制药实操篇_配图")
res7 = json.load(open(f"{PROC}/07_baseline.json"))

# ============================================================
# m2: random vs scaffold split（工程口径：0.5095 / 0.6015）
# ============================================================
r1, r2 = res7["random"], res7["scaffold"]
labels = ["ECFP4+RF\nrandom split", "ECFP4+RF\nscaffold split"]
vals_mae = [r1["MAE"], r2["MAE"]]
vals_r2 = [r1["R2"], r2["R2"]]
vals_ov = [r1["scaffold_overlap"], r2["scaffold_overlap"]]
inc = (vals_mae[1] - vals_mae[0]) / vals_mae[0]

fig, axes = plt.subplots(3, 1, figsize=(8.6, 8.0))
for ax, vv, tt, ym, fmt in [
        (axes[0], vals_mae, "MAE（越低越好）", 0.80, "{:.4f}"),
        (axes[1], vals_r2, "R2（越高越好）", 0.85, "{:.4f}")]:
    ax.bar(labels, vv, color=[GREEN, RED], width=0.55)
    for i, v in enumerate(vv):
        ax.text(i, v + ym * 0.015, fmt.format(v), ha="center", fontsize=10.5, fontweight="bold", color=INK)
    ax.set_title(tt, fontsize=11, color=INK)
    ax.set_ylim(0, ym); ax.grid(axis="y", ls=":", alpha=0.4); ax.set_axisbelow(True)
    ax.tick_params(labelsize=8.5)
axes[0].annotate("", xy=(1, vals_mae[1]), xytext=(1, vals_mae[0]),
                 arrowprops=dict(arrowstyle="-|>", color=RED, lw=1.6))
axes[0].text(1.03, (vals_mae[0] + vals_mae[1]) / 2, f"+{inc:.0%}", fontsize=9.5, color=RED, fontweight="bold")

axes[2].bar(labels, vals_ov, color=[GREEN, RED], width=0.55)
for i, v in enumerate(vals_ov):
    axes[2].text(i, v + 15, str(v), ha="center", fontsize=10.5, fontweight="bold", color=INK)
axes[2].set_title("测试集骨架在训练集出现过的数量", fontsize=10.5, color=INK)
axes[2].set_ylim(0, 780); axes[2].grid(axis="y", ls=":", alpha=0.4); axes[2].set_axisbelow(True)
axes[2].tick_params(labelsize=8.5)

fig.suptitle("同一个模型、同一批分子，只换切分方式：指标掉了一截（aidd-mvp 工程实跑，n=10 468）",
             fontsize=12.5, fontweight="bold", color=INK, y=0.985)
plt.tight_layout(rect=[0, 0, 1, 0.965])
plt.savefig(f"{OUT}/m2_split.png", dpi=170, bbox_inches="tight", facecolor="white")
plt.close(); print("m2 done")

# ============================================================
# m4: 活性-QED 张力 + 帕累托（admet_panel.csv）
# ============================================================
full = pd.read_csv(f"{PROC}/admet_panel.csv")
fig, axes = plt.subplots(2, 1, figsize=(8.6, 7.2))
ax = axes[0]
sub = full.dropna(subset=["pIC50_median", "QED"])
sub = sub[(sub.pIC50_median >= 3) & (sub.pIC50_median <= 12)]
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
ax.annotate(f"QED 均值\n强 {strong.QED.mean():.3f} / 弱 {weak.QED.mean():.3f}", xy=(9.55, 0.08),
            fontsize=9, color=RED, ha="center",
            bbox=dict(boxstyle="round,pad=0.35", fc="#FBF8F1", ec=RED, lw=1))
ax.annotate(f"MW 均值\n强 {strong.MW.mean():.0f} / 弱 {weak.MW.mean():.0f}", xy=(4.3, 0.90),
            fontsize=9, color=BLUE, ha="center",
            bbox=dict(boxstyle="round,pad=0.35", fc="#F2F4F6", ec=BLUE, lw=1))

ax = axes[1]
s2 = full[full.clean].dropna(subset=["pIC50_median", "QED"])   # 与 10_admet.py 口径一致：clean 全体
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
ax.set_title("没有唯一的「最佳分子」，只有一条前沿", fontsize=11.5, color=INK, fontweight="bold")
ax.legend(fontsize=9.5, framealpha=0.9); ax.grid(ls=":", alpha=0.35); ax.set_axisbelow(True)
plt.tight_layout(rect=[0, 0, 1, 0.965])
plt.savefig(f"{OUT}/m4_qed_pareto.png", dpi=170, bbox_inches="tight", facecolor="white")
plt.close(); print("m4 done | 强活性组 QED", round(strong.QED.mean(),3), "弱", round(weak.QED.mean(),3),
                   "| pareto", int(keep.sum()))

# ============================================================
# m5: 主动学习三策略（iteration_history.csv）
# ============================================================
hist = pd.read_csv(f"{PROC}/iteration_history.csv")
fig, ax = plt.subplots(figsize=(8.6, 5.2))
style = {"uncertainty": (GREEN, "o", "挑模型最不确定的（uncertainty）"),
         "random": (GREY, "s", "随机挑（random）"),
         "greedy": (RED, "^", "只挑预测最强的（greedy）")}
for k, (c, m, lab) in style.items():
    d = hist[hist.strategy == k].sort_values("n_labeled")
    xs, ys = d.n_labeled.values, d.MAE.values
    ax.plot(xs, ys, marker=m, color=c, lw=2.0, ms=7, label=lab, zorder=3)
    ax.annotate(f"{ys[-1]:.3f}", (xs[-1], ys[-1]), textcoords="offset points", xytext=(8, -3),
                fontsize=9.5, color=c, fontweight="bold")
ax.set_xlabel("已标记分子数（每轮新增 300）", fontsize=10.5, color="#3F3F3F")
ax.set_ylabel("测试集 MAE（越低越好）", fontsize=10.5, color="#3F3F3F")
ax.set_title("同样做 6 轮实验，挑谁去测，结果差很多（aidd-mvp 工程实跑）",
             fontsize=12, color=INK, fontweight="bold", pad=28)
ax.legend(fontsize=9.5, framealpha=0.95, loc=(0.33, 0.55))
ax.grid(ls=":", alpha=0.35); ax.set_axisbelow(True)
ax.annotate("uncertainty 前 4 轮落后，\n第 6 轮才反超（1800 样本）",
            xy=(1800, 0.747), xytext=(1150, 0.70), fontsize=9, color=GREEN,
            arrowprops=dict(arrowstyle="->", color=GREEN, lw=1.3),
            bbox=dict(boxstyle="round,pad=0.35", fc="#F0F9F5", ec=GREEN, lw=1))
ax.annotate("greedy 全程垫底：\nMAE 几乎不动（0.90 → 0.88）",
            xy=(1500, 0.8945), xytext=(700, 0.955), fontsize=9, color=RED,
            arrowprops=dict(arrowstyle="->", color=RED, lw=1.3),
            bbox=dict(boxstyle="round,pad=0.35", fc="#FFF5F5", ec=RED, lw=1))
ax.set_ylim(0.70, 1.00); ax.set_xlim(280, 2080)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig(f"{OUT}/m5_active_learning.png", dpi=170, bbox_inches="tight", facecolor="white")
plt.close(); print("m5 done")

# ============================================================
# m6: 数据体检（raw 分布 + 测量分歧）
# ============================================================
raw = pd.read_csv(f"{RAW}/raw_activity.csv")
v = raw[(raw.standard_units == "nM") & (raw.standard_relation == "=")]["standard_value"].values
p_raw = 9 - np.log10(v[v > 0])
p_raw = p_raw[~np.isnan(p_raw)]
n_dirty = int((p_raw > 11).sum())

df = pd.read_csv(f"{PROC}/dataset.csv")
m2_ = df[df.n_measurements >= 2]
spread = (m2_.pIC50_max - m2_.pIC50_min).values
med_spread = float(np.median(spread))

fig, axes = plt.subplots(2, 1, figsize=(8.6, 6.4))
ax = axes[0]
ax.hist(p_raw, bins=70, color=GREEN, alpha=0.75, edgecolor="white", linewidth=0.4)
ax.set_xlabel("pIC50（清洗前，全部 18 381 条等值记录）", fontsize=10.5, color="#3F3F3F")
ax.set_ylabel("记录数", fontsize=10.5, color="#3F3F3F")
ax.set_title("pIC50 分布：右边那条尾巴有问题", fontsize=11.5, color=INK, fontweight="bold")
ax.axvspan(11, 18, color=RED, alpha=0.16)
ax.annotate(f"pIC50 > 11\n= IC50 < 0.01 nM\n物理上不可能", xy=(12.6, 300), fontsize=9, color=RED,
            ha="center", bbox=dict(boxstyle="round,pad=0.35", fc="#FFF5F5", ec=RED, lw=1))
ax.text(14.3, 90, f"{n_dirty} 条脏数据\n藏在这里（最高一条\npIC50 = 17.30）\nTop1 排序会毁在这",
        fontsize=8.5, color=RED, ha="center", linespacing=1.6)
ax.grid(axis="y", ls=":", alpha=0.35); ax.set_axisbelow(True)

ax = axes[1]
ax.hist(np.clip(spread, 0, 4), bins=50, color=BLUE, alpha=0.75, edgecolor="white", linewidth=0.4)
ax.axvline(med_spread, color=RED, lw=1.6, ls="--")
ax.text(med_spread + 0.07, ax.get_ylim()[1] * 0.75,
        f"中位数分歧 {med_spread:.2f}\n个对数单位", fontsize=9, color=RED)
ax.set_xlabel("同一分子多次测量，pIC50 的最大-最小跨度", fontsize=10, color="#3F3F3F")
ax.set_ylabel("分子数", fontsize=10.5, color="#3F3F3F")
ax.set_title(f"{len(m2_)/len(df):.1%} 的分子被测过不止一次，而且测得不一样", fontsize=11.5, color=INK, fontweight="bold")
ax.grid(axis="y", ls=":", alpha=0.35); ax.set_axisbelow(True)
plt.tight_layout(rect=[0, 0, 1, 0.965])
plt.savefig(f"{OUT}/m6_data_qc.png", dpi=170, bbox_inches="tight", facecolor="white")
plt.close(); print("m6 done | dirty", n_dirty, "| med spread", round(med_spread, 2))
