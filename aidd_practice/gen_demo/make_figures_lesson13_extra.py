# -*- coding: utf-8 -*-
"""⑬ 补充配图：基于 e27 综合诊断结果，借鉴参考稿的审计视角，全部本机实跑。"""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results")
OUT = "/Users/dujiepeng/WorkBuddy/Work/AI制药第13篇_配图"
os.makedirs(OUT, exist_ok=True)
DARK = "#1A1A1A"; GREEN = "#0FA87E"; BLUE = "#2B5F8F"; RED = "#C0392B"; GREY = "#9A9A9A"; ORANGE = "#E67E22"
plt.rcParams["font.sans-serif"] = ["Hiragino Sans GB"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.dpi"] = 165

def save(fig, name):
    fig.savefig(os.path.join(OUT, name), bbox_inches="tight", facecolor="white")
    plt.close(fig); print("wrote", name, flush=True)

x = json.load(open(os.path.join(RES, "e27_audit_extras.json")))

# ---------- 图9 动作偏置：分子量与 SA 的逐轮变化（借鉴参考稿：动作本身也有偏置） ----------
runs = json.load(open(os.path.join(RES, "e25_closedloop.json")))["runs"]
fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.3))
for ax, key, ylab in zip(axes, ["MW_med", "SA_med"], ["保留候选分子量中位", "保留候选 SA 中位"]):
    for arm, col, lab in [("A", BLUE, "A 随机"), ("B", RED, "B 按分数")]:
        sub = [r for r in runs if r["arm"] == arm and len(r["rounds"]) == 8]
        med = [np.median([r["rounds"][i][key] for r in sub]) for i in range(8)]
        ax.plot(range(1, 9), med, "o-", color=col, lw=1.7, ms=4, label=lab)
    ax.set_xlabel("轮次", fontsize=9.5); ax.set_ylabel(ylab, fontsize=9.5, color="#3F3F3F")
    ax.grid(alpha=.25); ax.legend(fontsize=8.5)
axes[0].set_title("动作偏置：按分数臂分子量明显上升", fontsize=10.5, color=DARK, loc="left")
axes[1].set_title("SA：追高分臂同步抬高", fontsize=10.5, color=DARK, loc="left")
fig.suptitle("分数上涨伴随分子量与合成难度上升：这是结构代价，不是活性证据", fontsize=11.5, color=DARK, y=1.02)
fig.tight_layout(); save(fig, "图09_动作偏置.png")

# ---------- 图10 终点身份回查（借鉴参考稿：最高分终点可能回到已知分子） ----------
byarm = x["endpoint_identity"]["by_arm"]
fig, ax = plt.subplots(figsize=(9.0, 2.9)); ax.axis("off")
ax.set_xlim(0, 10.1); ax.set_ylim(0, 2.9)
ax.text(0.1, 2.62, "终点身份回查：规范 SMILES 精确匹配拟合池", fontsize=12, color=DARK, weight="bold")
ax.text(0.1, 2.28, "高 Tanimoto 只是线索；是否已知分子，要看结构身份。", fontsize=9.5, color=GREY)
cols = [(0.3, "A 随机保留", BLUE), (3.7, "B 按分数保留", RED), (7.1, "C 带 AD 约束", GREEN)]
for cx, name, col in cols:
    d = byarm[name[0]]
    ax.add_patch(FancyBboxPatch((cx, .45), 2.75, 1.55, boxstyle="round,pad=0.05,rounding_size=0.09",
                                fc="#F7F9FA", ec=col, lw=1.3))
    ax.text(cx + 1.375, 1.77, name, fontsize=10, color=col, ha="center", weight="bold")
    ax.text(cx + 1.375, 1.30, f"{d['n_runs']} 次运行", fontsize=9.5, color="#3F3F3F", ha="center")
    ax.text(cx + 1.375, 0.92, f"终点重现拟合池：{d['runs_with_fit_endpoint']} 次", fontsize=9.5, color=RED if d['runs_with_fit_endpoint'] else GREEN, ha="center")
    ax.text(cx + 1.375, 0.60, f"终点不同结构：{d['endpoint_unique']} 个", fontsize=9, color=GREY, ha="center")
save(fig, "图10_终点身份回查.png")

# ---------- 图11 属性筛选的代价（借鉴参考稿：能筛查不能认证） ----------
f = x["attr_filter_demo"]["by_arm"]
fig, ax = plt.subplots(figsize=(9.0, 3.2))
labels = ["终点去重", "MW 150-600", "再加 SA<=4.5", "再加 QED>=0.4"]
xx = np.arange(4); w = .25
for i, (arm, col, lab) in enumerate([("A", BLUE, "A 随机"), ("B", RED, "B 按分数"), ("C", GREEN, "C 带约束")]):
    vals = [f[arm]["endpoint_unique"], f[arm]["mw"], f[arm]["mw_sa"], f[arm]["mw_sa_qed"]]
    bars = ax.bar(xx + (i-1)*w, vals, w, color=col, alpha=.85, label=lab)
    for b, v in zip(bars, vals): ax.text(b.get_x()+b.get_width()/2, v+.25, str(v), ha="center", fontsize=8.5, color="#3F3F3F")
ax.set_xticks(xx); ax.set_xticklabels(labels, fontsize=10)
ax.set_ylabel("保留的不同结构数", fontsize=10, color="#3F3F3F"); ax.legend(fontsize=9)
ax.set_title("事后属性筛选：每加一道条件，留下多少（搜索中未使用这些条件）", fontsize=11, color=DARK, loc="left")
ax.grid(alpha=.2, axis="y"); fig.tight_layout(); save(fig, "图11_属性筛选代价.png")

print("ALL extra figures done")
