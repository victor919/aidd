# -*- coding: utf-8 -*-
"""⑬ 八张配图：全部基于 e23-e26 固定 json 结果，无预填趋势。"""
import os, json, random
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
OUT = "/Users/dujiepeng/WorkBuddy/Work/AI制药第13篇_配图"
os.makedirs(OUT, exist_ok=True)
DARK = "#1A1A1A"; GREEN = "#0FA87E"; BLUE = "#2B5F8F"; RED = "#C0392B"; GREY = "#9A9A9A"
plt.rcParams["font.sans-serif"] = ["Hiragino Sans GB"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["axes.edgecolor"] = "#D5DBDF"
plt.rcParams["figure.dpi"] = 165

def save(fig, name):
    p = os.path.join(OUT, name)
    fig.savefig(p, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", name, flush=True)

e24 = json.load(open(os.path.join(RES, "e24_offline.json")))
e25 = json.load(open(os.path.join(RES, "e25_closedloop.json")))
e26 = json.load(open(os.path.join(RES, "e26_disagreement.json")))
pool = json.load(open(os.path.join(RES, "e23_pool.json")))

# ---------- 图1 开环 vs 闭环 ----------
fig, ax = plt.subplots(figsize=(9.0, 3.4)); ax.axis("off")
ax.set_xlim(0, 10); ax.set_ylim(0, 4)
ax.text(0.2, 3.25, "开环评分（第 11、12 篇）", fontsize=12, color=BLUE, weight="bold")
ax.text(0.2, 1.15, "闭环搜索（第 13 篇）", fontsize=12, color=RED, weight="bold")
def box(x, y, w, h, txt, fc):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.06", fc=fc, ec="none"))
    ax.text(x+w/2, y+h/2, txt, ha="center", va="center", fontsize=10, color="white")
def arrow(x1, y1, x2, y2, color):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=13,
                                 color=color, lw=1.6))
box(.4, 2.5, 2.0, .55, "已知分子库", BLUE); box(3.3, 2.5, 2.0, .55, "评分器", BLUE)
box(6.2, 2.5, 2.0, .55, "排序结果", BLUE)
arrow(2.45, 2.78, 3.25, 2.78, BLUE); arrow(5.35, 2.78, 6.15, 2.78, BLUE)
box(.4, .4, 2.0, .55, "生成器", RED); box(3.3, .4, 2.0, .55, "候选", RED)
box(6.2, .4, 2.0, .55, "评分器", RED); box(8.65, .4, 1.1, .55, "保留", RED)
arrow(2.45, .68, 3.25, .68, RED); arrow(5.35, .68, 6.15, .68, RED); arrow(8.25, .68, 8.6, .68, RED)
ax.add_patch(FancyArrowPatch((9.2, .95), (1.5, .95), connectionstyle="arc3,rad=-0.18",
                             arrowstyle="-|>", mutation_scale=13, color=RED, lw=1.8))
ax.text(5.0, 1.55, "反馈：评分器在闭环内，且是自家的", ha="center", fontsize=10.5, color=RED)
save(fig, "图01_开环vs闭环.png")

# ---------- 图2 三部件 ----------
fig, ax = plt.subplots(figsize=(9.0, 2.6)); ax.axis("off")
ax.set_xlim(0, 10.1); ax.set_ylim(0, 2.3)
cols = [("#E8F5F0", "表示层", "SMILES / SELFIES / 片段 / 图\n本篇：BRICS 片段"),
        ("#EAF0F8", "提议动作", "随机突变 / 重组 / 爬山 / RL\n本篇：局部片段替换"),
        ("#FBEFEC", "保留规则", "随机 / 按预测分数 / 带约束\n本篇：A、B、C 三臂")]
for i, (fc, t, s) in enumerate(cols):
    x = .3 + i*3.25
    ax.add_patch(FancyBboxPatch((x, .3), 3.0, 1.75, boxstyle="round,pad=0,rounding_size=0.08", fc=fc, ec="none"))
    ax.text(x+.15, 1.72, t, fontsize=11.5, color=DARK, weight="bold")
    ax.text(x+.15, 1.05, s, fontsize=9.5, color="#444")
    if i < 2:
        arrow(x+3.02, 1.15, x+3.23, 1.15, GREY)
ax.text(4.5, .03, "三臂共用表示层与提议动作，只有保留规则不同", ha="center", fontsize=10, color=RED)
save(fig, "图02_三部件.png")

# ---------- 图3/图4 离线对照：上下对齐，共享横轴（禁双纵轴） ----------
runs = e24["runs"]
fig, axes = plt.subplots(2, 1, figsize=(9.2, 6.2), sharex=True)
for ax, key, ylab, title in zip(
        axes, ["pred_med", "nn_med"],
        ["RF 预测 pIC50 中位", "最近邻 Tanimoto 中位"],
        ["A｜按分数选择的候选，预测分数更高", "B｜但最近邻相似度：有升有降，不是稳定剪刀差"]):
    for arm, col, lab in [("A", BLUE, "A 随机保留"), ("B", RED, "B 按分数保留")]:
        vals = [r[arm][key] for r in runs]
        xs = np.arange(len(vals))
        ax.plot(xs, vals, "o-", color=col, lw=1.7, ms=5, label=lab)
        ax.plot(xs, vals, "o", color=col, ms=5)
    ax.set_ylabel(ylab, fontsize=10, color=DARK)
    ax.set_title(title, fontsize=11, color=DARK, loc="left")
    ax.grid(alpha=.25)
axes[0].legend(fontsize=9, loc="upper left", ncol=2)
axes[1].set_xticks(np.arange(len(runs)))
axes[1].set_xticklabels([f"{r['seed']}/{r['start']}" for r in runs], fontsize=8.5)
axes[1].set_xlabel("随机种子 / 起点（9 组）", fontsize=10)
fig.suptitle("离线对照：同池随机选择 vs 按分数排序（保留 100/池）", fontsize=12.5, color=DARK, y=.995)
fig.tight_layout()
save(fig, "图03_离线分数.png")  # 两面板合图：正文图3/图4合并为一张上下对齐图
# 同时单独存两张，正文按需取用
fig, ax = plt.subplots(figsize=(9.2, 3.1))
for arm, col, lab in [("A", BLUE, "A 随机保留"), ("B", RED, "B 按分数保留")]:
    ax.plot(np.arange(len(runs)), [r[arm]["pred_med"] for r in runs], "o-",
            color=col, lw=1.7, ms=5, label=lab)
ax.set_ylabel("RF 预测 pIC50 中位", fontsize=10); ax.grid(alpha=.25)
ax.set_xticks(np.arange(len(runs))); ax.set_xticklabels([f"{r['seed']}/{r['start']}" for r in runs], fontsize=8.5)
ax.set_xlabel("随机种子 / 起点（9 组）", fontsize=10); ax.legend(fontsize=9)
ax.set_title("分数在涨：9 组里 B 臂全部高于 A 臂", fontsize=11, color=DARK, loc="left")
save(fig, "图03b_离线分数.png")
fig, ax = plt.subplots(figsize=(9.2, 3.1))
for arm, col, lab in [("A", BLUE, "A 随机保留"), ("B", RED, "B 按分数保留")]:
    ax.plot(np.arange(len(runs)), [r[arm]["nn_med"] for r in runs], "o-",
            color=col, lw=1.7, ms=5, label=lab)
ax.axhline(.5, color=RED, ls="--", lw=1, alpha=.7)
ax.text(.1, .505, "最近邻 0.5（第11篇 域内/域外的诊断分箱）", fontsize=8.5, color=RED)
ax.set_ylabel("最近邻 Tanimoto 中位", fontsize=10); ax.grid(alpha=.25)
ax.set_xticks(np.arange(len(runs))); ax.set_xticklabels([f"{r['seed']}/{r['start']}" for r in runs], fontsize=8.5)
ax.set_xlabel("随机种子 / 起点（9 组）", fontsize=10); ax.legend(fontsize=9)
ax.set_title("相似度却有升有降：本次未观察到稳定剪刀差", fontsize=11, color=DARK, loc="left")
save(fig, "图04_离线相似度.png")

# ---------- 图5 闭环：按轮次对齐（每轮提议数相同），跨 9 次运行的中位数 + 四分位带 ----------
runs25 = e25["runs"]
fig, axes = plt.subplots(2, 1, figsize=(9.2, 6.4), sharex=True)
rounds = list(range(1, 9))
for ax, key, ylab in zip(axes, ["pred_med", "nn_med"],
                         ["保留候选预测 pIC50 中位", "保留候选最近邻 Tanimoto 中位"]):
    for arm, col, lab in [("A", BLUE, "A 随机保留"), ("B", RED, "B 按分数保留"), ("C", GREEN, "C 带 AD 约束")]:
        sub = [r for r in runs25 if r["arm"] == arm and len(r["rounds"]) == 8]
        meds = [np.median([r["rounds"][i][key] for r in sub]) for i in [x-1 for x in rounds]]
        q25 = [np.quantile([r["rounds"][i][key] for r in sub], .25) for i in [x-1 for x in rounds]]
        q75 = [np.quantile([r["rounds"][i][key] for r in sub], .75) for i in [x-1 for x in rounds]]
        xs = np.array(rounds) + {"A": -.12, "B": 0, "C": .12}[arm]
        ax.plot(xs, meds, "o-", color=col, lw=1.8, ms=4.5, label=lab)
        ax.fill_between(xs, q25, q75, color=col, alpha=.10, lw=0)
    ax.set_ylabel(ylab, fontsize=10); ax.grid(alpha=.25)
    ax.set_xticks(rounds)
axes[0].legend(fontsize=9, ncol=3, loc="upper left")
axes[0].set_title("闭环三臂：同一提议函数下，只改保留规则（9 起点-种子 × 3 臂）",
                  fontsize=11.5, color=DARK, loc="left")
axes[1].set_xlabel("搜索轮次（每轮同预算：5 束 × 10 候选；阴影为四分位区间）", fontsize=10)
fig.tight_layout()
save(fig, "图05_闭环轨迹.png")

# ---------- 图6 分歧—误差（测试集有标签） ----------
t = e26["test_set"]
fig, ax = plt.subplots(figsize=(9.0, 3.8))
labs = [b["bin"] for b in t["bins"]]
vals = [b["max_err_med"] for b in t["bins"]]
frac = [b["frac_max_err_gt1"] for b in t["bins"]]
x = np.arange(len(labs))
b = ax.bar(x, vals, color=[GREEN, "#7FC7AE", "#E9B384", RED], width=.55)
for xi, v, f in zip(x, vals, frac):
    ax.text(xi, v+.025, f"{v:.3f}", ha="center", fontsize=9, color=DARK)
    ax.text(xi, .06, f"误差>1: {f:.0%}", ha="center", fontsize=8.5, color="white")
ax.set_xticks(x); ax.set_xticklabels(labs, fontsize=10)
ax.set_ylabel("两模型中较大误差的中位数", fontsize=10)
ax.set_title(f"测试集（n={t['n']}，有真值）：分歧分箱 vs 误差  ρ={t['spearman_disagree_vs_max_err']}",
             fontsize=11, color=DARK, loc="left")
ax.grid(axis="y", alpha=.25)
save(fig, "图06_分歧与误差.png")

# ---------- 图7 四道闸门 ----------
fig, ax = plt.subplots(figsize=(9.0, 3.2)); ax.axis("off")
ax.set_xlim(0, 10.1); ax.set_ylim(0, 2.5)
gates = [("1 化学有效性", "RDKit 解析、去重\n拦：乱码 / 非法结构", "#E8F5F0"),
         ("2 分布与属性", "最近邻、SA、QED\n拦：极端外推 / 明显类药性差", "#EAF0F8"),
         ("3 正交证据", "对接姿势审阅（本篇不跑）\n拦：结构冲突；不给活性加分", "#FBF3E8"),
         ("4 实验", "合成 + 活性测定\n拦：前面所有自欺", "#FBEFEC")]
for i, (t, s, fc) in enumerate(gates):
    x = .2 + i*2.45
    ax.add_patch(FancyBboxPatch((x, .5), 2.25, 1.8, boxstyle="round,pad=0,rounding_size=0.08", fc=fc, ec="none"))
    ax.text(x+.12, 1.98, t, fontsize=10.5, color=DARK, weight="bold")
    ax.text(x+.12, 1.28, s, fontsize=8.8, color="#444", va="top")
    if i < 3: arrow(x+2.26, 1.4, x+2.43, 1.4, GREY)
ax.text(4.5, .12, "每道闸门都要写明：拦得住什么，拦不住什么", ha="center", fontsize=10.5, color=RED)
save(fig, "图07_四道闸门.png")

# ---------- 图8 系列地图 ----------
fig, ax = plt.subplots(figsize=(9.0, 2.0)); ax.axis("off")
ax.set_xlim(0, 10.1); ax.set_ylim(0, 1.6)
steps = [("9\n对象", "#E8F5F0"), ("10\n摆放与评分", "#E8F5F0"), ("11\n学习与可信", "#EAF0F8"),
         ("12\n表示与代价", "#EAF0F8"), ("13\n提出与验收", "#FBEFEC")]
for i, (t, fc) in enumerate(steps):
    x = .2 + i*1.95
    ax.add_patch(FancyBboxPatch((x, .35), 1.7, 1.1, boxstyle="round,pad=0,rounding_size=0.08", fc=fc, ec="none"))
    ax.text(x+.85, .9, t, ha="center", va="center", fontsize=10, color=DARK)
    if i < 4: arrow(x+1.72, .9, x+1.93, .9, GREY)
ax.text(4.5, .05, "共同纪律：每个输出先问它回答什么问题、在哪失效", ha="center", fontsize=10, color=RED)
save(fig, "图08_系列地图.png")
print("ALL DONE")
