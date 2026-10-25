# -*- coding: utf-8 -*-
# 文章关联：AI制药实践篇_公众号文章.html
# 脚本职责：配图：v4：修订 fig1，新增 fig7 三层支撑
# 使用注解：历史增量版本，不是五套完整图；同名图会被后运行的版本覆盖。
# 完整对应关系与运行顺序：README.md
"""实践篇 v4 新图：七站口诀(fig1) + 一条主线三层支撑(fig7)"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
import numpy as np

OUT = "/Users/dujiepeng/WorkBuddy/Work/AI制药实践篇_配图/"
plt.rcParams["font.family"] = "Hiragino Sans GB"
plt.rcParams["axes.unicode_minus"] = False
GREEN, BLUE, INK, GREY = "#0FA87E", "#2B5F8F", "#1A1A1A", "#9A9A9A"
BG, CARD, RED, AMBER = "#FFFFFF", "#F7F9FA", "#E53E3E", "#F5A623"

# ============ fig1：七站 + 口诀 ============
STATIONS = [
    ("① 疾病致病机制", "病是哪一步坏了？", "找因果", GREEN),
    ("② 靶点发现与验证", "动哪个蛋白能好？", "找可干预节点", GREEN),
    ("③ 蛋白结构与分析", "长什么样？口袋在哪？", "找口袋", BLUE),
    ("④ 苗头化合物发现", "谁能结合上去？", "找起点", BLUE),
    ("⑤ 命中化合物优化", "值得做吗？", "找平衡", AMBER),
    ("⑥ ADMET 与成药性", "活得下来吗？", "找风险", RED),
    ("⑦ 干湿闭环验证", "模型说的对吗？", "找证据", RED),
]

fig, ax = plt.subplots(figsize=(7.2, 7.6), dpi=150)
fig.patch.set_facecolor(BG)
ax.set_xlim(0, 10); ax.set_ylim(0, 10.6); ax.axis("off")
ax.text(0.2, 10.25, "七站流水线：一站只做好一件事", fontsize=16.5, fontweight="bold", color=INK)
ax.text(0.2, 9.88, "每行右端的绿字是这一站的核心动作——七个动词串起来就是整条链", fontsize=10.5, color=GREY)

y = 9.32
for i, (t, q, act, accent) in enumerate(STATIONS):
    ax.add_patch(FancyBboxPatch((1.0, y - 0.85), 7.0, 0.82, boxstyle="round,pad=0.01,rounding_size=0.06",
                                linewidth=1.2, edgecolor="#DDE3E6", facecolor=CARD, zorder=1))
    ax.add_patch(Rectangle((1.0, y - 0.85), 0.09, 0.82, facecolor=accent, zorder=2))
    ax.text(1.35, y - 0.30, t, fontsize=12.5, fontweight="bold", color=INK, va="center", zorder=3)
    ax.text(1.35, y - 0.60, q, fontsize=10.5, color="#5B6B70", va="center", zorder=3)
    # 核心动作胶囊
    ax.add_patch(FancyBboxPatch((6.30, y - 0.68), 1.60, 0.28, boxstyle="round,pad=0.008,rounding_size=0.12",
                                linewidth=0, facecolor="#EAF7F2", zorder=2))
    ax.text(7.10, y - 0.54, act, fontsize=9.6, color=GREEN, va="center", ha="center",
            fontweight="bold", zorder=3)
    if i < len(STATIONS) - 1:
        ax.annotate("", xy=(2.0, y - 0.92), xytext=(2.0, y - 0.82),
                    arrowprops=dict(arrowstyle="-|>", color="#C7D2D6", lw=1.5))
    y -= 1.05

ax.plot([0.72, 0.72], [2.6, 9.15], color=RED, lw=1.8, ls="--", zorder=4)
ax.annotate("", xy=(0.72, 9.25), xytext=(0.72, 2.6),
            arrowprops=dict(arrowstyle="-|>", color=RED, lw=1.8, ls="--"))
ax.text(0.30, 5.9, "真实数据回流\n更新模型与假说", fontsize=9.5, color=RED, rotation=90,
        va="center", ha="center")
ax.text(0.2, 1.05, "每一站都有各自的 Go / No-Go：证据不够就停下，别把误差带进下一站",
        fontsize=10.5, color=GREEN, fontweight="bold")
ax.text(0.2, 0.40, "作者自制 · 七站顺序依据本系列前序课程内容整理", fontsize=9, color="#B8B8B8")
fig.savefig(OUT + "fig1_seven_stations.png", facecolor=BG, bbox_inches="tight", pad_inches=0.12)
plt.close(fig)

# ============ fig7：一条主线 + 三层支撑 ============
fig = plt.figure(figsize=(7.2, 5.6), dpi=150)
fig.patch.set_facecolor(BG)
axt = fig.add_axes([0, 0, 1, 1]); axt.axis("off")
axt.set_xlim(0, 10); axt.set_ylim(0, 10)
axt.text(0.25, 9.55, "一条主线，三层支撑", fontsize=16, fontweight="bold", color=INK)
axt.text(0.25, 9.10, "七站是决策链；下面三层能力托住它——任何一站塌了，先查是哪层没搭好",
         fontsize=10.5, color=GREY)

# 主线：七个小方块
axt.text(0.25, 8.45, "主线 · 七站决策链", fontsize=11.5, fontweight="bold", color=INK)
labels = ["① 机制", "② 靶点", "③ 结构", "④ 苗头", "⑤ 优化", "⑥ ADMET", "⑦ 闭环"]
short = ["找因果", "找节点", "找口袋", "找起点", "找平衡", "找风险", "找证据"]
accents = [GREEN, GREEN, BLUE, BLUE, AMBER, RED, RED]
x = 0.25
for i, (lb, act, ac) in enumerate(zip(labels, short, accents)):
    w = 1.30
    axt.add_patch(FancyBboxPatch((x, 7.55), w, 0.70, boxstyle="round,pad=0.006,rounding_size=0.08",
                                 linewidth=1.2, edgecolor="#DDE3E6", facecolor=CARD, zorder=1))
    axt.add_patch(Rectangle((x, 7.55), w, 0.06, facecolor=ac, zorder=2))
    axt.text(x + w / 2, 7.94, lb, fontsize=9.8, color=INK, ha="center", va="center",
             fontweight="bold", zorder=3)
    axt.text(x + w / 2, 7.66, act, fontsize=8.4, color=ac, ha="center", va="center", zorder=3)
    if i < 6:
        axt.annotate("", xy=(x + w + 0.075, 7.90), xytext=(x + w + 0.005, 7.90),
                     arrowprops=dict(arrowstyle="-|>", color="#C7D2D6", lw=1.4))
    x += w + 0.09

# 三层
layers = [
    ("第一层 · 生物语义层", "让疾病、基因、蛋白、通路、化合物能被机器连接",
     "MONDO / HPO / GO 本体 · UniProt / Ensembl / ChEMBL 标识 · 知识图谱 · 文献挖掘",
     GREEN),
    ("第二层 · 工具算法层", "让每一站的计算任务真的能执行",
     "序列/结构/化学工具 · AI 预测模型 · 分子生成 · 对接 · ADMET 预测 · 主动学习",
     BLUE),
    ("第三层 · 实验决策层", "最终裁判，也是模型迭代的数据来源",
     "生化/细胞/动物实验 · PK/PD 与毒理 · LIMS/ELN 数据回流 · Go / No-Go 决策",
     RED),
]
y = 6.90
for i, (title, sub, tools, ac) in enumerate(layers):
    h = 1.55
    axt.add_patch(FancyBboxPatch((0.25, y - h), 9.5, h, boxstyle="round,pad=0.012,rounding_size=0.05",
                                 linewidth=1.3, edgecolor=ac, facecolor="#FBFCFC", zorder=1))
    axt.add_patch(Rectangle((0.25, y - h), 0.10, h, facecolor=ac, zorder=2))
    axt.text(0.62, y - 0.34, title, fontsize=12.2, fontweight="bold", color=ac, va="center", zorder=3)
    axt.text(0.62, y - 0.72, sub, fontsize=10.5, color="#3F3F3F", va="center", zorder=3)
    axt.text(0.62, y - 1.10, tools, fontsize=9.2, color="#6B7A80", va="center", zorder=3,
             fontfamily="Hiragino Sans GB")
    # 向上的支撑箭头
    axt.annotate("", xy=(4.9, y + 0.16), xytext=(4.9, y + 0.02),
                 arrowprops=dict(arrowstyle="->", color="#C7D2D6", lw=1.2))
    y -= h + 0.30

axt.text(0.25, 0.42, "生物语义篇解决第一层，工具篇解决第二层，实践篇把三层串成一条能跑的链",
         fontsize=10.5, color=GREEN, fontweight="bold")
axt.text(0.25, 0.10, "作者自制 · 三层划分用于体系统筹，非行业标准术语", fontsize=8.5, color="#B8B8B8")
fig.savefig(OUT + "fig7_three_layers.png", facecolor=BG, bbox_inches="tight", pad_inches=0.12)
plt.close(fig)
print("done v4")
