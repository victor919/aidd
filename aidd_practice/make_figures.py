# -*- coding: utf-8 -*-
# 文章关联：AI制药实践篇_公众号文章.html
# 脚本职责：配图：早期七站图、苗头到先导图、干湿闭环图
# 使用注解：历史增量版本，不是五套完整图；同名图会被后运行的版本覆盖。
# 完整对应关系与运行顺序：README.md
"""AI 制药实践篇 · 配图生成（matplotlib，本机实测数据）"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

OUT = "/Users/dujiepeng/WorkBuddy/Work/AI制药实践篇_配图/"
plt.rcParams["font.family"] = "Hiragino Sans GB"
plt.rcParams["axes.unicode_minus"] = False

GREEN, BLUE, GREY, INK, BG = "#0FA87E", "#2B5F8F", "#9A9A9A", "#1A1A1A", "#FFFFFF"
CARD = "#F7F9FA"

def card(ax, x, y, w, h, num, title, question, io, tools, accent=GREEN):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.05",
                                linewidth=1.2, edgecolor="#DDE3E6", facecolor=CARD, zorder=1))
    ax.add_patch(FancyBboxPatch((x, y), 0.055, h, boxstyle="round,pad=0,rounding_size=0.02",
                                linewidth=0, facecolor=accent, zorder=2))
    ax.text(x + 0.12, y + h * 0.74, num, fontsize=15, fontweight="bold", color=accent, va="center", zorder=3)
    ax.text(x + 0.30, y + h * 0.75, title, fontsize=15.5, fontweight="bold", color=INK, va="center", zorder=3)
    ax.text(x + 0.30, y + h * 0.44, question, fontsize=11.5, color="#4A4A4A", va="center", zorder=3)
    ax.text(x + 0.30, y + h * 0.20, io, fontsize=10, color="#5B6B70", va="center", zorder=3)
    ax.text(x + 0.30, y + h * 0.03, tools, fontsize=9.5, color=GREY, va="center", zorder=3)

STATIONS = [
    ("①", "疾病致病机制分析", "这个病，到底是哪一步坏了？",
     "多组学/文献/临床样本 → 机制假说 + 关键基因通路", "GEO / TCGA · STRING · KEGG/GO · PandaOmics", GREEN),
    ("②", "靶点发现与验证", "动哪个蛋白，病能好？",
     "机制假说 → 优先靶点 + 可成药性 + 生物学验证", "DepMap · Open Targets · ChEMBL · CRISPR/类器官", GREEN),
    ("③", "蛋白质结构预测与分析", "靶点长什么样？口袋在哪？",
     "氨基酸序列 → 3D 结构 + 口袋 + 置信度(pLDDT)", "PDB · AlphaFold DB · ColabFold · fpocket", BLUE),
    ("④", "苗头化合物发现", "海量分子里，谁可能结合上去？",
     "口袋 + 化合物库 → Hit 名单（数百个）", "Vina/DiffDock · RDKit · ZINC / Enamine REAL", BLUE),
    ("⑤", "命中化合物优化", "怎么把『能结合』变成『值得做』？",
     "Hit + SAR 数据 → 先导化合物（Lead）", "FEP · Chemprop · REINVENT · RDKit SAR", "#F5A623"),
    ("⑥", "ADMET 与成药性评估", "它在人体里活得下来吗？",
     "Lead 结构 → 吸收/代谢/毒性画像 + 风险清单", "RDKit 描述符 · TDC ADMET · Caco-2/hERG/CYP", "#E5892B"),
    ("⑦", "干湿闭环验证", "模型说的对吗？下一轮该做什么？",
     "湿实验真实数据 → 更新后的模型 + 下一批待做分子", "主动学习 · 贝叶斯优化 · LIMS/ELN", "#E53E3E"),
]

# ---------- 图 1：七站流水线全景 ----------
fig, ax = plt.subplots(figsize=(7.2, 11.6), dpi=150)
fig.patch.set_facecolor(BG)
ax.set_xlim(0, 7.2); ax.set_ylim(0, 11.6); ax.axis("off")
ax.text(0.05, 11.30, "AI 制药七站流水线：每一站在回答什么问题", fontsize=17, fontweight="bold", color=INK)
ax.text(0.05, 11.05, "从疾病机制到干湿闭环 —— 每一站的「输入 → 输出」与代表工具", fontsize=11, color=GREY)

y = 10.55
h = 1.16
for i, (num, title, q, io, tools, accent) in enumerate(STATIONS):
    card(ax, 0.28, y - h, 6.64, h, num, title, q, io, tools, accent)
    if i < len(STATIONS) - 1:
        ax.annotate("", xy=(3.6, y - h - 0.10), xytext=(3.6, y - h + 0.10),
                    arrowprops=dict(arrowstyle="-|>", color="#C7D2D6", lw=1.6))
    y -= h + 0.28

ax.text(0.28, 0.32, "左色条：上游(绿·机制与靶点) → 中游(蓝·结构与筛选) → 下游(橙红·优化与验证)",
        fontsize=9.5, color=GREY)
ax.text(0.28, 0.10, "作者自制 · 数据与方法均为公开来源（RDKit / TDC / AlphaFold DB / Open Targets 等）",
        fontsize=9, color="#B8B8B8")
fig.savefig(OUT + "fig1_seven_stations.png", facecolor=BG, bbox_inches="tight", pad_inches=0.12)
plt.close(fig)

# ---------- 图 2：命中优化前后对比 ----------
metrics = [
    ("氢键供体 HBD", 5, 3, "{:.0f}"),
    ("可旋转键 RotB", 11, 6, "{:.0f}"),
    ("TPSA（极性表面积）", 149.6, 131.4, "{:.1f}"),
    ("cLogP", 0.89, 1.10, "{:.2f}"),
    ("Fsp3 (三维性)", 0.50, 0.78, "{:.2f}"),
    ("QED (成药性)", 0.328, 0.504, "{:.3f}"),
]
fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.4), dpi=150)
fig.patch.set_facecolor(BG)
fig.suptitle("命中优化到底改了什么：PF-00835231 → Nirmatrelvir", fontsize=14.5, fontweight="bold", color=INK, y=0.99)
for ax, (name, a, b, fmt) in zip(axes.ravel(), metrics):
    bars = ax.bar(["PF-00835231\n(苗头)", "Nirmatrelvir\n(上市药)"], [a, b],
                  color=["#C7D2D6", GREEN], width=0.55)
    for rect, v in zip(bars, [a, b]):
        ax.text(rect.get_x() + rect.get_width() / 2, v + (max(a, b) * 0.04),
                fmt.format(v), ha="center", fontsize=10.5, color=INK, fontweight="bold")
    ax.set_title(name, fontsize=11.5, color="#3F3F3F", pad=6)
    ax.set_ylim(0, max(a, b) * 1.32)
    ax.tick_params(labelsize=8.5, colors=GREY)
    for s in ["top", "right"]: ax.spines[s].set_visible(False)
    for s in ["left", "bottom"]: ax.spines[s].set_color("#E0E3E6")
fig.text(0.02, 0.015, "代价：Ki 0.27 nM → 3.11 nM（弱约 11 倍）　　收益：大鼠口服生物利用度 1.4% → 约 50%",
         fontsize=10.5, color="#B7791F", fontweight="bold")
fig.text(0.02, -0.03, "分子结构来自 PubChem；描述符由 RDKit 2026.03.5 本机计算；活性与 PK 数据来自 Hoffman 2020 / Owen 2021",
         fontsize=8.5, color="#B8B8B8")
fig.tight_layout(rect=[0, 0.06, 1, 0.93])
fig.savefig(OUT + "fig2_hit_to_lead.png", facecolor=BG, bbox_inches="tight", pad_inches=0.22)
plt.close(fig)

# ---------- 图 3：主动学习 vs 随机选样 ----------
n = [100, 150, 200, 250, 300, 350, 400]
rnd = [0.397, 0.437, 0.471, 0.495, 0.505, 0.545, 0.594]
act = [0.397, 0.448, 0.515, 0.545, 0.572, 0.589, 0.611]
fig, ax = plt.subplots(figsize=(7.2, 4.2), dpi=150)
fig.patch.set_facecolor(BG)
ax.plot(n, rnd, marker="o", lw=2.4, color="#C7D2D6", label="随机选样", markersize=6)
ax.plot(n, act, marker="o", lw=2.4, color=GREEN, label="主动学习（挑最不确定的）", markersize=6)
ax.axhline(0.545, ls="--", lw=1.2, color="#E5892B")
ax.annotate("同样的 R^2 ≈ 0.545：\n随机需 350 个，主动学习 250 个", xy=(250, 0.545), xytext=(150, 0.66),
            fontsize=10, color="#B7791F",
            arrowprops=dict(arrowstyle="->", color="#E5892B", lw=1.2))
ax.set_xlabel("已完成的「实验」次数（已标注分子数）", fontsize=11, color="#3F3F3F")
ax.set_ylabel("测试集 R^2", fontsize=11, color="#3F3F3F")
ax.set_title("干湿闭环：每一轮做哪些分子，比做多少个更重要", fontsize=14, fontweight="bold", color=INK, pad=12)
ax.legend(frameon=False, fontsize=10.5, loc="lower right")
ax.set_ylim(0.35, 0.70)
ax.tick_params(labelsize=10, colors=GREY)
for s in ["top", "right"]: ax.spines[s].set_visible(False)
for s in ["left", "bottom"]: ax.spines[s].set_color("#E0E3E6")
ax.text(0.0, -0.20, "数据：TDC Caco-2 渗透性数据集（910 个真实分子），ECFP4 + RandomForest，3 个随机种子平均",
        transform=ax.transAxes, fontsize=9, color="#B8B8B8")
fig.tight_layout()
fig.savefig(OUT + "fig3_dry_wet_loop.png", facecolor=BG, bbox_inches="tight", pad_inches=0.22)
plt.close(fig)

# ---------- 封面 ----------
def cover(path, w, h, title_size, sub_size):
    fig, ax = plt.subplots(figsize=(w / 150, h / 150), dpi=150)
    fig.patch.set_facecolor("#0FA87E")
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
    ax.add_patch(FancyBboxPatch((0.4, 0.4), 9.2, 9.2, boxstyle="round,pad=0,rounding_size=0.3",
                                linewidth=1.5, edgecolor="#7FD8BE", facecolor="none"))
    ax.text(5, 6.6, "AI 制药 · 实践篇", fontsize=sub_size, color="#C9F2E2", ha="center", fontweight="bold")
    ax.text(5, 5.0, title_size[0], fontsize=title_size[1], color="white", ha="center", fontweight="bold")
    ax.text(5, 3.6, "疾病机制 → 靶点 → 结构 → 苗头 → 优化 → ADMET → 干湿闭环",
            fontsize=sub_size * 0.72, color="#D9F6EC", ha="center")
    ax.text(5, 1.6, "七站流水线 · 全程实跑 · 零基础入门 AI 制药 ⑥", fontsize=sub_size * 0.66,
            color="#B7EBD8", ha="center")
    fig.savefig(path, facecolor="#0FA87E", bbox_inches="tight", pad_inches=0)
    plt.close(fig)

cover(OUT + "cover_900x383.png", 900, 383, ("从疾病到药物的七站流水线", 26), 15)
cover(OUT + "cover_500x500.png", 500, 500, ("从疾病到药物\n七站流水线", 21), 14)
print("done")
