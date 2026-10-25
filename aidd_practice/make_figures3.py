# -*- coding: utf-8 -*-
# 文章关联：AI制药实践篇_公众号文章.html
# 脚本职责：配图：v3：修订 fig1 / fig4 / fig5
# 使用注解：历史增量版本，不是五套完整图；同名图会被后运行的版本覆盖。
# 完整对应关系与运行顺序：README.md
"""实践篇 v3：修复 fig1 / fig4 / fig5 的错乱问题"""
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

# ================= fig1：七站闭环（底部注释下移） =================
fig, ax = plt.subplots(figsize=(7.2, 7.8), dpi=150)
fig.patch.set_facecolor(BG)
ax.set_xlim(0, 10); ax.set_ylim(0, 10.9); ax.axis("off")
ax.text(0.2, 10.55, "七站不是瀑布流，是螺旋上升的闭环", fontsize=16.5, fontweight="bold", color=INK)
ax.text(0.2, 10.15, "第七站的实验结果，可以推翻并回退到 ② ~ ⑥ 的任意一站", fontsize=10.5, color=GREY)

steps = [("① 疾病致病机制", "病是哪一步坏了？"),
         ("② 靶点发现与验证", "动哪个蛋白能好？"),
         ("③ 蛋白结构与分析", "长什么样？口袋在哪？"),
         ("④ 苗头化合物发现", "谁能结合上去？"),
         ("⑤ 命中化合物优化", "值得做吗？"),
         ("⑥ ADMET 与成药性", "活得下来吗？"),
         ("⑦ 干湿闭环验证", "模型说的对吗？")]
y = 9.5
for i, (t, q) in enumerate(steps):
    accent = GREEN if i < 2 else (BLUE if i < 4 else (AMBER if i < 6 else RED))
    ax.add_patch(FancyBboxPatch((1.0, y - 0.85), 7.0, 0.82, boxstyle="round,pad=0.01,rounding_size=0.06",
                                linewidth=1.2, edgecolor="#DDE3E6", facecolor=CARD, zorder=1))
    ax.add_patch(Rectangle((1.0, y - 0.85), 0.09, 0.82, facecolor=accent, zorder=2))
    ax.text(1.35, y - 0.36, t, fontsize=12.5, fontweight="bold", color=INK, va="center", zorder=3)
    ax.text(7.6, y - 0.36, q, fontsize=10.5, color="#5B6B70", va="center", ha="right", zorder=3)
    if i < len(steps) - 1:
        ax.annotate("", xy=(4.5, y - 0.93), xytext=(4.5, y - 0.82),
                    arrowprops=dict(arrowstyle="-|>", color="#C7D2D6", lw=1.6))
    y -= 1.12

ax.plot([0.72, 0.72], [2.4, 9.35], color=RED, lw=1.8, ls="--", zorder=4)
ax.annotate("", xy=(0.72, 9.45), xytext=(0.72, 2.4),
            arrowprops=dict(arrowstyle="-|>", color=RED, lw=1.8, ls="--"))
ax.text(0.30, 5.8, "真实数据\n回流更新", fontsize=10, color=RED, rotation=90, va="center", ha="center")
for yy in [9.25, 8.15, 7.05]:
    ax.plot([0.72, 0.95], [yy, yy], color=RED, lw=1.0, ls=":")
ax.text(0.2, 1.25, "螺旋上升：每一轮回环，模型与假说都应比上一轮更准", fontsize=10.5, color=RED)
ax.text(0.2, 0.45, "作者自制 · 七站顺序依据本系列前序课程内容整理", fontsize=9, color="#B8B8B8")
fig.savefig(OUT + "fig1_seven_stations.png", facecolor=BG, bbox_inches="tight", pad_inches=0.12)
plt.close(fig)

# ================= fig4：接缝（补站名） =================
seams = [("①→②", "列了 20 个候选靶点，却排不出优先级"),
         ("②→③", "靶点定了才发现没有可用结构"),
         ("③→④", "对接打分与真实亲和力相关性差"),
         ("④→⑤", "Hit 骨架太集中，优化很快撞墙"),
         ("⑤→⑥", "只优化活性，Lead 在 ADMET 阶段集体阵亡"),
         ("⑥→⑦", "低置信预测没安排实验，直接推进"),
         ("⑦→①", "结果与预期不符，却只是「换个分子再试」")]
fig, ax = plt.subplots(figsize=(7.2, 5.4), dpi=150)
fig.patch.set_facecolor(BG)
ax.set_xlim(0, 10); ax.set_ylim(0, 7.6); ax.axis("off")
ax.text(0.2, 7.25, "七站之间的六道接缝：项目最容易漏的地方", fontsize=15.5, fontweight="bold", color=INK)
ax.text(0.2, 6.82, "每一步的输出都是下一步的输入，接缝处的信息损耗是失败头号原因", fontsize=10.5, color=GREY)

xs = np.linspace(1.1, 8.9, 7)
station_labels = ["①", "②", "③", "④", "⑤", "⑥", "⑦"]
for i, x in enumerate(xs):
    ax.add_patch(FancyBboxPatch((x - 0.5, 5.45), 1.0, 0.62, boxstyle="round,pad=0.01,rounding_size=0.08",
                                linewidth=1.2, edgecolor="#DDE3E6", facecolor=CARD, zorder=1))
    ax.text(x, 5.76, station_labels[i], fontsize=13, fontweight="bold", color=INK,
            ha="center", va="center", zorder=3)
for i, x in enumerate(xs[:-1]):
    ax.plot([x + 0.52, x + 0.98], [5.76, 5.76], color=GREY, lw=1.2, zorder=0)
    ax.plot([x + 0.75, x + 0.75], [5.45, 6.07], color=RED, lw=1.6, ls="--", zorder=2)
    ax.text(x + 0.75, 6.18, seams[i][0], fontsize=8.5, color=RED, ha="center", va="bottom")

for i, (lab, text) in enumerate(seams):
    yy = 5.00 - i * 0.66
    ax.text(0.35, yy, lab, fontsize=11, fontweight="bold", color=RED, va="center")
    ax.text(1.45, yy, text, fontsize=11.5, color="#3F3F3F", va="center")
ax.text(0.35, 0.30, "记住：每一步都要问「我交出去的东西，下一站真的能直接用吗？」",
        fontsize=10.5, color=GREEN, fontweight="bold")
fig.savefig(OUT + "fig4_seams.png", facecolor=BG, bbox_inches="tight", pad_inches=0.14)
plt.close(fig)

# ================= fig5：两条更难的路（彻底重排） =================
fig = plt.figure(figsize=(7.2, 9.2), dpi=150)
fig.patch.set_facecolor(BG)

# --- 上半：两个案例卡片（gs 两列） ---
gs_top = fig.add_gridspec(1, 2, left=0.04, right=0.96, top=0.885, bottom=0.545, wspace=0.14)
ax.text = None
axt = fig.add_axes([0, 0, 1, 1]); axt.axis("off")
axt.set_xlim(0, 1); axt.set_ylim(0, 1)
axt.text(0.03, 0.955, "两条更难的路：同一条流水线，两种打法", fontsize=15.5, fontweight="bold",
         color=INK, transform=axt.transAxes)

def case_card(ax, tag, title, lines, accent):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.add_patch(FancyBboxPatch((0.02, 0.04), 0.96, 0.92, boxstyle="round,pad=0.01,rounding_size=0.03",
                                linewidth=1.4, edgecolor=accent, facecolor="#FBFCFC"))
    ax.text(0.08, 0.885, tag, fontsize=9.5, color=accent, fontweight="bold")
    ax.text(0.08, 0.785, title, fontsize=13, fontweight="bold", color=INK)
    n = len(lines)
    for i, ln in enumerate(lines):
        ax.text(0.08, 0.635 - i * (0.55 / max(n - 1, 1)), ln, fontsize=9.2, color="#4A4A4A", va="top")

axA = fig.add_subplot(gs_top[0, 0])
case_card(axA, "案例 A · 把「不可成药」变成可成药", "KRAS G12C", [
    "难点：表面平滑无深口袋，",
    "　　　与 GTP 亲和力皮摩尔级",
    "破局1：G12C 带来非天然 Cys12",
    "破局2：动态隐性口袋 switch-II",
    "打法：共价锁定失活态(OFF)",
    "结果：Sotorasib 2021-05 获批",
    "　　：ORR 37.1% / mPFS 6.8 月",
], GREEN)

axB = fig.add_subplot(gs_top[0, 1])
case_card(axB, "案例 B · 当「抓手」自己消失了", "EGFR C797S", [
    "背景：三代药共价结合 Cys797",
    "耐药：C797S 后共价弹头失效",
    "打法：强效非共价 / 变构抑制",
    "硬约束：对野生型必须高选择性",
    "现状：第四代 TKI 多在 I/II 期",
    "　　：截至 2026 无获批上市",
    "　　：C797S 约占耐药 10%~25%",
], AMBER)

# --- 下半：ΔΔG 选择性窗口 ---
ax = fig.add_axes([0.07, 0.06, 0.86, 0.40])
ax.set_facecolor(CARD)
for spine in ax.spines.values():
    spine.set_color("#DDE3E6")
ax.text(0.0, 1.10, "选择性窗口怎么算：同一个分子，在野生型和突变体上的结合能之差",
        transform=ax.transAxes, fontsize=11.5, fontweight="bold", color=INK)
ax.text(0.0, 1.005, "ΔΔG(选择性) = ΔG(野生型) - ΔG(突变体)　——　这个值越大，治疗窗口越宽",
        transform=ax.transAxes, fontsize=10, color="#B7791F")

names = ["野生型\n(要弱)", "突变体\n(要强)"]
vals = [-8.2, -12.5]
bars = ax.bar([0, 1], vals, width=0.5, color=["#C7D2D6", GREEN], zorder=3)
for b, v in zip(bars, vals):
    ax.text(b.get_x() + b.get_width() / 2, v + 0.45, "%.1f" % v, ha="center",
            fontsize=10.5, fontweight="bold", color=INK)
ax.annotate("", xy=(1.30, -12.5), xytext=(1.30, -8.2),
            arrowprops=dict(arrowstyle="<->", color=RED, lw=1.8))
ax.text(1.42, -10.3, "ΔΔG = -8.2 - (-12.5)\n= +4.3 kcal/mol\n(越大越好)",
        fontsize=9.5, color=RED, va="center")
ax.axhline(0, color="#C7D2D6", lw=1.0)
ax.set_xticks([0, 1]); ax.set_xticklabels(names, fontsize=10)
ax.set_ylabel("结合自由能 ΔG (kcal/mol)", fontsize=9.5, color="#3F3F3F")
ax.set_xlim(-0.55, 2.35)
ax.set_ylim(-14.5, 0.8)
ax.tick_params(labelsize=9, colors=GREY)
ax.text(0.0, -0.30, "FEP 最适合算「小基团局部微扰」的 ΔΔG；骨架跃迁请交给对接 + MM-GBSA，别迷信 FEP。",
        transform=ax.transAxes, fontsize=9.5, color="#3F3F3F")
ax.text(0.0, -0.44, "案例数据来源：CodeBreaK 100/200、KRYSTAL-1 等公开披露；ΔG 为示意值",
        transform=ax.transAxes, fontsize=8.5, color="#B8B8B8")

fig.savefig(OUT + "fig5_two_wars.png", facecolor=BG, bbox_inches="tight", pad_inches=0.14)
plt.close(fig)
print("done v3")
