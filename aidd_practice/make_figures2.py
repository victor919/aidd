# -*- coding: utf-8 -*-
# 文章关联：AI制药实践篇_公众号文章.html
# 脚本职责：配图：v2：七站回环、接缝、双案例、学习路径
# 使用注解：历史增量版本，不是五套完整图；同名图会被后运行的版本覆盖。
# 完整对应关系与运行顺序：README.md
"""实践篇 v2 新增配图：闭环回环 / 接缝 / 双案例 / 学习路径"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
import numpy as np

OUT = "/Users/dujiepeng/WorkBuddy/Work/AI制药实践篇_配图/"
plt.rcParams["font.family"] = "Hiragino Sans GB"
plt.rcParams["axes.unicode_minus"] = False
GREEN, BLUE, INK, GREY = "#0FA87E", "#2B5F8F", "#1A1A1A", "#9A9A9A"
BG, CARD, RED, AMBER = "#FFFFFF", "#F7F9FA", "#E53E3E", "#F5A623"

# ===== 图 1（v2）：七站 + 回环 =====
fig, ax = plt.subplots(figsize=(7.2, 7.6), dpi=150)
fig.patch.set_facecolor(BG)
ax.set_xlim(0, 10); ax.set_ylim(0, 10.6); ax.axis("off")
ax.text(0.2, 10.25, "七站不是瀑布流，是螺旋上升的闭环", fontsize=16.5, fontweight="bold", color=INK)
ax.text(0.2, 9.85, "第七站的实验结果，可以推翻并回退到 ② ~ ⑥ 的任意一站", fontsize=10.5, color=GREY)

steps = [("① 疾病致病机制", "病是哪一步坏了？"),
         ("② 靶点发现与验证", "动哪个蛋白能好？"),
         ("③ 蛋白结构与分析", "长什么样？口袋在哪？"),
         ("④ 苗头化合物发现", "谁能结合上去？"),
         ("⑤ 命中化合物优化", "值得做吗？"),
         ("⑥ ADMET 与成药性", "活得下来吗？"),
         ("⑦ 干湿闭环验证", "模型说的对吗？")]
y = 9.2
for i, (t, q) in enumerate(steps):
    accent = GREEN if i < 2 else (BLUE if i < 4 else (AMBER if i < 6 else RED))
    ax.add_patch(FancyBboxPatch((1.0, y - 0.85), 7.0, 0.82, boxstyle="round,pad=0.01,rounding_size=0.06",
                                linewidth=1.2, edgecolor="#DDE3E6", facecolor=CARD, zorder=1))
    ax.add_patch(Rectangle((1.0, y - 0.85), 0.09, 0.82, facecolor=accent, zorder=2))
    ax.text(1.35, y - 0.36, t, fontsize=12.5, fontweight="bold", color=INK, va="center", zorder=3)
    ax.text(7.6, y - 0.36, q, fontsize=10.5, color="#5B6B70", va="center", ha="right", zorder=3)
    if i < len(steps) - 1:
        ax.annotate("", xy=(4.5, y - 0.90), xytext=(4.5, y - 0.80),
                    arrowprops=dict(arrowstyle="-|>", color="#C7D2D6", lw=1.6))
    y -= 1.12

# 回环：从 ⑦ 左侧绕回上方
ax.plot([0.72, 0.72], [2.15, 9.05], color=RED, lw=1.8, ls="--", zorder=4)
ax.annotate("", xy=(0.72, 9.15), xytext=(0.72, 2.15),
            arrowprops=dict(arrowstyle="-|>", color=RED, lw=1.8, ls="--"))
ax.text(0.30, 5.4, "真实数据\n回流更新", fontsize=10, color=RED, rotation=90, va="center", ha="center")
for yy, lab in [(8.95, "回到②重新排靶点"), (7.85, "回到③换构象"), (6.75, "回到④换库筛")]:
    ax.plot([0.72, 0.95], [yy, yy], color=RED, lw=1.0, ls=":")
ax.text(1.05, 1.55, "螺旋上升：每一轮回环，模型与假说都应比上一轮更准", fontsize=10.5, color=RED)
ax.text(0.2, 0.5, "作者自制 · 七站顺序依据本系列前序课程内容整理", fontsize=9, color="#B8B8B8")
fig.savefig(OUT + "fig1_seven_stations.png", facecolor=BG, bbox_inches="tight", pad_inches=0.12)
plt.close(fig)

# ===== 图 4：接缝 =====
seams = [("①→②", "列了 20 个候选靶点，却排不出优先级"),
         ("②→③", "靶点定了才发现没有可用结构"),
         ("③→④", "对接打分与真实亲和力相关性差"),
         ("④→⑤", "Hit 骨架太集中，优化很快撞墙"),
         ("⑤→⑥", "只优化活性，Lead 在 ADMET 阶段集体阵亡"),
         ("⑥→⑦", "低置信预测没安排实验，直接推进"),
         ("⑦→①", "结果与预期不符，却只是「换个分子再试」")]
fig, ax = plt.subplots(figsize=(7.2, 5.0), dpi=150)
fig.patch.set_facecolor(BG)
ax.set_xlim(0, 10); ax.set_ylim(0, 7.2); ax.axis("off")
ax.text(0.2, 6.85, "七站之间的六道接缝：项目最容易漏的地方", fontsize=15.5, fontweight="bold", color=INK)
ax.text(0.2, 6.45, "每一步的输出都是下一步的输入，接缝处的信息损耗是失败头号原因", fontsize=10.5, color=GREY)
# 链条
xs = np.linspace(1.1, 8.9, 7)
for x in xs:
    ax.add_patch(FancyBboxPatch((x - 0.5, 5.55), 1.0, 0.62, boxstyle="round,pad=0.01,rounding_size=0.08",
                                linewidth=1.2, edgecolor="#DDE3E6", facecolor=CARD, zorder=1))
for i, x in enumerate(xs[:-1]):
    ax.plot([x + 0.5, x + 0.6], [5.86, 5.86], color=GREY, lw=1.2)
for i, x in enumerate(xs[:-1]):
    ax.plot([x + 0.62, x + 0.62], [5.55, 6.17], color=RED, lw=1.6, ls="--")
    ax.text(x + 0.62, 6.30, seams[i][0], fontsize=8.5, color=RED, ha="center", va="bottom")
for i, (lab, text) in enumerate(seams[:6]):
    yy = 5.05 - i * 0.72
    ax.text(0.35, yy, lab, fontsize=11, fontweight="bold", color=RED, va="center")
    ax.text(1.45, yy, text, fontsize=11.5, color="#3F3F3F", va="center")
ax.text(0.35, 5.05 - 6 * 0.72, "⑦→①", fontsize=11, fontweight="bold", color=RED, va="center")
ax.text(1.45, 5.05 - 6 * 0.72, seams[6][1], fontsize=11.5, color="#3F3F3F", va="center")
ax.text(0.35, 0.35, "记住：每一步都要问「我交出去的东西，下一站真的能直接用吗？」", fontsize=10.5, color=GREEN, fontweight="bold")
fig.savefig(OUT + "fig4_seams.png", facecolor=BG, bbox_inches="tight", pad_inches=0.14)
plt.close(fig)

# ===== 图 5：两条更难的路 =====
fig, ax = plt.subplots(figsize=(7.2, 6.0), dpi=150)
fig.patch.set_facecolor(BG)
ax.set_xlim(0, 10); ax.set_ylim(0, 12); ax.axis("off")
ax.text(0.2, 11.6, "两条更难的路：同一条流水线，两种打法", fontsize=15.5, fontweight="bold", color=INK)

def case(x, y, tag, title, lines, accent):
    ax.add_patch(FancyBboxPatch((x, y), 4.6, 4.45, boxstyle="round,pad=0.012,rounding_size=0.08",
                                linewidth=1.4, edgecolor=accent, facecolor="#FBFCFC", zorder=1))
    ax.text(x + 0.22, y + 4.05, tag, fontsize=10.5, color=accent, fontweight="bold")
    ax.text(x + 0.22, y + 3.55, title, fontsize=12.5, fontweight="bold", color=INK)
    for i, ln in enumerate(lines):
        ax.text(x + 0.22, y + 2.95 - i * 0.62, ln, fontsize=10, color="#4A4A4A", va="top", zorder=3)

case(0.25, 6.9, "案例 A · 把「不可成药」变成可成药", "KRAS G12C", [
    "难点：表面平滑无深口袋，与 GTP 亲和力极高",
    "破局 1：G12C 突变带来非天然半胱氨酸 Cys12",
    "破局 2：发现动态的隐性口袋 switch-II pocket",
    "打法：共价锁定失活态（GDP-bound OFF）",
    "结果：Sotorasib 2021-05 获 FDA 加速批准",
    "II 期 ORR 37.1% / mPFS 6.8 月（126 例）",
], GREEN)
case(5.15, 6.9, "案例 B · 当「抓手」自己消失了", "EGFR C797S", [
    "背景：三代药奥希替尼共价结合 Cys797",
    "耐药：C797S 突变，半胱氨酸没了，共价失效",
    "打法：转向强效非共价 / 变构抑制",
    "硬约束：对野生型 EGFR 必须高选择性",
    "现状：第四代 TKI 多在 I/II 期，尚无获批",
    "C797S 见于约 10%~25% 的奥希替尼耐药",
], AMBER)

# ΔΔG 示意
ax.add_patch(FancyBboxPatch((0.25, 0.35), 9.5, 6.2, boxstyle="round,pad=0.012,rounding_size=0.08",
                            linewidth=1.2, edgecolor="#DDE3E6", facecolor=CARD, zorder=1))
ax.text(0.6, 6.05, "选择性窗口怎么算：同一个分子，在突变体和野生型上的结合能之差", fontsize=11.5, fontweight="bold", color=INK)
ax.text(0.6, 5.55, "ΔΔG(选择性) = ΔG(野生型) - ΔG(突变体)　——　这个值越大，治疗窗口越宽", fontsize=10.5, color="#B7791F")
names = ["ΔG 突变体\n（要强）", "ΔG 野生型\n（要弱）"]
vals = [-12.5, -8.2]
bars = ax.bar([0.9, 2.1], vals, width=0.55, color=[GREEN, "#C7D2D6"], zorder=3)
for b, v in zip(bars, vals):
    ax.text(b.get_x() + b.get_width() / 2, v + 0.5, "%.1f" % v, ha="center", fontsize=10.5, fontweight="bold", color=INK)
ax.annotate("", xy=(2.1, -8.2), xytext=(2.1, -12.5), arrowprops=dict(arrowstyle="<->", color=RED, lw=1.8))
ax.text(2.35, -10.3, "ΔΔG = -8.2 -(-12.5)\n= +4.3 kcal/mol\n（越大越好）", fontsize=10.5, color=RED, va="center")
ax.set_xlim(0, 10); ax.set_ylim(-14.5, 12)
ax.text(4.6, -10.5, "FEP 最适合算的，就是这类「小基团局部微扰」带来的 ΔΔG；\n骨架跃迁请交给对接 + MM-GBSA，别迷信 FEP。", fontsize=10.5, color="#3F3F3F", va="center")
ax.axis("off")
fig.savefig(OUT + "fig5_two_wars.png", facecolor=BG, bbox_inches="tight", pad_inches=0.14)
plt.close(fig)

# ===== 图 6：学习路径 =====
mods = [("P1", "疾病机制与靶点", "①②", GREEN),
        ("P2", "结构生物学与 AI", "③", GREEN),
        ("P3", "分子生成与虚拟筛选", "④", BLUE),
        ("P4", "先导化合物优化", "⑤", BLUE),
        ("P5", "成药性建模与评估", "⑥", AMBER),
        ("P6", "实验设计与闭环迭代", "⑦", AMBER),
        ("P7", "端到端项目实战", "①–⑦", RED)]
fig, ax = plt.subplots(figsize=(7.2, 6.4), dpi=150)
fig.patch.set_facecolor(BG)
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
ax.text(0.2, 9.6, "后续学习路径：七个模块，按流水线顺序展开", fontsize=15.5, fontweight="bold", color=INK)
ax.text(0.2, 9.15, "每个模块统一节奏：原理回顾 → 工具实操 → 案例复现 → 开放练习", fontsize=10.5, color=GREEN)
y = 8.5
for i, (code, name, link, accent) in enumerate(mods):
    ax.add_patch(FancyBboxPatch((0.9, y - 0.78), 8.2, 0.72, boxstyle="round,pad=0.01,rounding_size=0.07",
                                linewidth=1.2, edgecolor="#DDE3E6", facecolor=CARD, zorder=1))
    ax.add_patch(Rectangle((0.9, y - 0.78), 0.08, 0.72, facecolor=accent, zorder=2))
    ax.text(1.28, y - 0.42, code, fontsize=12.5, fontweight="bold", color=accent, va="center")
    ax.text(2.05, y - 0.42, name, fontsize=12.5, fontweight="bold", color=INK, va="center")
    ax.text(8.75, y - 0.42, "对应环节 " + link, fontsize=10.5, color=GREY, va="center", ha="right")
    if i < len(mods) - 1:
        ax.annotate("", xy=(5.0, y - 0.84), xytext=(5.0, y - 0.74),
                    arrowprops=dict(arrowstyle="-|>", color="#C7D2D6", lw=1.5))
    y -= 1.02
ax.text(0.2, 0.75, "建议：先跟着 P1 走完一遍，再回头补齐自己最弱的模块", fontsize=10.5, color="#3F3F3F")
ax.text(0.2, 0.35, "作者自制 · 模块划分参考本系列课程安排", fontsize=9, color="#B8B8B8")
fig.savefig(OUT + "fig6_roadmap.png", facecolor=BG, bbox_inches="tight", pad_inches=0.14)
plt.close(fig)
print("done v2")
