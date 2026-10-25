# -*- coding: utf-8 -*-
# 文章关联：AI制药行业调研篇_公众号文章.html
# 脚本职责：配图：来源重绘 fig6 / fig7 / fig8
# 使用注解：当前文章引用 fig7；fig6 / fig8 为保留素材，不能把所有生成图都算作当前正文插图。
# 完整对应关系与运行顺序：README.md
"""来源图核对后的高清重绘：
fig6 传统制药 vs AI 制药（对应来源1 上半 / 3个Right）
fig7 AIDD 技术路线总览 + 干湿闭环治理层（对应来源4）
fig8 三条验证主线 + 商业化正向飞轮 + 关键数据（对应来源1 下半）
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch

OUT = "/Users/dujiepeng/WorkBuddy/Work/AI制药行业调研篇_配图"
os.makedirs(OUT, exist_ok=True)

GREEN = "#0FA87E"
GREEN_D = "#0B6E52"
GREEN_L = "#E6F7F1"
BLUE = "#2B5F8F"
BLUE_L = "#EAF1F8"
PURPLE = "#7C5CBF"
PURPLE_L = "#F0EBFA"
AMBER = "#F5A623"
AMBER_L = "#FDF1DC"
ORANGE = "#E8730C"
GRAY = "#6B7F94"
GRAY_L = "#F0F3F6"
INK = "#1A1A1A"
SUB = "#4A5B70"
LINE = "#D8DEE5"

plt.rcParams["font.sans-serif"] = ["Hiragino Sans GB"]
plt.rcParams["axes.unicode_minus"] = False


def blank(figsize):
    fig, ax = plt.subplots(figsize=figsize, dpi=150)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")
    return fig, ax


def box(ax, x, y, w, h, fc, ec, lw=1.6, r=0.9, z=2):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle=f"round,pad=0.2,rounding_size={r}",
                                facecolor=fc, edgecolor=ec, linewidth=lw, zorder=z))


def pill(ax, x, y, w, h, fc, ec, r=1.4):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle=f"round,pad=0.15,rounding_size={r}",
                                facecolor=fc, edgecolor=ec, linewidth=1.3, zorder=3))


def arrow(ax, x1, y1, x2, y2, color="#9FB0C2", lw=1.8):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                 mutation_scale=13, color=color, lw=lw, zorder=4))


# ============ fig6：传统制药 vs AI 制药（核心对比 + 3 个 Right） ============
def fig6():
    fig, ax = blank((9.4, 6.4))
    ax.text(50, 95.5, "传统制药 vs AI 制药", fontsize=17, color=INK,
            fontweight="bold", ha="center", va="center")
    ax.text(50, 89.5,
            "从「纸上谈兵」到「真金白银」：AI 制药进入商业与临床双重验证期",
            fontsize=10.6, color=SUB, ha="center", va="center")

    # 表头
    box(ax, 1.5, 78.0, 19, 8.0, GRAY_L, LINE, lw=1.2)
    ax.text(11, 82.0, "对比维度", fontsize=10.4, color=SUB, ha="center",
            va="center", fontweight="bold")
    box(ax, 21.5, 78.0, 38.5, 8.0, "#E9EDF2", LINE, lw=1.2)
    ax.text(40.75, 82.0, "传统制药 ｜ 试错法", fontsize=11.6, color="#3C5266",
            ha="center", va="center", fontweight="bold")
    box(ax, 61.0, 78.0, 37.5, 8.0, GREEN, GREEN, lw=1.2)
    ax.text(79.75, 82.0, "AI 制药 ｜ 预测法", fontsize=11.6, color="white",
            ha="center", va="center", fontweight="bold")

    rows = [
        ("核心逻辑", "物理筛选 ＋ 专家经验", "数据预测 ＋ 虚拟筛选"),
        ("候选药物提名周期", "传统：约 2.5 - 4 年", "AI：约 12 - 18 个月"),
        ("主要价值", "实验验证与全流程开发", "提升靶点发现与分子设计效率"),
        ("不可替代环节", "实验、临床与监管审批", "AI 仍需接受人体试验硬考核"),
    ]
    y = 68.5
    for i, (k, v1, v2) in enumerate(rows):
        fc = "#FFFFFF" if i % 2 == 0 else "#FAFBFC"
        box(ax, 1.5, y, 19, 9.6, GRAY_L, LINE, lw=1.1)
        ax.text(11, y + 4.8, k, fontsize=9.8, color=SUB, ha="center",
                va="center", fontweight="bold", linespacing=1.4)
        box(ax, 21.5, y, 38.5, 9.6, fc, LINE, lw=1.1)
        ax.text(40.75, y + 4.8, v1, fontsize=10.4, color=INK, ha="center",
                va="center", linespacing=1.5)
        box(ax, 61.0, y, 37.5, 9.6, GREEN_L if i % 2 == 0 else "#F3FBF8",
            "#BEE8DA", lw=1.1)
        ax.text(79.75, y + 4.8, v2, fontsize=10.4, color=GREEN_D,
                ha="center", va="center", linespacing=1.5, fontweight="bold")
        y -= 10.4

    # 3 个 Right
    ax.text(1.5, 26.5, "药物研发的 3 个 Right", fontsize=11.8, color=INK,
            fontweight="bold", va="center")
    items = [
        ("1", "Right Target", "正确靶点", "AI 系统性提升", BLUE, BLUE_L),
        ("2", "Right Molecule", "正确分子", "AI 系统性提升", PURPLE, PURPLE_L),
        ("3", "Right Clinical Trial", "正确临床试验", "临床硬考核", ORANGE, AMBER_L),
    ]
    x = 1.5
    for num, en, cn, tag, c, cl in items:
        box(ax, x, 10.5, 31.0, 13.5, cl, c, lw=1.6)
        ax.add_patch(plt.Circle((x + 3.6, 21.0), 1.7, facecolor=c,
                                edgecolor="none", zorder=4))
        ax.text(x + 3.6, 21.0, num, fontsize=10.5, color="white",
                ha="center", va="center", fontweight="bold", zorder=5)
        ax.text(x + 6.4, 21.4, en, fontsize=9.6, color=c, va="center",
                fontweight="bold")
        ax.text(x + 6.4, 18.6, cn, fontsize=10.6, color=INK, va="center",
                fontweight="bold")
        ax.text(x + 15.5, 13.6, tag, fontsize=9.6, color=c, va="center",
                ha="center", fontweight="bold")
        if x < 60:
            arrow(ax, x + 31.4, 17.2, x + 33.4, 17.2, color="#B7C2CE", lw=2.0)
        x += 33.0

    ax.add_patch(Rectangle((0.5, 1.5), 99, 6.6, facecolor=BLUE,
                           edgecolor="none", zorder=2))
    ax.text(50, 4.8, "AI 无法替代实验、临床试验和监管审批",
            fontsize=12.0, color="white", ha="center", va="center",
            fontweight="bold", zorder=5)

    fig.tight_layout(pad=0.4)
    fig.savefig(f"{OUT}/fig6_traditional_vs_ai.png", dpi=150, facecolor="white")
    plt.close(fig)


# ============ fig7：AIDD 技术路线 + 干湿闭环治理层 ============
def fig7():
    fig, ax = blank((9.2, 4.6))
    ax.text(0.5, 93, "AIDD 目前技术路线", fontsize=15, color=INK,
            fontweight="bold", va="center")

    steps = [
        ("数据与知识", "化学结构\n组学 文献\n实验 临床 RWD"),
        ("靶点发现", "NLP 知识图谱\n图学习\n因果推断"),
        ("结构与筛选", "结构预测 对接\n虚拟筛选\n分子动力学"),
        ("生成与优化", "扩散 流匹配\nVAE GAN\n强化学习"),
        ("成药性与合成", "QSAR ADMET\n逆合成\n可合成性"),
        ("实验与临床", "DMTA\n自动化实验\n患者分层"),
    ]
    n = len(steps)
    w, gap = 15.1, 1.0
    x0 = 0.5
    for i, (title, body) in enumerate(steps):
        x = x0 + i * (w + gap)
        box(ax, x, 52, w, 30, "#FBFDFF", "#C9D8E8", lw=1.4)
        ax.text(x + w / 2, 76.5, title, fontsize=10.6, color=BLUE,
                ha="center", va="center", fontweight="bold")
        ax.text(x + w / 2, 63.5, body, fontsize=8.6, color=SUB,
                ha="center", va="center", linespacing=1.7)
        if i < n - 1:
            arrow(ax, x + w + 0.1, 67, x + w + gap - 0.1, 67,
                  color="#B7C2CE", lw=1.6)

    # 干湿闭环与治理层
    box(ax, 0.5, 6, 99, 36, GREEN_L, GREEN, lw=1.8)
    ax.text(50, 36, "干湿闭环与治理层", fontsize=11.6, color=GREEN_D,
            ha="center", va="center", fontweight="bold")
    nodes = ["模型提出候选", "湿实验验证", "数据回流", "模型更新",
             "人工复核与变更控制"]
    nx, nw, ngap = 3.0, 16.6, 3.3
    for i, nd in enumerate(nodes):
        x = nx + i * (nw + ngap)
        pill(ax, x, 15, nw, 11, "white", GREEN, r=1.2)
        ax.text(x + nw / 2, 20.5, nd, fontsize=9.4, color=GREEN_D,
                ha="center", va="center", fontweight="bold")
        if i < len(nodes) - 1:
            arrow(ax, x + nw + 0.3, 20.5, x + nw + ngap - 0.3, 20.5,
                  color=GREEN, lw=1.8)

    fig.tight_layout(pad=0.4)
    fig.savefig(f"{OUT}/fig7_tech_roadmap.png", dpi=150, facecolor="white")
    plt.close(fig)


# ============ fig8：三条验证主线 + 商业飞轮 + 关键数据 ============
def fig8():
    fig, ax = blank((9.2, 5.4))
    ax.text(0.5, 96, "AI 制药的三条验证主线与商业化飞轮", fontsize=14,
            color=INK, fontweight="bold", va="center")

    # 三条验证主线
    ax.text(0.5, 88, "从技术到商业再到临床的价值验证主线", fontsize=10.6,
            color=INK, va="center", fontweight="bold")
    lines = [
        ("1", "技术验证", "模型能否找到更优候选分子", BLUE, BLUE_L),
        ("2", "商业验证", "药企是否签约、付款并持续复购", PURPLE, PURPLE_L),
        ("3", "临床验证", "安全性、有效性与真实患者获益", ORANGE, AMBER_L),
    ]
    x = 0.5
    for num, t, d, c, cl in lines:
        box(ax, x, 70, 31.2, 13, cl, c, lw=1.6)
        ax.add_patch(plt.Circle((x + 3.4, 80.5), 1.9, facecolor=c,
                                edgecolor="none", zorder=4))
        ax.text(x + 3.4, 80.5, num, fontsize=10.4, color="white", ha="center",
                va="center", fontweight="bold", zorder=5)
        ax.text(x + 6.4, 80.6, t, fontsize=10.8, color=INK, va="center",
                fontweight="bold")
        ax.text(x + 6.4, 74.5, d, fontsize=9.0, color=SUB, va="center")
        if x < 60:
            arrow(ax, x + 31.6, 76.5, x + 33.6, 76.5, color="#B7C2CE", lw=2.0)
        x += 33.4

    # 商业化正向飞轮
    ax.text(0.5, 62, "商业化正向飞轮：从候选资产到商业回报", fontsize=10.6,
            color=INK, va="center", fontweight="bold")
    fly = ["AI 平台输出候选资产", "药企负责临床与全球开发",
           "首付款 + 里程碑 + 授权", "复购与数据反馈"]
    fx, fw, fg = 0.5, 22.2, 3.4
    for i, t in enumerate(fly):
        x = fx + i * (fw + fg)
        pill(ax, x, 44, fw, 13, "#F7FBFA", GREEN, r=1.2)
        ax.text(x + fw / 2, 50.5, t, fontsize=8.9, color=GREEN_D, ha="center",
                va="center", fontweight="bold", linespacing=1.5)
        if i < len(fly) - 1:
            arrow(ax, x + fw + 0.4, 50.5, x + fw + fg - 0.4, 50.5,
                  color=GREEN, lw=1.8)
    rc = fx + 3 * (fw + fg) + fw / 2
    lc = fx + fw / 2
    ax.plot([rc, rc], [44, 42.2], color=GREEN, lw=1.6, zorder=3)
    ax.plot([lc, lc], [44, 42.2], color=GREEN, lw=1.6, zorder=3)
    ax.annotate("", xy=(lc, 42.2), xytext=(rc, 42.2),
                arrowprops=dict(arrowstyle="-|>", color=GREEN, lw=1.8,
                                mutation_scale=14))
    ax.text(50, 38.6, "数据与经验沉淀，持续提升 AI 平台能力", fontsize=9.4,
            color=GREEN_D, ha="center", va="center", fontweight="bold")

    # 关键数据条
    data = [
        ("2026 年上半年合作潜在总额", "近 70 亿美元", BLUE, BLUE_L),
        ("全球进入临床的 AI 药物管线", "170+ 条（推进至 Ⅲ 期 10+）", PURPLE, PURPLE_L),
        ("行业估算：上市时间缩短 20%-25%",
         "研发周期缩短约 3 年\n整体成本降低 25%-30%", ORANGE, AMBER_L),
    ]
    x = 0.5
    for t, v, c, cl in data:
        box(ax, x, 20, 31.2, 15, cl, c, lw=1.5)
        ax.text(x + 15.6, 30.6, t, fontsize=8.6, color=SUB, ha="center",
                va="center", fontweight="bold", linespacing=1.4)
        ax.text(x + 15.6, 24.4, v, fontsize=10.2, color=c, ha="center",
                va="center", fontweight="bold", linespacing=1.6)
        x += 33.4

    ax.add_patch(Rectangle((0.5, 6.5), 99, 9.5, facecolor="#FDF6E8",
                           edgecolor=AMBER, linewidth=1.3, zorder=2))
    ax.text(50, 11.2,
            "口径提醒：「170+ 条管线 / 10+ Ⅲ 期」与「20%-25%、约 3 年、25%-30%」转引自知乎文章（2026 年 7 月），"
            "属二级来源与行业估算，\n未在监管文件或企业公告中独立核实；潜在合同总额包含未来里程碑，不等于已实现收入。",
            fontsize=8.6, color="#8A5A00", ha="center", va="center",
            linespacing=1.7, zorder=5)

    fig.tight_layout(pad=0.4)
    fig.savefig(f"{OUT}/fig8_value_chain.png", dpi=150, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    fig6(); fig7(); fig8()
    for f in ["fig6_traditional_vs_ai.png", "fig7_tech_roadmap.png",
              "fig8_value_chain.png"]:
        p = os.path.join(OUT, f)
        print(f, f"{os.path.getsize(p)/1024:.0f} KB")
