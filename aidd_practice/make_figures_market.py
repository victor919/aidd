# -*- coding: utf-8 -*-
# 文章关联：AI制药行业调研篇_公众号文章.html
# 脚本职责：配图：fig1–fig5 市场、阶段、命中率、技术与巨头
# 使用注解：图表数据写在代码内，来自文章研究素材；重绘不代表重新获取最新行业数据。
# 完整对应关系与运行顺序：README.md
"""AI制药行业调研篇 配图生成（数据均来自《AI制药行业研究报告 2026》正文引用）"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle

OUT = "/Users/dujiepeng/WorkBuddy/Work/AI制药行业调研篇_配图"
os.makedirs(OUT, exist_ok=True)

GREEN = "#0FA87E"
GREEN_L = "#A8E5D3"
BLUE = "#2B5F8F"
BLUE_L = "#9FC0DE"
AMBER = "#F5A623"
AMBER_L = "#FBE2B7"
RED = "#E53E3E"
GRAY = "#9A9A9A"
INK = "#1A1A1A"
SUB = "#5D6C80"

plt.rcParams["font.sans-serif"] = ["Hiragino Sans GB"]
plt.rcParams["axes.unicode_minus"] = False


def style_ax(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#D8DEE5")
    ax.tick_params(colors=SUB, labelsize=10.5)


# ---------------- 图1：市场规模的两种口径 ----------------
def fig1():
    fig, ax = plt.subplots(figsize=(9.0, 5.2))
    labels = [
        "艾媒 2025 估算\n(中国窄口径市场研究)",
        "头豹 2026E 预测\n(模型外推)",
        "晶泰 2025\n药物发现分部",
        "英矽智能 2025\n全年(美元折算)",
        "剂泰科技 2025\n全年",
    ]
    values = [6.23, 20.4, 5.379, 4.05, 1.05]
    colors = [AMBER, AMBER_L, GREEN, GREEN, GREEN]
    bars = ax.bar(range(5), values, color=colors, width=0.62,
                  edgecolor="#ffffff", linewidth=1.2)
    for b, v in zip(bars, values):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.55, f"{v:.2f}",
                ha="center", va="bottom", fontsize=12, fontweight="bold",
                color=INK)
    ax.set_xticks(range(5))
    ax.set_xticklabels(labels, fontsize=9.8, color=SUB, linespacing=1.5)
    ax.set_ylabel("亿元人民币", fontsize=11, color=SUB)
    ax.set_ylim(0, 25.5)
    ax.set_title("同一个行业，两套口径：估算/预测 vs 公司财报已实现收入",
                 fontsize=14, color=INK, pad=14, fontweight="bold")
    # 图例（右上角空白区）
    ax.add_patch(Rectangle((0.60, 0.885), 0.045, 0.055, transform=ax.transAxes,
                           facecolor=AMBER, edgecolor="none"))
    ax.text(0.665, 0.912, "市场研究估算 / 预测（尚未发生）", transform=ax.transAxes,
            fontsize=10.2, color=SUB, va="center")
    ax.add_patch(Rectangle((0.60, 0.785), 0.045, 0.055, transform=ax.transAxes,
                           facecolor=GREEN, edgecolor="none"))
    ax.text(0.665, 0.812, "企业已实现收入（年报 / 中报披露）", transform=ax.transAxes,
            fontsize=10.2, color=SUB, va="center")
    ax.text(0.0, -0.24,
            "注：英矽智能 2025 年收入 5,623.9 万美元，按 1 美元 = 7.2 元人民币折算；三家公司业务边界与地域口径各不相同，\n"
            "不能简单相加代表行业总规模。仅晶泰一家公司的药物发现分部收入，就已接近艾媒口径下「全行业」估算值的 86%。",
            transform=ax.transAxes, fontsize=9.3, color=GRAY, ha="left",
            va="top", linespacing=1.7)
    style_ax(ax)
    fig.tight_layout(rect=[0, 0.10, 1, 1])
    fig.savefig(f"{OUT}/fig1_market_gap.png", dpi=150, facecolor="white")
    plt.close(fig)


# ---------------- 图2：设计被压缩，验证占比上升 ----------------
def fig2():
    fig, ax = plt.subplots(figsize=(9.0, 4.6))
    rows = [
        ("传统流程", 4.5, [("靶点到 PCC\n3-6 年", 4.5, GREEN, 10.2),
                          ("临床前\n1-2 年", 1.5, BLUE, 9.0),
                          ("I 期\n约 1 年", 1.0, "#4A7CFE", 8.6),
                          ("II 期\n约 2 年", 2.0, "#7FA6E8", 9.2),
                          ("III 期\n3-4 年", 3.5, BLUE_L, 9.6)]),
        ("AI 压缩设计后", 1.5, [("约\n1.5 年", 1.5, GREEN, 8.8),
                               ("临床前\n1-2 年", 1.5, BLUE, 9.0),
                               ("I 期\n约 1 年", 1.0, "#4A7CFE", 8.6),
                               ("II 期\n约 2 年", 2.0, "#7FA6E8", 9.2),
                               ("III 期\n3-4 年", 3.5, BLUE_L, 9.6)]),
    ]
    for idx, (name, design, segs) in enumerate(rows):
        y = 1.0 - idx * 1.0
        x = 0.0
        for label, w, c, fs in segs:
            ax.add_patch(Rectangle((x, y - 0.30), w, 0.60, facecolor=c,
                                   edgecolor="white", linewidth=2))
            ax.text(x + w / 2, y, label, ha="center", va="center",
                    fontsize=fs, color="white", fontweight="bold",
                    linespacing=1.35)
            x += w
        total = x
        ax.text(-0.28, y, name, ha="right", va="center", fontsize=11.5,
                color=INK, fontweight="bold")
        ax.text(x + 0.25, y,
                f"总计约 {total:g} 年\n验证环节占 {(total - design) / total * 100:.0f}%",
                ha="left", va="center", fontsize=10.6, color=SUB,
                fontweight="bold", linespacing=1.6)
    ax.set_xlim(-3.2, 18.6)
    ax.set_ylim(-0.85, 1.85)
    ax.axis("off")
    ax.set_title("AI 压缩的是前端设计，后端验证的时间不会自动缩短",
                 fontsize=14, color=INK, pad=12, fontweight="bold")
    ax.text(-3.2, -0.62,
            "注：阶段时长为行业常见区间的比例示意（取区间中位数），不是单一项目预测；实际周期受模态、适应症、终点、\n"
            "招募与监管沟通影响。绿色为设计与发现环节，蓝色为实验与临床验证环节。",
            fontsize=9.3, color=GRAY, ha="left", va="center", linespacing=1.7)
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    fig.savefig(f"{OUT}/fig2_phase_shift.png", dpi=150, facecolor="white")
    plt.close(fig)


# ---------------- 图3：Claude 蛋白设计命中率 ----------------
def fig3():
    fig, ax = plt.subplots(figsize=(9.0, 4.4))
    total, hit = 1320, 354
    miss = total - hit
    ax.add_patch(Rectangle((0, 1.05), hit / total * 100, 0.60, facecolor=AMBER,
                           edgecolor="white", linewidth=2))
    ax.add_patch(Rectangle((hit / total * 100, 1.05), miss / total * 100, 0.60,
                           facecolor="#E5E9EE", edgecolor="white", linewidth=2))
    ax.text(hit / total * 44, 1.35, f"实验确认结合 {hit} 个\n26.8%",
            ha="center", va="center", fontsize=12, color="white",
            fontweight="bold", linespacing=1.4)
    ax.text(hit / total * 100 + miss / total * 50, 1.35,
            f"未确认结合 {miss} 个\n73.2%", ha="center", va="center",
            fontsize=12, color=SUB, fontweight="bold", linespacing=1.4)
    ax.add_patch(Rectangle((10, 0.15), 5, 0.45, facecolor=GREEN_L,
                           edgecolor=GREEN, linewidth=1.2, hatch="//"))
    ax.text(12.5, -0.05, "常见水平\n10% - 15%", ha="center", va="top",
            fontsize=9.6, color=GREEN, fontweight="bold", linespacing=1.4)
    ax.text(21, 0.375, "本案例整体命中率 26.8%；不同运行设置在 22% - 35% 之间波动",
            fontsize=10.4, color=SUB, va="center")
    ax.plot([26.8, 26.8], [0.72, 1.00], color=AMBER, lw=2, ls="--")
    ax.set_xlim(0, 100)
    ax.set_ylim(-0.55, 1.80)
    ax.set_xticks([0, 20, 40, 60, 80, 100])
    ax.set_xticklabels(["0", "20%", "40%", "60%", "80%", "100%"],
                       fontsize=10.2, color=GRAY)
    ax.tick_params(axis="x", length=3, pad=4)
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title("1,320 个 AI 设计的蛋白结合体，354 个被湿实验确认",
                 fontsize=14, color=INK, pad=12, fontweight="bold")
    fig.text(0.04, 0.035,
             "数据来源：Anthropic 蛋白设计项目，由 Adaptyv Bio 与 Twist Bioscience 完成外部湿实验验证。\n"
             "「结合」不等于治疗有效，更不等于临床成功率。",
             fontsize=9.3, color=GRAY, ha="left", va="bottom", linespacing=1.7)
    fig.tight_layout(rect=[0, 0.16, 1, 1])
    fig.savefig(f"{OUT}/fig3_hit_rate.png", dpi=150, facecolor="white")
    plt.close(fig)


# ---------------- 图4：READ-DMTA 框架 ----------------
def fig4():
    fig, ax = plt.subplots(figsize=(9.2, 5.2))
    boxes = [
        ("READ", "读取疾病", "这个病发生了什么？\n该干预哪个机制？",
         "组学 / 病理 / 高内涵成像\n因果推断与靶点排序",
         "Recursion、10x Genomics\nOwkin、英矽 PandaOmics", GREEN),
        ("DESIGN", "设计", "下一批最值得制造和\n测试的分子是什么？",
         "结构预测、虚拟筛选\n生成模型、FEP、MPO",
         "Isomorphic、Schrodinger\n英矽、晶泰、深势", BLUE),
        ("MAKE", "制造", "数字设计能变成\n真实样品吗？",
         "基因 / 小分子合成\n蛋白表达纯化、自动化",
         "Twist、金斯瑞\nChemify、晶泰", "#4A7CFE"),
        ("TEST", "测试", "样品结合吗？有效吗？\n安全、可成药吗？",
         "SPR / BLI、细胞表型\nDMPK、毒理、可开发性",
         "Carterra、金斯瑞\n专业 CRO", AMBER),
        ("ANALYZE", "分析与学习", "本轮学到了什么？\n下一轮测什么？",
         "ELN / LIMS、主动学习\n贝叶斯优化、失败归档",
         "端到端 AIDD 平台\n自驱实验室", "#7C5CBF"),
    ]
    w, gap = 1.74, 0.10
    for i, (en, cn, q, tech, comp, color) in enumerate(boxes):
        x = i * (w + gap)
        ax.add_patch(Rectangle((x, 1.55), w, 2.35, facecolor="white",
                               edgecolor=color, linewidth=1.8))
        ax.add_patch(Rectangle((x + 0.09, 3.60), w - 0.18, 0.27,
                               facecolor=color, edgecolor="none"))
        ax.text(x + w / 2, 3.735, en, ha="center", va="center", fontsize=11.2,
                color="white", fontweight="bold")
        ax.text(x + w / 2, 3.30, cn, ha="center", va="center", fontsize=12.4,
                color=INK, fontweight="bold")
        ax.text(x + w / 2, 2.85, q, ha="center", va="center", fontsize=9.2,
                color=SUB, linespacing=1.55)
        ax.text(x + w / 2, 2.20, tech, ha="center", va="center", fontsize=8.8,
                color=color, linespacing=1.6, fontweight="bold")
        ax.text(x + w / 2, 1.80, comp, ha="center", va="center", fontsize=8.4,
                color=GRAY, linespacing=1.6)
        if i < len(boxes) - 1:
            ax.annotate("", xy=(x + w + gap + 0.01, 2.72),
                        xytext=(x + w - 0.01, 2.72),
                        arrowprops=dict(arrowstyle="-|>", color="#9FB0C2", lw=1.8))
    ax.annotate("", xy=(0.05, 1.26), xytext=(8.90, 1.26),
                arrowprops=dict(arrowstyle="-|>", color=GREEN, lw=2.0,
                                connectionstyle="arc3,rad=-0.16"))
    ax.text(4.47, 0.85, "正负结果回流，模型更新，进入下一轮 DMTA 循环",
            ha="center", va="center", fontsize=10.6, color=GREEN,
            fontweight="bold")
    ax.set_xlim(-0.30, 9.25)
    ax.set_ylim(0.55, 4.10)
    ax.axis("off")
    ax.set_title("READ - DMTA：把一家 AI 制药公司放回真实研发工作流里看",
                 fontsize=14, color=INK, pad=10, fontweight="bold")
    fig.tight_layout(rect=[0, 0.02, 1, 1])
    fig.savefig(f"{OUT}/fig4_read_dmta.png", dpi=150, facecolor="white")
    plt.close(fig)


# ---------------- 图5：科技巨头三层版图 ----------------
def fig5():
    fig, ax = plt.subplots(figsize=(9.6, 5.5))
    layers = [
        ("上层 · 候选药资产与科研智能体（直接承担临床成败）", 3.10,
         [("Alphabet / Isomorphic Labs", "自有药物设计引擎", GREEN),
          ("OpenAI / Anthropic",
           "不自持管线，做工具编排", "#7C5CBF"),
          ("字节跳动 / Anew Labs",
           "2026 分拆，4 条管线临床前", AMBER)]),
        ("中层 · 生物基础模型与科学工具链", 2.00,
         [("Microsoft Discovery", "科研智能体平台", BLUE),
          ("NVIDIA BioNeMo", "训练与推理微服务", BLUE),
          ("Meta ESM 系列", "开放研究为主", GRAY),
          ("华为云盘古药物分子", "分子预训练模型", AMBER)]),
        ("底层 · GPU、云与 HPC（收益与算力采用率相关）", 0.90,
         [("NVIDIA", "算力与加速库", "#4A7CFE"),
          ("AWS HealthOmics", "组学数据底座", "#4A7CFE"),
          ("Microsoft Azure", "企业云与 HPC", "#4A7CFE"),
          ("华为 昇腾 / 鲲鹏", "国产算力框架", AMBER)]),
    ]
    for title, y, items in layers:
        ax.text(0.0, y + 0.60, title, fontsize=11.4, color=INK,
                fontweight="bold", va="center")
        n = len(items)
        w = 3.00 if n == 3 else (2.30 if n == 4 else 4.47)
        gap = 0.08 if n == 4 else 0.15
        for i, (name, desc, color) in enumerate(items):
            x = i * (w + gap)
            ax.add_patch(FancyBboxPatch((x, y - 0.25), w, 0.72,
                                        boxstyle="round,pad=0.01,rounding_size=0.08",
                                        facecolor="white", edgecolor=color,
                                        linewidth=1.7))
            ax.add_patch(Rectangle((x, y - 0.25), 0.055, 0.72, facecolor=color,
                                   edgecolor="none"))
            ax.text(x + 0.17, y + 0.16, name, fontsize=9.8, color=INK,
                    fontweight="bold", va="center")
            ax.text(x + 0.17, y - 0.12, desc, fontsize=8.8, color=SUB, va="center")
    ax.text(0.0, 0.20,
            "Tesla：公开 AI 布局集中于自动驾驶与机器人，未发现 AIDD 直接布局；Apple 的 SimpleFold 属研究与端侧计算信号，不是药企级平台。\n"
            "国内两条相反的路：字节跳动把 AI 制药分拆为 Anew Labs 自持管线（2026 年 9 月首轮 2.9 亿美元、估值约 15 亿美元，字节持股 56%），\n"
            "但 4 条管线全部处于临床前；华为则明确不做医疗产品、不碰管线，以昇腾/鲲鹏算力 + 盘古药物分子大模型 + 伙伴生态做底座。",
            fontsize=9.4, color=RED, va="center", linespacing=1.7)
    ax.set_xlim(-0.35, 9.30)
    ax.set_ylim(-0.45, 3.95)
    ax.axis("off")
    ax.set_title("科技巨头在 AIDD 的三层仓位：不是七家同质化的 AI 药企",
                 fontsize=14, color=INK, pad=10, fontweight="bold")
    fig.tight_layout(rect=[0, 0.02, 1, 1])
    fig.savefig(f"{OUT}/fig5_bigtech.png", dpi=150, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    fig1(); fig2(); fig3(); fig4(); fig5()
    for f in sorted(os.listdir(OUT)):
        p = os.path.join(OUT, f)
        print(f, f"{os.path.getsize(p)/1024:.0f} KB")
