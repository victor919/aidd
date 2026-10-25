# -*- coding: utf-8 -*-
# 文章关联：AI制药行业调研篇_公众号文章.html
# 脚本职责：配图：fig9 市场预测口径对照
# 使用注解：数据为代码内的历史研究摘要；当前网页引用 fig9_global_forecast.png。
# 完整对应关系与运行顺序：README.md
"""图：全球 AIDD 市场规模——五家机构，五个终点（口径校准主线）"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.sans-serif"] = ["Hiragino Sans GB"]
plt.rcParams["axes.unicode_minus"] = False

OUT = "/Users/dujiepeng/WorkBuddy/Work/AI制药行业调研篇_配图"
INK, GRAY = "#1A1A1A", "#9A9A9A"
GREEN, GREEN_D = "#0FA87E", "#0B7A5B"
BLUE, AMBER, LINE = "#4A7CFE", "#F5A623", "#D8DEE6"


def fig9():
    fig, ax = plt.subplots(figsize=(9.0, 5.5))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    rows = [
        ("QYResearch",             2024, 22.8, 2031, 144.7, "30.7%", False),
        ("Arizton",                2024, 17.1, 2030,  85.2, "30.6%", False),
        ("Research and Markets",   2025, 23.5, 2033, 137.7, "24.8%", False),
        ("中信建投",                2025, 24.9, 2035, 460.0, "33%+",  False),
        ("Expert Market Research", 2025, 25.5, None, None,   "28.7%", True),
    ]

    xmax, y0, dy = 480.0, 3.95, 0.92

    for i, (name, by, bv, ty, tv, cagr, missing) in enumerate(rows):
        y = y0 - i * dy
        end = xmax if missing else tv
        ax.plot([bv, end], [y, y], color=LINE, lw=3.0, zorder=2,
                linestyle="--" if missing else "-",
                solid_capstyle="round")

        # 起点
        ax.plot([bv], [y], marker="o", ms=8.5, color=GREEN, zorder=5,
                markeredgecolor="white", markeredgewidth=1.5)
        ax.text(bv, y + 0.26, f"{bv:g} 亿$", fontsize=9.5, color=GREEN_D,
                ha="center", va="bottom", fontweight="bold")

        # 终点
        if missing:
            ax.text(end + 8, y, f"未披露终值（仅 CAGR {cagr}）", fontsize=9.2,
                    color=GRAY, ha="left", va="center", style="italic")
        else:
            col = AMBER if tv >= 300 else BLUE
            ax.plot([tv], [y], marker="o", ms=12.0, color=col, zorder=5,
                    markeredgecolor="white", markeredgewidth=1.7)
            if tv < 130:
                ax.text(tv + 16, y + 0.29, f"{tv:g} 亿$ / {ty} 年", fontsize=9.7,
                        color=INK, ha="left", va="bottom", fontweight="bold")
            else:
                ax.text(tv, y + 0.29, f"{tv:g} 亿$ / {ty} 年", fontsize=9.7,
                        color=INK, ha="center", va="bottom", fontweight="bold")

        # 机构名
        ax.text(-172, y + 0.06, name, fontsize=10.6, color=INK, ha="left",
                va="center", fontweight="bold")
        ax.text(-172, y - 0.26, f"基准 {by} 年 · CAGR {cagr}", fontsize=8.7,
                color=GRAY, ha="left", va="center")

    # 起点几乎一致的提示
    ax.axvspan(15.0, 27.0, color=GREEN, alpha=0.10, zorder=1)
    ax.annotate("起点几乎一致：17-26 亿$", xy=(21, 0.02), xytext=(21, -0.42),
                fontsize=9.5, color=GREEN_D, ha="center", va="top",
                fontweight="bold")

    ax.set_xlim(-175, 495)
    ax.set_ylim(-0.62, y0 + 0.72)
    ax.set_xticks([0, 100, 200, 300, 400, 480])
    ax.set_xticklabels(["0", "100", "200", "300", "400", "480 亿$"],
                       fontsize=9.3, color=GRAY)
    ax.set_yticks([])
    ax.tick_params(axis="x", length=0, pad=6)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color("#E5EAEF")
    ax.spines["bottom"].set_bounds(0, 480)
    ax.grid(axis="x", color="#F1F4F7", lw=1.0)
    ax.set_axisbelow(True)

    fig.text(0.012, 0.975, "全球 AI 制药市场规模：起点几乎一致，终点相差 5 倍",
             fontsize=14.5, color=INK, ha="left", va="top", fontweight="bold")
    fig.text(0.012, 0.895,
             "五家机构口径并不相同：QYResearch 与 Arizton 更聚焦 AIDD 软件与平台市场；Research and Markets 与\n"
             "Expert Market Research 口径更宽，包含 AI 在药物发现全流程的服务与工具。这张图只用于观察分歧幅度，\n"
             "不能把任何一条线当作「行业实际规模」。",
             fontsize=8.9, color=GRAY, ha="left", va="top", linespacing=1.7)

    fig.tight_layout(rect=[0, 0.02, 1, 0.845])
    p = f"{OUT}/fig9_global_forecast.png"
    fig.savefig(p, dpi=190, facecolor="white")
    plt.close(fig)
    print("saved", p)


if __name__ == "__main__":
    fig9()
