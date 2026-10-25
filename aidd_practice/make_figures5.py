# -*- coding: utf-8 -*-
# 文章关联：AI制药实践篇_公众号文章.html
# 脚本职责：配图：v5：fig8 DMTA / LoRA 与 fig9 更新代价
# 使用注解：历史增量版本，不是五套完整图；同名图会被后运行的版本覆盖。数值为代码内实测摘要，重跑实验后须人工同步。
# 完整对应关系与运行顺序：README.md
"""实践篇 v5 配图：DMTA 闭环 + 增量微调实测（数据来自本机实跑 lora_loop.py）"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Circle
import matplotlib.gridspec as gridspec
import numpy as np

OUT = "/Users/dujiepeng/WorkBuddy/Work/AI制药实践篇_配图/"
plt.rcParams["font.family"] = "Hiragino Sans GB"
plt.rcParams["axes.unicode_minus"] = False
GREEN, BLUE, INK, GREY = "#0FA87E", "#2B5F8F", "#1A1A1A", "#9A9A9A"
BG, CARD, RED, AMBER = "#FFFFFF", "#F7F9FA", "#E53E3E", "#F5A623"

fig = plt.figure(figsize=(7.2, 8.2), dpi=150)
fig.patch.set_facecolor(BG)
gs = gridspec.GridSpec(2, 1, height_ratios=[1.05, 1.0], hspace=0.48)

# ============ 上半：DMTA 循环 ============
ax = fig.add_subplot(gs[0])
ax.set_xlim(0, 10); ax.set_ylim(0, 9.2); ax.axis("off")
ax.text(0.05, 8.85, "第七站的真实形态：DMTA 循环", fontsize=15.5, fontweight="bold", color=INK)
ax.text(0.05, 8.35, "干端算力再强，也要过湿端这一关；每转一圈，模型与假说都应更准", fontsize=10.5, color=GREY)

cx, cy = 4.55, 4.45
steps = [("Design", "设计", "AI 推荐 10-50 个分子", 90, 2.55, GREEN),
         ("Make", "合成", "自动化 / 外包合成", 0, 3.05, AMBER),
         ("Test", "测试", "生化 + 细胞 + ADMET", 270, 2.55, BLUE),
         ("Analyze", "分析", "归因：设计 / 模型 / 实验？", 180, 3.05, RED)]
for eng, zh, desc, ang, R, col in steps:
    a = np.deg2rad(ang)
    x, y = cx + R * np.cos(a), cy + R * np.sin(a)
    ax.add_patch(FancyBboxPatch((x - 1.36, y - 0.66), 2.72, 1.32,
                                boxstyle="round,pad=0.012,rounding_size=0.14",
                                linewidth=1.4, edgecolor=col, facecolor=CARD, zorder=2))
    ax.text(x, y + 0.36, "%s %s" % (eng, zh), fontsize=11.5, fontweight="bold", color=col,
            ha="center", va="center", zorder=3)
    ax.text(x, y - 0.20, desc, fontsize=9.0, color="#4A4A4A", ha="center", va="center", zorder=3)

Rarc = 2.75
for i, (_, _, _, ang, _, col) in enumerate(steps):
    a0, a1 = np.deg2rad(ang - 22), np.deg2rad(ang - 90 + 22)
    x0, y0 = cx + Rarc * np.cos(a0), cy + Rarc * np.sin(a0)
    x1, y1 = cx + Rarc * np.cos(a1), cy + Rarc * np.sin(a1)
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="-|>", color="#B9C4C8", lw=1.6,
                                connectionstyle="arc3,rad=-0.28"), zorder=1)

ax.text(cx, cy + 0.28, "转一圈 = 一轮迭代", fontsize=11, fontweight="bold", color=INK,
        ha="center", va="center", zorder=4)
ax.text(cx, cy - 0.30, "周期越短，迭代越多", fontsize=9.5, color=GREY,
        ha="center", va="center", zorder=4)
ax.text(cx, cy - 0.88, "回传的是带误差的真值", fontsize=9.5, color=RED,
        ha="center", va="center", zorder=4)

ax.text(0.05, 0.55, "瓶颈通常不在 AI：一圈的实际周期 = 排队合成 + 检测等待 + 数据整理，而计算只占几小时",
        fontsize=10, color=AMBER)
ax.text(0.05, 0.10, "数据来源：作者根据 industry DMTA 通用实践整理", fontsize=9, color="#B8B8B8")

# ============ 下半：增量微调实测 ============
n_list = [8, 20, 50, 120, 250]
ft_delta = [-0.002, 0.031, -0.020, -0.123, -0.208]      # 全参微调 新系列 RMSE 变化
lo_delta = [0.029, 0.029, 0.018, -0.074, -0.145]        # LoRA r=8 新系列 RMSE 变化

ax2 = fig.add_subplot(gs[1])
ax2.text(-0.11, 1.10, "实测：累计回传多少个分子，才值得更新模型？", fontsize=15.5,
         fontweight="bold", color=INK, transform=ax2.transAxes)
ax2.text(-0.11, 1.035, "纵轴 = 新系列 RMSE 相对「不更新」的变化（向下为改善）；3 个随机种子平均",
         fontsize=10, color=GREY, transform=ax2.transAxes)

ax2.axhline(0, color="#C7D2D6", lw=1.2)
ax2.plot(n_list, ft_delta, marker="o", lw=2.3, color=BLUE, markersize=6, label="全参微调（34.3 万参数）")
ax2.plot(n_list, lo_delta, marker="s", lw=2.3, color=GREEN, markersize=6, label="LoRA 适配器（1.4 万参数）")
ax2.axvspan(4, 62, color="#F5F7F9", zorder=0)
ax2.annotate("少于约 50 个数据时，两种做法\n都在噪声里打转（甚至更差）", xy=(22, 0.030),
             xytext=(58, 0.062), fontsize=9.8, color=RED,
             arrowprops=dict(arrowstyle="->", color=RED, lw=1.1))
ax2.annotate("累计到 120 个以上，收益才稳定出现", xy=(165, -0.150), xytext=(95, -0.255),
             fontsize=9.8, color=GREEN,
             arrowprops=dict(arrowstyle="->", color=GREEN, lw=1.1))
ax2.set_xlabel("湿实验累计回传的自家分子数", fontsize=11, color="#3F3F3F", labelpad=6)
ax2.set_ylabel("新系列 RMSE 变化量", fontsize=11, color="#3F3F3F")
ax2.set_ylim(-0.30, 0.12)
ax2.set_xlim(-8, 268)
ax2.tick_params(labelsize=10, colors=GREY)
ax2.legend(frameon=False, fontsize=10, loc="upper right")
for s in ["top", "right"]:
    ax2.spines[s].set_visible(False)
for s in ["left", "bottom"]:
    ax2.spines[s].set_color("#E0E3E6")

ax2.text(-0.11, -0.30, "数据：TDC Lipophilicity（4200 真实分子，logD7.4），骨架切分出新系列；PyTorch 本机实跑，3 个随机种子平均",
         fontsize=9, color="#B8B8B8", transform=ax2.transAxes)
ax2.text(-0.11, -0.38, "教训：拿到十几个数据就急着微调，多半是在给模型灌噪声",
         fontsize=10, color=RED, fontweight="bold", transform=ax2.transAxes)

fig.savefig(OUT + "fig8_dmta_lora.png", facecolor=BG, bbox_inches="tight", pad_inches=0.16)
plt.close(fig)

# ==================== fig9：更新方式取舍（含遗忘代价） ====================
fig = plt.figure(figsize=(7.2, 4.3), dpi=150)
fig.patch.set_facecolor(BG)
fig.suptitle("同一批湿实验数据，三种更新方式的代价账（实测）", fontsize=14.5, fontweight="bold", color=INK, y=0.97)
gs9 = gridspec.GridSpec(1, 2, width_ratios=[1.25, 1.0], wspace=0.32,
                        left=0.10, right=0.96, top=0.80, bottom=0.16)

ax = fig.add_subplot(gs9[0])
labels = ["不更新", "全参\n微调", "LoRA\nlr=1e-3", "LoRA\nlr=3e-4"]
new_gain = [0.0, 0.195, 0.168, 0.098]
old_loss = [0.0, 0.023, 0.090, -0.011]
x = np.arange(len(labels))
bars = ax.bar(x - 0.19, new_gain, width=0.36, color=GREEN, label="新系列 RMSE 改善")
bars2 = ax.bar(x + 0.19, old_loss, width=0.36, color=RED, label="历史系列变化（正=退化）")
for b, v in zip(bars, new_gain):
    if abs(v) > 1e-9:
        ax.text(b.get_x() + b.get_width() / 2, v + 0.006, "%+.3f" % v, ha="center", fontsize=8.5, color=GREEN)
for b, v in zip(bars2, old_loss):
    if abs(v) > 1e-9:
        ax.text(b.get_x() + b.get_width() / 2, v + (0.006 if v >= 0 else -0.024), "%+.3f" % v,
                ha="center", fontsize=8.5, color=RED)
ax.text(x[0], 0.006, "基线", ha="center", fontsize=8.5, color=GREY)
ax.axhline(0, color="#C7D2D6", lw=1.0)
ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=9.5)
ax.set_title("累计回传 250 个分子时", fontsize=11, color="#3F3F3F", pad=6)
ax.legend(frameon=False, fontsize=8.5, loc="upper right")
ax.set_ylim(-0.06, 0.27)
ax.tick_params(labelsize=9, colors=GREY)
for s in ["top", "right"]:
    ax.spines[s].set_visible(False)
for s in ["left", "bottom"]:
    ax.spines[s].set_color("#E0E3E6")

ax2 = fig.add_subplot(gs9[1])
ax2.set_xlim(0, 1); ax2.set_ylim(0, 1); ax2.axis("off")
ax2.set_title("一个底座 + 多个系列适配器", fontsize=11.5, fontweight="bold", color=INK, pad=8)
rows = [("底座权重", "1342.1 KB", "全项目共用，冻结不动"),
        ("单个 LoRA 适配器", "57.3 KB", "每个系列一份，随时插拔"),
        ("可训练参数占比", "4.18%", "全微调需动 100% 参数"),
        ("多个系列互不干扰", "实测：是", "底座权重前后完全一致")]
y = 0.86
for name, val, note in rows:
    ax2.add_patch(FancyBboxPatch((0.02, y - 0.135), 0.96, 0.24, boxstyle="round,pad=0.01,rounding_size=0.06",
                                 linewidth=1.1, edgecolor="#DDE3E6", facecolor=CARD))
    ax2.text(0.06, y + 0.035, name, fontsize=10.5, color=INK, fontweight="bold")
    ax2.text(0.94, y + 0.035, val, fontsize=11, color=GREEN, fontweight="bold", ha="right")
    ax2.text(0.06, y - 0.072, note, fontsize=9, color=GREY)
    y -= 0.195
ax2.text(0.02, 0.035, "低秩不等于天然防遗忘：lr 调到 1e-3 时，适配器对历史\n知识的破坏 (+0.090) 反而比全参微调 (+0.023) 更大。",
         fontsize=9.5, color=RED, linespacing=1.6)

fig.savefig(OUT + "fig9_lora_tradeoff.png", facecolor=BG, bbox_inches="tight", pad_inches=0.18)
plt.close(fig)
print("done")
