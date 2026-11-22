# -*- coding: utf-8 -*-
# 文章关联：AI制药第13篇_分子生成_公众号文章.html
# 脚本职责：封面：公众号（1920x817，2.35:1）
# 使用注解：沿用行业篇/细节篇②/第11篇 封面模板（背景渐变 / 分子网络 / 胶囊样式），
#           主轴 = 第13篇反转结论：追分臂没有把候选推出适用域
# 数字来源：e27_audit_extras.json（逐对 +0.745 / +0.071；MW +83.4；SA +0.909）、
#           e25_closedloop.json（4342 查询）、e26_disagreement.json（914 候选）
# 完整对应关系与运行顺序：README.md
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch

plt.rcParams["font.family"] = ["Hiragino Sans GB", "Arial Unicode MS"]
plt.rcParams["font.sans-serif"] = ["Hiragino Sans GB", "Arial Unicode MS"]
plt.rcParams["axes.unicode_minus"] = False

OUT = "/Users/dujiepeng/WorkBuddy/Work/AI制药第13篇_公众号封面_1920x817.jpg"

fig = plt.figure(figsize=(19.2, 8.17), dpi=100)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")

# ---------- 背景渐变（与行业篇同配方） ----------
c1 = np.array([9, 46, 40]) / 255.0
c2 = np.array([15, 168, 126]) / 255.0
nx = ny = 256
gx = np.linspace(0, 1, nx).reshape(1, nx)
gy = np.linspace(0, 1, ny).reshape(ny, 1)
t = np.clip(0.75 * gx + 0.55 * gy, 0, 1)
img = c1.reshape(1, 1, 3) * (1 - t)[:, :, None] + c2.reshape(1, 1, 3) * t[:, :, None]
ax.imshow(img, extent=[0, 1, 0, 1], aspect="auto", origin="lower", zorder=1)

# ---------- 右上抽象分子网络（限制在右侧，避开左侧标题区） ----------
rng = np.random.default_rng(13)
pts = [(0.745 + rng.uniform(0, 0.20), 0.42 + rng.uniform(0, 0.50)) for _ in range(17)]
for i, (x1, y1) in enumerate(pts):
    for j in range(i + 1, len(pts)):
        x2, y2 = pts[j]
        d = ((x1 - x2) ** 2 + (y1 - y2) ** 2) ** 0.5
        if d < 0.11:
            ax.plot([x1, x2], [y1, y2], color="white", alpha=0.18, lw=1.4, zorder=2)
for k, (x, y) in enumerate(pts):
    face = "#F5A623" if k % 6 == 0 else "white"
    ax.add_patch(Circle((x, y), 0.0125, facecolor=face, alpha=0.55 if face == "white" else 0.95,
                        edgecolor="none", zorder=3))

# ---------- 编号标 ----------
ax.add_patch(FancyBboxPatch((0.055, 0.878), 0.012, 0.055,
                            boxstyle="round,pad=0.002,rounding_size=0.004",
                            facecolor="#F5A623", edgecolor="none", zorder=5))
ax.text(0.085, 0.905, "零基础入门 AI 制药 · ⑬ 分子生成", fontsize=25,
        color="#CFF3E7", va="center", zorder=5)

# ---------- 两行大标题（白 + 黄） ----------
ax.text(0.055, 0.775, "高分是发现，", fontsize=68, color="white",
        fontweight="bold", va="center", zorder=5)
ax.text(0.055, 0.625, "还是骗过了裁判？", fontsize=68, color="#F5C542",
        fontweight="bold", va="center", zorder=5)

# ---------- 副标题（文章思路一句话） ----------
ax.text(0.058, 0.495, "把评分器装进生成循环：同一起点、同一动作，只改保留规则", fontsize=24,
        color="#DFF7EF", va="center", zorder=5)

# ---------- 右侧注释（反转结论，对齐细节篇②的引线注释） ----------
ax.plot([0.700, 0.745], [0.330, 0.330], color="#F5C542", lw=2, zorder=5)
ax.text(0.755, 0.330, "逐对实测：追分臂相似度不降反升", fontsize=19,
        color="#F5C542", va="center", zorder=5)
ax.text(0.985, 0.115, "4342 次评分查询 · 914 个候选 · 全部本机实跑",
        fontsize=15, color="#9FDCC8", ha="right", va="center", zorder=5)

# ---------- 关键胶囊（核心摘要三件事：反转 / 代价 / 分歧） ----------
items = ["追分臂 分数 +0.745", "相似度反升 +0.071", "结构代价 MW +83"]
x0 = 0.058
for it in items:
    w = 0.0125 * len(it) + 0.045
    ax.add_patch(FancyBboxPatch((x0, 0.135), w, 0.105,
                                boxstyle="round,pad=0.004,rounding_size=0.02",
                                facecolor="white", alpha=0.16, edgecolor="none",
                                zorder=5))
    ax.text(x0 + w / 2, 0.1875, it, fontsize=20.5, color="white",
            ha="center", va="center", zorder=6)
    x0 += w + 0.022

fig.savefig(OUT, format="jpg", dpi=100, facecolor="#0B3B32")
print("saved", OUT)
