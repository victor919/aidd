# -*- coding: utf-8 -*-
# 文章关联：AI制药实践篇B_实操篇_公众号文章.html
# 脚本职责：封面：公众号初版
# 使用注解：与 v2 写入同一 JPG；历史初版应先运行，避免覆盖 v2。
# 完整对应关系与运行顺序：README.md
"""AI制药实践篇B（实操篇）公众号封面（1920x817，2.35:1）——沿用行业篇系列模板"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch

plt.rcParams["font.sans-serif"] = ["Hiragino Sans GB"]
plt.rcParams["axes.unicode_minus"] = False

OUT = "/Users/dujiepeng/WorkBuddy/Work/AI制药实践篇B_公众号封面_1920x817.jpg"

fig = plt.figure(figsize=(19.2, 8.17), dpi=100)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")

# 背景渐变：与行业篇同一配方（左上深墨绿 → 右下亮绿）
c1 = np.array([9, 46, 40]) / 255.0     # #092E28
c2 = np.array([15, 168, 126]) / 255.0  # #0FA87E
nx, ny = 256, 256
gx = np.linspace(0, 1, nx).reshape(1, nx)
gy = np.linspace(0, 1, ny).reshape(ny, 1)
t = np.clip(0.75 * gx + 0.55 * gy, 0, 1)
img = c1.reshape(1, 1, 3) * (1 - t)[:, :, None] + c2.reshape(1, 1, 3) * t[:, :, None]
ax.imshow(img, extent=[0, 1, 0, 1], aspect="auto", origin="lower", zorder=1)

# 右上抽象分子网络（同一配方，限制在右上，避开数据条）
rng = np.random.default_rng(7)
pts = [(0.755 + rng.uniform(0, 0.185), 0.55 + rng.uniform(0, 0.34)) for _ in range(15)]
for i, (x1, y1) in enumerate(pts):
    for j in range(i + 1, len(pts)):
        x2, y2 = pts[j]
        d = ((x1 - x2) ** 2 + (y1 - y2) ** 2) ** 0.5
        if d < 0.10:
            ax.plot([x1, x2], [y1, y2], color="white", alpha=0.18, lw=1.4, zorder=2)
for x, y in pts:
    ax.add_patch(Circle((x, y), 0.0125, facecolor="white", alpha=0.55,
                        edgecolor="none", zorder=3))
ax.add_patch(Circle((0.83, 0.68), 0.016, facecolor="#F5A623", alpha=0.95,
                    edgecolor="none", zorder=4))

# 编号标（橙色竖条 + 系列名）
ax.add_patch(FancyBboxPatch((0.055, 0.735), 0.012, 0.055,
                            boxstyle="round,pad=0.002,rounding_size=0.004",
                            facecolor="#F5A623", edgecolor="none", zorder=5))
ax.text(0.085, 0.762, "零基础入门 AI 制药 · ⑥下 · 实操篇", fontsize=25, color="#CFF3E7",
        va="center", zorder=5)

# 大标题
ax.text(0.055, 0.585, "亲手跑通 AIDD", fontsize=104, color="white",
        fontweight="bold", va="center", zorder=5)

# 副标题
ax.text(0.058, 0.405, "六个坑 · 三处硬伤 · 一条真跑通的流水线",
        fontsize=34, color="#EAFBF5", va="center", zorder=5)

# 数据胶囊（全部来自本机实跑数字）
items = ["写死的坐标，偏了 31.4 埃", "对接第 1 名，真实活性第 18",
         "10,468 个分子，全程可复现"]
x0 = 0.058
for i, it in enumerate(items):
    w = 0.0125 * len(it) + 0.045
    ax.add_patch(FancyBboxPatch((x0, 0.185), w, 0.115,
                                boxstyle="round,pad=0.004,rounding_size=0.02",
                                facecolor="white", alpha=0.16, edgecolor="none",
                                zorder=5))
    ax.text(x0 + w / 2, 0.2425, it, fontsize=20.5, color="white",
            ha="center", va="center", zorder=6)
    x0 += w + 0.022

fig.savefig(OUT, format="jpg", dpi=100, facecolor="#0B3B32")
print("saved", OUT)
