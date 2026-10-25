# -*- coding: utf-8 -*-
# 文章关联：AI制药实践篇B_实操篇_公众号文章.html
# 脚本职责：封面：小红书竖版
# 使用注解：本篇跨平台封面，不是实验或正文配图。
# 完整对应关系与运行顺序：README.md
"""AI制药实践篇B（实操篇）小红书封面（1080x1440，3:4 竖版）
沿用系列视觉：深绿渐变 + 白色大字 + 胶囊
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch

plt.rcParams["font.sans-serif"] = ["Hiragino Sans GB"]
plt.rcParams["axes.unicode_minus"] = False

OUT = "/Users/dujiepeng/WorkBuddy/Work/AI制药实践篇B_小红书封面_1080x1440.jpg"

fig = plt.figure(figsize=(10.8, 14.4), dpi=100)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")

# 背景渐变（与系列同配方）
c1 = np.array([9, 46, 40]) / 255.0
c2 = np.array([15, 168, 126]) / 255.0
nx = ny = 256
gx = np.linspace(0, 1, nx).reshape(1, nx)
gy = np.linspace(0, 1, ny).reshape(ny, 1)
t = np.clip(0.75 * gx + 0.55 * gy, 0, 1)
img = c1.reshape(1, 1, 3) * (1 - t)[:, :, None] + c2.reshape(1, 1, 3) * t[:, :, None]
ax.imshow(img, extent=[0, 1, 0, 1], aspect="auto", origin="lower", zorder=1)

# 右下抽象分子网络（右下为空白区，避开标题与胶囊）
rng = np.random.default_rng(7)
pts = [(0.72 + rng.uniform(0, 0.23), 0.05 + rng.uniform(0, 0.16)) for _ in range(16)]
for i, (x1, y1) in enumerate(pts):
    for j in range(i + 1, len(pts)):
        x2, y2 = pts[j]
        d = ((x1 - x2) ** 2 + (y1 - y2) ** 2) ** 0.5
        if d < 0.135:
            ax.plot([x1, x2], [y1, y2], color="white", alpha=0.18, lw=1.6, zorder=2)
for x, y in pts:
    ax.add_patch(Circle((x, y), 0.017, facecolor="white", alpha=0.55,
                        edgecolor="none", zorder=3))
ax.add_patch(Circle((0.87, 0.175), 0.022, facecolor="#F5A623", alpha=0.95,
                    edgecolor="none", zorder=4))

# 编号标
ax.add_patch(FancyBboxPatch((0.075, 0.925), 0.014, 0.040,
                            boxstyle="round,pad=0.002,rounding_size=0.005",
                            facecolor="#F5A623", edgecolor="none", zorder=5))
ax.text(0.105, 0.945, "零基础入门 AI 制药 · ⑥下", fontsize=30,
        color="#CFF3E7", va="center", zorder=5)

# 三行大标题
for y, s in [(0.845, "一台笔记本，"), (0.745, "跑通 AI 制药"), (0.645, "流水线")]:
    ax.text(0.075, y, s, fontsize=96, color="white", fontweight="bold",
            va="center", zorder=5)

# 摘要五行
lines = [
    "本文以 EGFR 为例，从靶点与口袋、数据清洗、",
    "骨架切分起步，走到模型对照、分子对接",
    "与主动学习模拟，逐个拆解六个实操陷阱，",
    "并附上工程结构与关键代码——把零散",
    "工具串成可运行、可测试的研发实践地图。",
]
for i, ln in enumerate(lines):
    ax.text(0.075, 0.50 - i * 0.058, ln, fontsize=28, color="#DFF7EF",
            va="center", zorder=5)

# 竖排关键词胶囊
items = ["六个实操陷阱", "10468 个分子实跑", "全程 seed=42 可复现"]
for i, it in enumerate(items):
    y = 0.185 - i * 0.070
    ax.add_patch(FancyBboxPatch((0.075, y), 0.62, 0.058,
                                boxstyle="round,pad=0.004,rounding_size=0.028",
                                facecolor="white", alpha=0.16, edgecolor="none",
                                zorder=5))
    ax.text(0.105, y + 0.029, it, fontsize=30, color="white",
            va="center", zorder=6)

fig.savefig(OUT, format="jpg", dpi=100, facecolor="#0B3B32")
print("saved", OUT)
