# -*- coding: utf-8 -*-
# 文章关联：AI制药第13篇_分子生成_公众号文章.html
# 脚本职责：封面 v2：公众号（1920x817，2.35:1），微信绿风格（深墨绿+金色+衬线+分子轨迹）
# 版式来源：用户提供参考图（深墨绿底/宋体标题/金色轨迹线/球棍分子/叶片装饰）
#           文案不用参考图原句——其"最高分终点回到训练集"与实测相反；
#           本封面文案全部来自 e25/e26/e27 实测（27 组、终点回查 0 命中）
# 数字来源：e27_audit_extras.json（endpoint_identity 0 命中）、e25_closedloop.json（27 组运行）
# 渲染注意：非等比轴下 patches.Circle 会变形 → 原子/轨迹点一律 ax.plot(marker="o")（A3 坑）
#           原子位置用像素坐标定义，经 px() 转数据坐标；markersize 单位 pt = px*0.72
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, FancyBboxPatch, FancyArrowPatch
from matplotlib.path import Path
import matplotlib.patches as mpatches

plt.rcParams["font.family"] = ["Songti SC", "STKaiti", "Hiragino Sans GB", "Arial Unicode MS"]
plt.rcParams["font.sans-serif"] = ["Hiragino Sans GB", "Arial Unicode MS"]
plt.rcParams["axes.unicode_minus"] = False

W, H = 1920, 817
OUT = "/Users/dujiepeng/WorkBuddy/Work/AI制药第13篇_公众号封面_1920x817.jpg"

GOLD = "#D9A441"
GOLD_L = "#E9C46A"
CREAM = "#F3EFE3"
INK_N = "#4A90D9"   # 氮
INK_O = "#D9503F"   # 氧
CARBON = "#E8E4DA"  # 碳

def px(xp, yp):
    """像素坐标(左上原点 y 向下) → 数据坐标(左下原点)"""
    return xp / W, 1.0 - yp / H

fig = plt.figure(figsize=(19.2, 8.17), dpi=100)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")

# ---------- 背景渐变（深墨绿） + 暗角（直接混入像素，避免灰雾蒙层） ----------
c1 = np.array([7, 30, 25]) / 255.0
c2 = np.array([22, 56, 44]) / 255.0
nx = ny = 256
gx = np.linspace(0, 1, nx).reshape(1, nx)
gy = np.linspace(0, 1, ny).reshape(ny, 1)
t = np.clip(0.7 * gx + 0.5 * gy, 0, 1)
img = c1.reshape(1, 1, 3) * (1 - t)[:, :, None] + c2.reshape(1, 1, 3) * t[:, :, None]
cx_, cy_ = nx / 2, ny / 2
jg, ig = np.meshgrid(np.arange(nx), np.arange(ny))
d = np.sqrt(((jg - cx_) / cx_) ** 2 + ((ig - cy_) / cy_) ** 2)
vig = np.clip((d - 0.55) / 0.75, 0, 1) * 0.45
img = img * (1 - vig)[..., None]          # 暗角：边缘乘法压暗，不发灰
ax.imshow(img, extent=[0, 1, 0, 1], aspect="auto", origin="lower", zorder=1)

# ---------- 叶片装饰（角落，线描风：轮廓+中脉+侧脉细线，不画叶肉色块） ----------
def leaf(xp, yp, wp, hp, angle_deg, alpha=0.8):
    xd, yd = px(xp, yp)
    a = np.deg2rad(angle_deg)
    ux, uy = np.cos(a) * wp / 2 / W, np.sin(a) * wp / 2 / H
    vx, vy = -np.sin(a) * hp / 2 / W, np.cos(a) * hp / 2 / H
    # 轮廓（两条对称二次弧近似叶缘）
    n = 40
    tt = np.linspace(-1, 1, n)
    for sgn in (1, -1):
        ex = xd + tt * ux + sgn * (1 - tt ** 2) * vx * 1.9
        ey = yd + tt * uy + sgn * (1 - tt ** 2) * vy * 1.9
        ax.plot(ex, ey, color=GOLD, alpha=alpha * 0.55, lw=1.1, zorder=3)
    # 中脉 + 侧脉
    ax.plot([xd - ux, xd + ux], [yd - vy, yd + vy], color=GOLD,
            alpha=alpha * 0.5, lw=1.0, zorder=3)
    for f in (-0.5, -0.2, 0.1, 0.4):
        mx, my = xd + f * ux, yd + f * vy
        for sgn in (1, -1):
            ax.plot([mx, mx + (vx * 1.35 + ux * 0.22) * sgn],
                    [my, my + (vy * 1.35 + uy * 0.22) * sgn],
                    color=GOLD, alpha=alpha * 0.35, lw=0.9, zorder=3)

leaf(-30, 780, 320, 110, 20)
leaf(140, 850, 250, 90, -6, alpha=0.6)
leaf(1950, 790, 340, 115, 162)
leaf(1840, 850, 240, 85, 148, alpha=0.6)
leaf(1985, 130, 290, 100, 198, alpha=0.5)

# ---------- 球棍分子 ----------
def molecule(cx_p, cy_p, ring_r, atoms_spec, bonds, sub_labels):
    """cx_p/cy_p 中心像素坐标; atoms_spec: [(dx,dy,color,size_pt,depth)]; bonds: [(i,j)]"""
    pts = []
    for dx, dy, color, size, z in atoms_spec:
        xd, yd = px(cx_p + dx, cy_p - dy)
        pts.append((xd, yd))
    for i, j in bonds:
        (x1, y1), (x2, y2) = pts[i], pts[j]
        ax.plot([x1, x2], [y1, y2], color="#BFC7C2", alpha=0.75, lw=2.6,
                solid_capstyle="round", zorder=6)
    for (dx, dy, color, size, z), (xd, yd) in zip(atoms_spec, pts):
        ax.plot(xd, yd, "o", ms=size, mfc=color, mec="#0B2B23", mew=1.0, zorder=z)

def ring_atoms(cx_p, cy_p, r, colors, sizes):
    spec = []
    for k in range(6):
        a = np.deg2rad(90 - k * 60)
        spec.append((r * np.cos(a), r * np.sin(a), colors[k], sizes[k], 7))
    return spec

# 分子 1（左侧小，六元环 + 1 个红 O 取代基）
col1 = [CARBON, INK_N, CARBON, CARBON, CARBON, INK_N]
sz1 = [20, 18, 20, 20, 20, 18]
spec1 = ring_atoms(1265, 400, 50, col1, sz1)
spec1.append((0, 50 + 32, INK_O, 17, 7))          # 顶部红 O
spec1.append((50 * np.cos(np.deg2rad(30)), 50 * np.sin(np.deg2rad(30)) + 36, CARBON, 14, 6))  # 右上甲基
b1 = [(i, (i + 1) % 6) for i in range(6)] + [(0, 6), (2, 7)]
molecule(1265, 400, 50, spec1, b1, None)

# 分子 2（右侧大，六元环 + 红 O 与碳取代基）
col2 = [CARBON, CARBON, INK_N, CARBON, CARBON, INK_N]
sz2 = [25, 25, 22, 25, 25, 22]
spec2 = ring_atoms(1782, 535, 64, col2, sz2)
spec2.append((64 * np.cos(np.deg2rad(60)), 64 * np.sin(np.deg2rad(60)) + 40, INK_O, 20, 8))
spec2.append((-64 * np.cos(np.deg2rad(60)) - 36, -64 * np.sin(np.deg2rad(60)) * 0.6, CARBON, 16, 6))
b2 = [(i, (i + 1) % 6) for i in range(6)] + [(1, 6), (4, 7)]
molecule(1782, 535, 64, spec2, b2, None)

# ---------- 金色平台光环 ----------
for (cxp, cyp, w, h) in [(1265, 496, 170, 44), (1782, 650, 210, 54)]:
    xd, yd = px(cxp, cyp)
    ax.add_patch(Ellipse((xd, yd), w / W, h / H, facecolor="none",
                         edgecolor=GOLD, alpha=0.9, lw=2.4, zorder=5))
    ax.add_patch(Ellipse((xd, yd), w / W, h / H, facecolor=GOLD, alpha=0.16,
                         edgecolor="none", zorder=5))

# ---------- 金色搜索轨迹（贝塞尔平滑曲线 + 圆点 + 箭头） ----------
verts = [px(1318, 468), px(1450, 350), px(1620, 330), px(1700, 420),
         px(1745, 505), px(1640, 545), px(1610, 460), px(1700, 400)]
codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
         Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CURVE4]
path = Path(verts, codes)
# 手工采样 bezier（三段四次? CURVE4 为三次贝塞尔：MOVETO + 7 个 CURVE4 = 7 段? 实为 7/3 不整除）
# 改为手动构造：两段三次贝塞尔
def bez(p0, p1, p2, p3, n=40):
    t = np.linspace(0, 1, n)[:, None]
    return ((1 - t) ** 3 * p0 + 3 * (1 - t) ** 2 * t * p1 +
            3 * (1 - t) * t ** 2 * p2 + t ** 3 * p3)

seg1 = bez(np.array(verts[0]), np.array(verts[1]), np.array(verts[2]), np.array(verts[3]), 60)
seg2 = bez(np.array(verts[3]), np.array(verts[4]), np.array(verts[5]), np.array(verts[6]), 60)
seg3 = bez(np.array(verts[6]), np.array(verts[7]), np.array(px(1740, 430)), np.array(px(1738, 470)), 40)
curve = np.vstack([seg1, seg2, seg3])
ax.plot(curve[:, 0], curve[:, 1], color=GOLD_L, alpha=0.9, lw=2.0, ls=(0, (4, 3)), zorder=6)
for k in range(0, len(curve), 5):
    ax.plot(curve[k, 0], curve[k, 1], "o", ms=6, mfc=GOLD_L, mec="none", alpha=0.95, zorder=7)

arrow = FancyArrowPatch(px(1700, 400), px(1738, 458), arrowstyle="-|>",
                        mutation_scale=30, color=GOLD, lw=2.4, zorder=8)
ax.add_patch(arrow)

# ---------- 编号标（金竖条 + 衬线小标） ----------
xb, yb = px(96, 62)
ax.add_patch(FancyBboxPatch((xb, yb - 0.028), 0.005, 0.056,
                            boxstyle="round,pad=0.001,rounding_size=0.002",
                            facecolor=GOLD, edgecolor="none", zorder=5))
ax.text(0.0705, yb, "AI 制药 · 第 13 篇", fontsize=25, color=GOLD,
        va="center", zorder=5, family=["Songti SC", "STKaiti"])

# ---------- 两行大标题（宋体，米白） ----------
ax.text(0.050, 0.745, "追着高分找分子，", fontsize=74, color=CREAM,
        fontweight="bold", va="center", zorder=5, family=["Songti SC", "STKaiti"])
ax.text(0.050, 0.565, "为什么没越跑越偏？", fontsize=74, color=CREAM,
        fontweight="bold", va="center", zorder=5, family=["Songti SC", "STKaiti"])

# ---------- 副标题 + 金线 ----------
ax.text(0.0525, 0.405, "同一起点，随机选 vs 按分数选", fontsize=27, color="#D8D2C0",
        va="center", zorder=5, family=["Songti SC", "STKaiti"])
ax.plot([0.0525, 0.415], [0.345, 0.345], color=GOLD, lw=1.6, alpha=0.9, zorder=5)

# ---------- 底部金色圆角框横条（实测结论） ----------
bx, by = px(96, 700)
bw, bh = 0.472, 0.115
ax.add_patch(FancyBboxPatch((bx, by), bw, bh,
                            boxstyle="round,pad=0.004,rounding_size=0.018",
                            facecolor="#0E2C24", alpha=0.55, edgecolor=GOLD,
                            lw=2.0, zorder=5))
ax.text(bx + bw / 2, by + bh / 2, "27 组对照搜索，追分终点无一回到训练集",
        fontsize=31, color=GOLD_L, ha="center", va="center", zorder=6,
        family=["Songti SC", "STKaiti"])

fig.savefig(OUT, format="jpg", dpi=100, facecolor="#0B2B23")
print("saved", OUT)
