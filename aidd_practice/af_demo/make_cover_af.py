# -*- coding: utf-8 -*-
# 文章关联：AI制药细节篇1_靶点与口袋_公众号文章.html
# 脚本职责：封面：公众号横版与小红书竖版
# 使用注解：使用真实 AF Cα 坐标投影；依赖本地结构，不是正文实验。
# 完整对应关系与运行顺序：../README.md
"""《细节篇①：给蛋白质「开天眼」》公众号封面（1920x817，2.35:1）+ 小红书封面（1080x1440，3:4）

沿用系列模板：深墨绿 → 亮绿渐变背景、金色编号标、数据胶囊条。
本篇专属视觉：右侧蛋白点云 = **AlphaFold P00533 v6 的真实 Cα 坐标**（SVD 投影到二维），
按 pLDDT 分档给不同透明度——高置信区亮、低置信区如雾；橙色光圈 + 黄色小分子标记口袋位置。
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

plt.rcParams["font.sans-serif"] = ["Hiragino Sans GB"]
plt.rcParams["axes.unicode_minus"] = False

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.abspath(os.path.join(HERE, "..", ".."))
PDB_AF = os.path.join(HERE, "struct", "egfr_af_v6.pdb")

DEEP = np.array([9, 46, 40]) / 255.0      # #092E28
MAIN = np.array([15, 168, 126]) / 255.0   # #0FA87E
AMBER = "#F5A623"


# ---------------------------------------------------------------- 真实蛋白点云
def load_ca(path):
    xyz, b = [], []
    for l in open(path):
        if l.startswith("ATOM") and l[12:16].strip() == "CA":
            xyz.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
            b.append(float(l[60:66]))
    return np.array(xyz), np.array(b)


def project2d(xyz):
    """SVD 取前两个主成分方向做正交投影；归一化用统一缩放因子，保持原始纵横比"""
    x = xyz - xyz.mean(0)
    _, _, vt = np.linalg.svd(x, full_matrices=False)
    p = x @ vt[:2].T
    p = p - p.min(0)
    p = p / p.max()          # 关键：同一标量缩放，否则点云会被拉变形
    return p


def fit_box(p, box, axes_w, axes_h):
    """把保持比例的点云等比缩放并居中放进归一化框 box=(x0,y0,x1,y1)"""
    x0, y0, x1, y1 = box
    bw_in = (x1 - x0) * axes_w      # 框的物理尺寸（inch）
    bh_in = (y1 - y0) * axes_h
    du = float(p[:, 0].max() - p[:, 0].min())
    dv = float(p[:, 1].max() - p[:, 1].min())
    ar = du / dv                     # 点云自身纵横比
    if bw_in / bh_in > ar:           # 框更"胖" → 以高为限
        h_in, w_in = bh_in, bh_in * ar
    else:                            # 框更"瘦" → 以宽为限
        w_in, h_in = bw_in, bw_in / ar
    cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    w_n, h_n = w_in / axes_w, h_in / axes_h
    u = (p[:, 0] - p[:, 0].min()) / du * w_n + (cx - w_n / 2)
    v = (p[:, 1] - p[:, 1].min()) / dv * h_n + (cy - h_n / 2)
    return u, v


def text_w_n(s, fs, axes_w):
    """估算一段文字在归一化坐标下的宽度：全角 1.0 / ASCII 字母数字 0.55 / 空格 0.3"""
    n = sum(1.0 if ord(c) > 0x2E80 else (0.30 if c == " " else 0.55) for c in s)
    return n * fs / 72.0 / axes_w


def draw_protein(ax, p, b, box, axes_w, axes_h, scale=1.0, halo=None,
                 halo_ms=(112, 68)):
    """按 pLDDT 分档绘制白色点云；halo 为 (x,y) 时在该处加橙色光圈 + 黄色小分子"""
    u, v = fit_box(p, box, axes_w, axes_h)
    bands = [(90, 101, 0.88, 11), (70, 90, 0.62, 8.5),
             (50, 70, 0.44, 7), (0, 50, 0.22, 6)]
    for lo, hi, alpha, ms in bands:
        m = (b >= lo) & (b < hi)
        ax.plot(u[m], v[m], linestyle="none", marker="o", ms=ms * scale,
                color="white", alpha=alpha, zorder=3)
    if halo is not None:
        hx, hy = halo
        ms_out, ms_in = halo_ms
        ax.plot([hx], [hy], marker="o", ms=ms_out * scale, color=AMBER,
                alpha=0.20, zorder=4)
        ax.plot([hx], [hy], marker="o", ms=ms_in * scale, color=AMBER,
                alpha=0.42, zorder=5)
        # 口袋里的小分子（示意：三原子折线）
        d = 0.021
        xs = [hx - d, hx, hx + d * 0.8]
        ys = [hy - d * 0.75, hy + d * 0.15, hy + d * 0.95]
        ax.plot(xs, ys, "-", color="#FFD489", lw=3.4 * scale,
                solid_capstyle="round", zorder=6)
        ax.plot(xs, ys, marker="o", ms=13 * scale, color="#FFD489",
                linestyle="none", zorder=7)
    return u, v


def gradient(ax, c1=DEEP, c2=MAIN, k=(0.75, 0.55)):
    nx = ny = 256
    gx = np.linspace(0, 1, nx).reshape(1, nx)
    gy = np.linspace(0, 1, ny).reshape(ny, 1)
    t = np.clip(k[0] * gx + k[1] * gy, 0, 1)
    img = c1.reshape(1, 1, 3) * (1 - t)[:, :, None] + c2.reshape(1, 1, 3) * t[:, :, None]
    ax.imshow(img, extent=[0, 1, 0, 1], aspect="auto", origin="lower", zorder=1)


def tag(ax, x, y, text, fs, color="#CFF3E7", bar=AMBER, z=5):
    ax.add_patch(FancyBboxPatch((x, y - 0.026), 0.011, 0.052,
                                boxstyle="round,pad=0.002,rounding_size=0.004",
                                facecolor=bar, edgecolor="none", zorder=z))
    ax.text(x + 0.030, y, text, fontsize=fs, color=color, va="center", zorder=z)


def pill(ax, fig, text, x, y, fs, axes_w, pad=0.034, h=0.105, z=5):
    """画一个半透明数据胶囊。宽度 = 渲染实测文字宽 + 两侧 pad，绝不溢出也不空旷。
    返回胶囊右端 x。"""
    dpi = fig.dpi
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    probe = ax.text(0.5, 0.5, text, fontsize=fs, color="white",
                    ha="center", va="center")
    w = probe.get_window_extent(renderer=r).width / dpi / axes_w + 2 * pad
    probe.remove()
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle="round,pad=0.004,rounding_size=0.02",
                                facecolor="white", alpha=0.16, edgecolor="none",
                                zorder=z))
    ax.text(x + w / 2, y + h / 2, text, fontsize=fs, color="white",
            ha="center", va="center", zorder=z + 1)
    return x + w


def pills_row(ax, fig, items, y, x0, fs, axes_w, gap=0.024, **kw):
    """横排胶囊条"""
    for it in items:
        x0 = pill(ax, fig, it, x0, y, fs, axes_w, **kw) + gap
    return x0 - gap


# ---------------------------------------------------------------- 公众号封面 2.35:1
def cover_wechat(p, b):
    W, H = 19.2, 8.17
    out = os.path.join(WORK, "AI制药细节篇1_公众号封面_1920x817.jpg")
    fig = plt.figure(figsize=(W, H), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    gradient(ax)

    # 右侧：真实 AlphaFold 点云 + 口袋标记（口袋位置取点云右侧 28% 区域的中位高度，示意）
    # 框下沿抬到 0.32，给左下角胶囊条留出净空，避免胶囊压在点云上
    box = (0.585, 0.32, 0.98, 0.95)
    u, v = fit_box(p, box, W, H)
    hx = float(np.quantile(u, 0.86))
    hy = float(np.median(v[u > np.quantile(u, 0.72)]))
    draw_protein(ax, p, b, box, W, H, scale=1.0, halo=(hx, hy),
                 halo_ms=(100, 60))

    # 左侧文字
    tag(ax, 0.055, 0.880, "零基础入门 AI 制药 · ⑨ 细节篇①", 25)
    ax.text(0.055, 0.660, "给蛋白质「开天眼」", fontsize=82, color="white",
            fontweight="bold", va="center", zorder=6)
    ax.text(0.058, 0.470, "从 AlphaFold 结构预测，到活性口袋识别", fontsize=34,
            color="#EAFBF5", va="center", zorder=6)
    ax.text(0.058, 0.372, "全长模型与实验结构差在哪？pLDDT 能指出口袋吗？",
            fontsize=22.5, color="#DFF7EF", va="center", zorder=6)

    end = pills_row(ax, fig, ["pLDDT 前 200 名命中 0 个", "口袋是配体邻域的 6.4 倍",
                              "口袋内部 RMSE 2.858 埃"], 0.155, 0.058, 20.5, W)
    print("    胶囊条右端 x = %.3f（应 < 0.95）" % end)

    fig.savefig(out, format="jpg", dpi=100, facecolor="#0B3B32")
    plt.close(fig)
    print("saved", out, os.path.getsize(out) // 1024, "KB")


# ---------------------------------------------------------------- 小红书封面 3:4
def cover_xhs(p, b):
    W, H = 10.8, 14.4
    out = os.path.join(WORK, "AI制药细节篇1_小红书封面_1080x1440.jpg")
    fig = plt.figure(figsize=(W, H), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    gradient(ax, k=(0.35, 0.85))

    # 上部文字，下部点云
    tag(ax, 0.075, 0.900, "零基础入门 AI 制药 · ⑨ 细节篇①", 20)
    ax.text(0.075, 0.755, "给蛋白质", fontsize=72, color="white",
            fontweight="bold", va="center", zorder=6)
    ax.text(0.075, 0.635, "「开天眼」", fontsize=72, color="white",
            fontweight="bold", va="center", zorder=6)
    ax.text(0.078, 0.545, "从 AlphaFold 结构预测", fontsize=30,
            color="#EAFBF5", va="center", zorder=6)
    ax.text(0.078, 0.487, "到活性口袋识别", fontsize=30,
            color="#EAFBF5", va="center", zorder=6)

    box = (0.10, 0.195, 0.92, 0.465)
    u, v = fit_box(p, box, W, H)
    hx = float(np.quantile(u, 0.86))
    hy = float(np.median(v[u > np.quantile(u, 0.72)]))
    draw_protein(ax, p, b, box, W, H, scale=1.35, halo=(hx, hy),
                 halo_ms=(96, 58))

    items = ["pLDDT 前 200 名命中 0 个", "口袋是配体邻域的 6.4 倍"]
    y = 0.132
    for it in items:
        end = pill(ax, fig, it, 0.075, y, 17.0, W, pad=0.036, h=0.058)
        print("    小红书胶囊右端 %.3f（应 < 0.93）  %s" % (end, it))
        y -= 0.080

    fig.savefig(out, format="jpg", dpi=100, facecolor="#0B3B32")
    plt.close(fig)
    print("saved", out, os.path.getsize(out) // 1024, "KB")


if __name__ == "__main__":
    xyz, b = load_ca(PDB_AF)
    p = project2d(xyz)
    print("点云 %d 个 Cα，pLDDT 分档：>=90 %d / 70-90 %d / 50-70 %d / <50 %d"
          % (len(b), (b >= 90).sum(), ((b >= 70) & (b < 90)).sum(),
             ((b >= 50) & (b < 70)).sum(), (b < 50).sum()))
    cover_wechat(p, b)
    cover_xhs(p, b)
