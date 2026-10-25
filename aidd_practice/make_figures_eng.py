# 文章关联：AI制药工程篇_公众号文章.html
# 脚本职责：配图：e1 工程结构 / e2 坑与断言
# 使用注解：真正的流水线、预检和自检在 ../aidd-mvp/；本脚本只绘制说明图。
# 完整对应关系与运行顺序：README.md
"""工程篇配图：
e1_structure.png  工程结构：11 个脚本 + 三层数据区 + 两端守门
e2_guard.png      坑 -> 断言 映射：每一个断言都对应一个真实踩过的坑
"""
import os
import warnings

warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

plt.rcParams["font.family"] = "Hiragino Sans GB"
plt.rcParams["axes.unicode_minus"] = False

GREEN, RED, BLUE = "#0FA87E", "#C0392B", "#2B5F8F"
GREY, AMBER = "#9A9A9A", "#E8A33D"
FIG = os.path.expanduser("~/WorkBuddy/Work/AI制药工程篇_配图")
os.makedirs(FIG, exist_ok=True)

# ============================================================ 图 1：工程结构
fig, ax = plt.subplots(figsize=(9.2, 6.8))
ax.set_xlim(0, 11.1); ax.set_ylim(0, 10); ax.axis("off")

STEPS = [
    ("01", "靶点排序"), ("02", "下载结构"), ("03", "口袋+编号"), ("04", "拉活性"),
    ("05", "清洗聚合"), ("06", "切分+指纹"), ("07", "RF 基线"), ("08", "D-MPNN"),
    ("09", "对接自证"), ("10", "成药性"), ("11", "闭环"),
]
BW, BH = 1.35, 1.05
rows = [STEPS[:6], STEPS[6:]]
X0, DX = 1.75, 1.55
y0 = 8.4
for ri, row in enumerate(rows):
    y = y0 - ri * 1.75
    for i, (num, name) in enumerate(row):
        x = X0 + i * DX
        ax.add_patch(Rectangle((x, y), BW, BH, fc="#eaf6f1", ec=GREEN, lw=1.6, zorder=3))
        ax.text(x + BW / 2, y + BH / 2 + 0.19, num, fontsize=12, fontweight="bold",
                color=GREEN, ha="center", va="center", zorder=4)
        ax.text(x + BW / 2, y + BH / 2 - 0.32, name, fontsize=9.5, color="#333",
                ha="center", va="center", zorder=4)
        if i < len(row) - 1:
            ax.annotate("", xy=(x + BW + 0.14, y + BH / 2), xytext=(x + BW + 0.01, y + BH / 2),
                        arrowprops=dict(arrowstyle="-|>", color=GREEN, lw=1.4))
    if ri == 0:
        # 换行连接：首列虚线向下
        ax.annotate("", xy=(X0 + BW / 2, y - 0.64), xytext=(X0 + BW / 2, y - 0.02),
                    arrowprops=dict(arrowstyle="-|>", color=GREEN, lw=1.4, ls=":"))

# 两端守门：00 在第一行左端，99 在第二行右端
ax.add_patch(Rectangle((0.05, y0), 1.3, BH, fc="#eaf0fb", ec=BLUE, lw=1.6, zorder=3))
ax.text(0.05 + 0.65, y0 + BH / 2, "00\n预检", fontsize=11, fontweight="bold", color=BLUE,
        ha="center", va="center", zorder=4, linespacing=1.3)
ax.annotate("", xy=(X0 - 0.06, y0 + BH / 2), xytext=(1.38, y0 + BH / 2),
            arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=1.5))

y2 = y0 - 1.75
ax.add_patch(Rectangle((9.65, y2), 1.3, BH, fc="#eaf0fb", ec=BLUE, lw=1.6, zorder=3))
ax.text(9.65 + 0.65, y2 + BH / 2, "99\n自检", fontsize=11, fontweight="bold", color=BLUE,
        ha="center", va="center", zorder=4, linespacing=1.3)
ax.annotate("", xy=(9.62, y2 + BH / 2), xytext=(X0 + 4 * DX + BW + 0.01, y2 + BH / 2),
            arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=1.5))

ax.text(0.0, 5.95, "数据三层（目录即契约，谁写谁负责）", fontsize=11.5, fontweight="bold", color="#333")
bands = [
    ("data/raw/", "原始下载物，永不手改", "#eef1f4", "#8a9299"),
    ("data/interim/", "中间产物：pdbqt / pose / 盒子", "#e8f2ee", GREEN),
    ("data/processed/", "可交付：表格、指标、结论", "#fdf3e3", AMBER),
]
for i, (name, desc, fc, ec) in enumerate(bands):
    y = 4.95 - i * 1.15
    ax.add_patch(Rectangle((0.05, y), 10.9, 0.92, fc=fc, ec=ec, lw=1.5, zorder=3))
    ax.text(0.22, y + 0.46, name, fontsize=11, fontweight="bold", color="#333",
            va="center", zorder=4, family="monospace")
    ax.text(3.05, y + 0.46, desc, fontsize=10.5, color="#555", va="center", zorder=4)

ax.text(0.0, 1.35, "config.yaml：唯一参数源（脚本内禁止硬编码 P00533 / 2048 / 42）",
        fontsize=11, color=RED, fontweight="bold")
ax.text(0.0, 0.55, "run_all.sh：预检 → 01-11 → 报告 → 产物自洽 → 不变量自检（任一 FAIL 即非零退出）",
        fontsize=10.5, color="#555")
ax.set_title("aidd-mvp：把 11 个脚本串成一条能复跑、能自检的流水线", fontsize=12.5,
             fontweight="bold", color="#1a1a1a", pad=12, loc="left")
fig.tight_layout()
fig.savefig(f"{FIG}/e1_structure.png", dpi=170, bbox_inches="tight", facecolor="#fff")
plt.close(fig)
print("e1_structure.png")

# ============================================================ 图 2：坑 -> 断言
fig, ax = plt.subplots(figsize=(9.2, 7.0))
ax.set_xlim(0, 10); ax.set_ylim(0, 11.4); ax.axis("off")

PAIRS = [
    ("编号差 24（C773≠C797）", "C1/C2/C3 偏移自检测 + 序列匹配 100%"),
    ("口袋坐标写死", "C4 盒中心 == 配体质心（偏差 0.000）"),
    ("受体里混进共晶配体", "D1 晶体自打分 < 0（-6.02，异常时 +438）"),
    ("对接没做自证", "D2 重对接 RMSD < 2 埃（1.411）"),
    ("随机切分导致「默写」", "B1/B2 scaffold 零重叠 + 对照组有重叠"),
    ("脏数据 pIC50 = 17.3", "A3 上界 ≤ 11（截断确已生效）"),
    ("chemprop 归一化写反", "E1/E2 MAE 区间 + 切分方向断言"),
    ("docking 与活性正交", "E3 |r| < 0.5（0.059）"),
    ("greedy 六轮原地踏步", "E4 终点 uncertainty ≤ greedy"),
    ("参数漂移 / 结果对不上", "F1/F2 重算指纹与切分，逐位一致"),
]
ax.text(0.1, 10.95, "踩过的坑", fontsize=12, fontweight="bold", color=RED, ha="center")
ax.text(9.9, 10.95, "锁住它的断言", fontsize=12, fontweight="bold", color=GREEN, ha="center")

y = 10.3
for i, (pit, guard) in enumerate(PAIRS):
    # 左侧坑
    ax.add_patch(Rectangle((0.1, y - 0.42), 4.25, 0.84, fc="#fdecec", ec="#e0b4b4", lw=1.2, zorder=3))
    ax.text(0.32, y, pit, fontsize=9.8, color="#8c2b2b", va="center", zorder=4)
    # 右侧断言
    ax.add_patch(Rectangle((5.55, y - 0.42), 4.35, 0.84, fc="#eaf6f1", ec=GREEN, lw=1.2, zorder=3))
    ax.text(5.77, y, guard, fontsize=9.5, color="#14624c", va="center", zorder=4)
    # 连线
    ax.annotate("", xy=(5.5, y), xytext=(4.4, y),
                arrowprops=dict(arrowstyle="-|>", color=GREY, lw=1.2, ls=":"))
    y -= 1.02

ax.text(5.0, 0.28, "文档会过时，注释会被忽略，只有断言会被执行：19 项快速自检 / 21 项含复现性重算",
        fontsize=10.5, color="#555", ha="center")
ax.set_title("99_selftest.py：把坑固化成断言", fontsize=12.5, fontweight="bold",
             color="#1a1a1a", pad=10, loc="left")
fig.tight_layout()
fig.savefig(f"{FIG}/e2_guard.png", dpi=170, bbox_inches="tight", facecolor="#fff")
plt.close(fig)
print("e2_guard.png")
