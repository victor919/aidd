#!/usr/bin/env python
"""第 11 篇配图（8 张）

数据来源（全部本机实跑）：
  results/e9_split_bench.json      random / scaffold 全指标、域内域外、92 分子闭环
  results/e9_predictions.csv       逐分子预测（画散点）
  results/e9_domain.csv            最大 Tanimoto vs 绝对误差
  results/e10_leak_audit.json      随机划分分解（seen / unseen）、近重复、重复测量
  results/e11_noise_floor.json     噪声底与可达上限
  results/e13_leak_ladder.json     92 分子泄漏阶梯

纯示意（非数据）：图 1 / 图 2(左、中) / 图 3(左)，均在图上标"方法示意"。
字体：Hiragino Sans GB（缺 Å / ² / −，一律写"埃"、用 ASCII 减号）
"""
import os
import json
import io

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Circle

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
RES = os.path.join(HERE, "results")
FIG = os.path.join(ROOT, "AI制药第11篇_配图")
os.makedirs(FIG, exist_ok=True)

plt.rcParams["font.sans-serif"] = ["Hiragino Sans GB", "PingFang SC", "Arial Unicode MS"]
plt.rcParams["axes.unicode_minus"] = False

GREEN = "#0FA87E"
DARK = "#1A3C34"
RED = "#E53E3E"
AMBER = "#F5A623"
BLUE = "#4A7CFE"
GREY = "#9A9A9A"
PURPLE = "#7C5CBF"
LGREEN = "#E6F4F0"


def load(p):
    with open(p, encoding="utf-8") as f:
        return json.load(f)


def save(fig, name):
    p = os.path.join(FIG, name)
    fig.savefig(p, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", p, os.path.getsize(p) // 1024, "KB")


def box(ax, x, y, w, h, text, fc="#ffffff", ec=GREEN, fs=9.5, bold=False,
        tc=DARK, lw=1.4):
    ax.add_patch(FancyBboxPatch((x, y), w, h, transform=ax.transAxes,
                                boxstyle="round,pad=0.006,rounding_size=0.03",
                                fc=fc, ec=ec, lw=lw, zorder=3))
    ax.text(x + w / 2, y + h / 2, text, transform=ax.transAxes, ha="center",
            va="center", fontsize=fs, color=tc,
            fontweight="bold" if bold else "normal", zorder=4, linespacing=1.55)


def arrow(ax, x1, y1, x2, y2, color=DARK, lw=1.6):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), transform=ax.transAxes,
                                 arrowstyle="->", color=color, lw=lw,
                                 mutation_scale=13, zorder=2))


E9 = load(os.path.join(RES, "e9_split_bench.json"))
E10 = load(os.path.join(RES, "e10_leak_audit.json"))
E11 = load(os.path.join(RES, "e11_noise_floor.json"))
E13 = load(os.path.join(RES, "e13_leak_ladder.json"))
PRED = pd.read_csv(os.path.join(RES, "e9_predictions.csv"))
DOM = pd.read_csv(os.path.join(RES, "e9_domain.csv"))

R, S = E9["splits"]["random"], E9["splits"]["scaffold"]
DEC = E10["decompose_random"]
MAE_FLOOR = E11["ceilings"]["median"]["MAE_floor"]
R2_CAP = E11["ceilings"]["median"]["R2_cap"]
SIGMA = E11["sigma_median"]


# ---------------- 图 1 三站地图 ----------------
def fig1():
    fig, ax = plt.subplots(figsize=(9.2, 3.6))
    ax.axis("off")
    box(ax, 0.02, 0.28, 0.26, 0.5, "第 9 篇  结构与口袋\n\n研究对象准备好了吗",
        fc="#F7FBFA", ec=GREY, fs=10)
    box(ax, 0.37, 0.28, 0.26, 0.5, "第 10 篇  空间与对接\n\n分子怎么摆，评分怎么看",
        fc="#F7FBFA", ec=GREY, fs=10)
    box(ax, 0.72, 0.28, 0.26, 0.5, "第 11 篇  表征与预测\n\n模型怎样学，凭什么可信",
        fc=LGREEN, ec=GREEN, fs=10, bold=True)
    arrow(ax, 0.285, 0.53, 0.365, 0.53)
    arrow(ax, 0.635, 0.53, 0.715, 0.53)
    ax.text(0.5, 0.10, "三站同一个人物：厄洛替尼、同一个靶点：EGFR　｜　本篇从蛋白口径转回分子口径",
            transform=ax.transAxes, ha="center", fontsize=9, color=GREY)
    ax.set_title("从评价给定的分子，到预测分子的活性", fontsize=13, color=DARK, pad=14)
    save(fig, "图01_三站地图.png")


# ---------------- 图 2 三种表示 ----------------
def fig2():
    from rdkit import Chem
    from rdkit.Chem import Draw, rdFingerprintGenerator
    smi = "COCCOc1cc2ncnc(Nc3cccc(C#C)c3)c2cc1OCCOC"
    mol = Chem.MolFromSmiles(smi)

    fig = plt.figure(figsize=(9.4, 4.4))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 0.85], hspace=0.35, wspace=0.22)

    ax0 = fig.add_subplot(gs[0, 0])
    ax0.axis("off")
    ax0.set_title("① SMILES：一行文本", fontsize=11, color=DARK)
    txt = "COCCOc1cc2ncnc\n(Nc3cccc(C#C)c3)\nc2cc1OCCOC"
    ax0.text(0.02, 0.55, txt, family="Menlo", fontsize=8.2, color="#C0392B",
             va="center", transform=ax0.transAxes, linespacing=1.7)
    ax0.text(0.02, 0.05, "字符顺序即连接顺序\n换一个字符就是另一个分子",
             fontsize=8.5, color=GREY, transform=ax0.transAxes)

    ax1 = fig.add_subplot(gs[0, 1])
    ax1.axis("off")
    ax1.set_title("② 分子图：原子与键", fontsize=11, color=DARK)
    d = Draw.MolDraw2DCairo(560, 380)
    Draw.rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
    d.FinishDrawing()
    ax1.imshow(plt.imread(io.BytesIO(d.GetDrawingText()), format="png"))

    ax2 = fig.add_subplot(gs[0, 2])
    ax2.axis("off")
    ax2.set_title("③ ECFP4：2048 位里的一行 0/1", fontsize=11, color=DARK)
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fp = gen.GetFingerprint(mol)
    nz = np.array(fp.GetOnBits(), dtype=int)
    grid = np.zeros((32, 64))
    for b in nz:
        grid[b // 64, b % 64] = 1
    ax2.imshow(grid, cmap="Greens", aspect="auto", vmin=0, vmax=1.3)
    ax2.text(0.5, -0.12, f"2048 位中 {len(nz)} 位为 1（示意图，每格 32 位）",
             fontsize=8.5, color=GREY, ha="center", transform=ax2.transAxes)

    ax3 = fig.add_subplot(gs[1, :])
    ax3.axis("off")
    ax3.text(0.0, 0.95, "三种表示各自丢掉了什么", fontsize=11, color=DARK, fontweight="bold",
             transform=ax3.transAxes, va="top")
    rows = [
        ("SMILES", "丢了图结构写成文本的规则性，需解析才能算相似", GREY),
        ("分子图", "保留连接与手性，但位置、朝向仍未编码（⑩ 讲的等变问题）", BLUE),
        ("ECFP4", "只留局部子结构是否出现；丢三维、丢手性、丢绝对位置", GREEN),
    ]
    for i, (a, b, c) in enumerate(rows):
        y = 0.62 - i * 0.28
        ax3.text(0.0, y, a, fontsize=9.5, color=c, fontweight="bold", transform=ax3.transAxes)
        ax3.text(0.13, y, b, fontsize=9, color="#3F3F3F", transform=ax3.transAxes)
    save(fig, "图02_三种表示.png")


# ---------------- 图 3 半径与手性 ----------------
def fig3():
    from rdkit import Chem, DataStructs
    from rdkit.Chem import rdFingerprintGenerator

    # 厄洛替尼本身没有手性中心（镜像即自身，相似度恒为 1，不能用作例子）。
    # 这里用一对真实的对映体：布洛芬 R / S。
    s_smi = "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O"
    r_smi = "CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O"

    def fp(smi, chir):
        g = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048,
                                                      includeChirality=chir)
        return g.GetFingerprint(Chem.MolFromSmiles(smi))

    t_no = DataStructs.TanimotoSimilarity(fp(s_smi, False), fp(r_smi, False))
    t_yes = DataStructs.TanimotoSimilarity(fp(s_smi, True), fp(r_smi, True))
    print(f"  布洛芬 R/S Tanimoto: 不含手性 {t_no:.4f} / 含手性 {t_yes:.4f}")

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.9), gridspec_kw={"width_ratios": [1.25, 1]})

    ax = axes[0]
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("ECFP4 只看半径 2 以内的邻域（方法示意）", fontsize=11, color=DARK)
    cx, cy = 0.5, 0.48
    for rr, lab, cc in [(0.30, "r = 2：两跳", "#CFE9E1"), (0.18, "r = 1：一跳", "#A8DCCD")]:
        ax.add_patch(Circle((cx, cy), rr, transform=ax.transAxes, fc=cc, ec=GREEN, lw=1.2, alpha=.55))
    ax.plot([cx], [cy], marker="o", ms=7, color=DARK, transform=ax.transAxes)
    ax.plot([cx + 0.13], [cy + 0.06], marker="o", ms=6, color=BLUE, transform=ax.transAxes)
    ax.plot([cx + 0.24], [cy + 0.19], marker="o", ms=6, color=AMBER, transform=ax.transAxes)
    ax.plot([cx - 0.16], [cy - 0.12], marker="o", ms=6, color=BLUE, transform=ax.transAxes)
    ax.plot([cx - 0.27], [cy - 0.25], marker="o", ms=6, color=GREY, transform=ax.transAxes)
    ax.plot([cx, cx + 0.13], [cy, cy + 0.06], color=DARK, lw=1.2, transform=ax.transAxes)
    ax.plot([cx + 0.13, cx + 0.24], [cy + 0.06, cy + 0.19], color=DARK, lw=1.2, transform=ax.transAxes)
    ax.plot([cx, cx - 0.16], [cy, cy - 0.12], color=DARK, lw=1.2, transform=ax.transAxes)
    ax.plot([cx - 0.16, cx - 0.27], [cy - 0.12, cy - 0.25], color=GREY, lw=1.2,
            linestyle="--", transform=ax.transAxes)
    ax.text(cx - 0.30, cy - 0.34, "r = 3 以外：看不见", fontsize=8.5, color=GREY,
            transform=ax.transAxes)
    ax.text(0.5, 0.02, "同一个原子，不同半径下被编码进不同的位", fontsize=8.5,
            color=GREY, ha="center", transform=ax.transAxes)

    ax = axes[1]
    ax.set_title("一对真实对映体：默认配置分不出左右手", fontsize=11, color=DARK)
    bars = ax.bar(["ECFP4\n默认", "ECFP4\n含手性"], [t_no, t_yes],
                  color=[RED, GREEN], width=0.5, edgecolor="white")
    for b, v in zip(bars, [t_no, t_yes]):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.02, f"{v:.3f}", ha="center",
                fontsize=11, color=DARK, fontweight="bold")
    ax.set_ylim(0, 1.12)
    ax.set_ylabel("布洛芬 R 型 vs S 型的 Tanimoto", fontsize=9.5, color="#3F3F3F")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=9)
    ax.text(0.5, 0.42, "默认配置下相似度 = 1.000\n两个分子对模型是同一个", fontsize=9,
            color=RED, ha="center", transform=ax.transAxes)
    save(fig, "图03_半径与手性.png")


# ---------------- 图 4 预测 vs 真实散点 ----------------
def fig4():
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.4), sharex=True, sharey=True)
    for ax, name, title, col in [
        (axes[0], "random", "随机划分", BLUE),
        (axes[1], "scaffold", "骨架划分", AMBER)]:
        d = PRED[PRED["split"] == name]
        ax.scatter(d["true_pIC50"], d["pred_pIC50"], s=9, alpha=.35, color=col,
                   edgecolors="none")
        lim = [1, 11.5]
        ax.plot(lim, lim, color=GREY, lw=1, linestyle="--")
        ax.fill_between(lim, [l - SIGMA for l in lim], [l + SIGMA for l in lim],
                        color=GREY, alpha=.12, lw=0)
        m = E9["splits"][name]
        ax.set_title(f"{title}　MAE {m['test_MAE']:.3f}　R2 {m['test_R2']:.3f}",
                     fontsize=10.5, color=DARK)
        ax.set_xlim(1, 11.5)
        ax.set_ylim(1, 11.5)
        ax.set_xlabel("实验 pIC50", fontsize=9.5, color="#3F3F3F")
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=9)
        ax.set_aspect("equal")
    axes[0].set_ylabel("预测 pIC50", fontsize=9.5, color="#3F3F3F")
    axes[0].text(1.4, 10.4, "灰带 = ±0.43\n（重复测量中位标准差）", fontsize=8.5, color=GREY)
    fig.suptitle("同一个模型、同一套指纹，两种划分的测试集表现", fontsize=12.5,
                 color=DARK, y=1.02)
    save(fig, "图04_预测vs真实.png")


# ---------------- 图 5 核心对照 ----------------
def fig5():
    labels = ["随机划分\n全部 2093", "其中骨架见过\n1536 个",
              "其中没见过\n557 个", "骨架划分\n全部 1770"]
    maes = [DEC["all_MAE"], DEC["seen_MAE"], DEC["unseen_MAE"], S["test_MAE"]]
    r2s = [DEC["all_R2"], DEC["seen_R2"], DEC["unseen_R2"], S["test_R2"]]
    cols = [BLUE, "#7FA8FF", "#9FBDFF", AMBER]

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.2))
    ax = axes[0]
    bars = ax.bar(range(4), maes, color=cols, width=.62, edgecolor="white")
    for b, v in zip(bars, maes):
        ax.text(b.get_x() + b.get_width() / 2, v + .012, f"{v:.3f}", ha="center",
                fontsize=10, color=DARK, fontweight="bold")
    ax.axhline(MAE_FLOOR, color=RED, lw=1.2, linestyle="--")
    ax.text(0.02, MAE_FLOOR + .015, f"噪声底 {MAE_FLOOR:.2f}", fontsize=8.5,
            color=RED, ha="left", transform=ax.get_yaxis_transform())
    ax.set_xticks(range(4))
    ax.set_xticklabels(labels, fontsize=8.2)
    ax.set_ylabel("测试集 MAE（pIC50）", fontsize=9.5)
    ax.set_ylim(0, .72)
    ax.set_title("误差：把随机划分拆开看", fontsize=11, color=DARK)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8.5)

    ax = axes[1]
    bars = ax.bar(range(4), r2s, color=cols, width=.62, edgecolor="white")
    for b, v in zip(bars, r2s):
        ax.text(b.get_x() + b.get_width() / 2, v + .015, f"{v:.3f}", ha="center",
                fontsize=10, color=DARK, fontweight="bold")
    ax.axhline(R2_CAP, color=RED, lw=1.2, linestyle="--")
    ax.text(0.02, R2_CAP - .045, f"噪声上限 {R2_CAP:.2f}", fontsize=8.5, color=RED, ha="left",
            transform=ax.get_yaxis_transform())
    ax.set_xticks(range(4))
    ax.set_xticklabels(labels, fontsize=8.2)
    ax.set_ylabel("测试集 R2", fontsize=9.5)
    ax.set_ylim(0, 1.0)
    ax.set_title("R2：跌下去的那截是什么", fontsize=11, color=DARK)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8.5)

    fig.suptitle("随机划分的 0.509 里，藏着 0.595 的另一半", fontsize=12.5, color=DARK, y=1.01)
    save(fig, "图05_核心对照.png")


# ---------------- 图 6 相似度分布 ----------------
def fig6():
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.9))
    ax = axes[0]
    bins = np.linspace(0, 1, 26)
    for name, col, lab in [("random", BLUE, "随机划分"), ("scaffold", AMBER, "骨架划分")]:
        v = DOM[DOM["split"] == name]["maxsim"]
        ax.hist(v, bins=bins, alpha=.55, color=col, label=lab, density=True)
    ax.axvline(0.7, color=GREY, lw=1, linestyle=":")
    ax.text(0.705, ax.get_ylim()[1] * .92, "0.7", fontsize=8.5, color=GREY)
    ax.set_xlabel("测试分子到训练集的最近邻 Tanimoto", fontsize=9.5)
    ax.set_ylabel("密度", fontsize=9.5)
    ax.legend(fontsize=9, frameon=False)
    ax.set_title("随机划分更容易撞上“几乎一样”的分子", fontsize=11, color=DARK)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8.5)

    ax = axes[1]
    d = DOM.copy()
    grp = [d[(d["split"] == "random") & (d["maxsim"] >= .7)]["abs_err"],
           d[(d["split"] == "random") & (d["maxsim"] < .7)]["abs_err"],
           d[(d["split"] == "scaffold") & (d["maxsim"] >= .7)]["abs_err"],
           d[(d["split"] == "scaffold") & (d["maxsim"] < .7)]["abs_err"]]
    labs = ["随机\n近邻≥0.7", "随机\n近邻<0.7", "骨架\n近邻≥0.7", "骨架\n近邻<0.7"]
    bp = ax.boxplot(grp, tick_labels=labs, showfliers=False, patch_artist=True, widths=.55)
    for patch, c in zip(bp["boxes"], [BLUE, "#B9CDFF", AMBER, "#FFD79A"]):
        patch.set_facecolor(c)
        patch.set_edgecolor("white")
    for med in bp["medians"]:
        med.set_color(DARK)
        med.set_linewidth(1.4)
    ax.set_ylabel("绝对误差（pIC50）", fontsize=9.5)
    ax.set_title("身边有没有“近亲”，决定误差大小", fontsize=11, color=DARK)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8.2)
    save(fig, "图06_相似度与误差.png")


# ---------------- 图 7 四把尺 ----------------
def fig7():
    panels = [
        ("MAE", [R["test_MAE"], S["test_MAE"]], "{:.3f}", "越低越好"),
        ("R2", [R["test_R2"], S["test_R2"]], "{:.3f}", "越高越好"),
        ("Spearman ρ", [R["spearman"], S["spearman"]], "{:.3f}", "越高越好"),
        ("EF@5%", [R["EF5"], S["EF5"]], "{:.2f}", "越高越好"),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(9.6, 3.5))
    for ax, (name, vals, fmt, hint) in zip(axes, panels):
        bars = ax.bar([0, 1], vals, color=[BLUE, AMBER], width=.62, edgecolor="white")
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, v * 1.02, fmt.format(v), ha="center",
                    fontsize=10, color=DARK, fontweight="bold")
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["随机", "骨架"], fontsize=9)
        ax.set_title(name, fontsize=11, color=DARK)
        ax.set_ylim(0, max(vals) * 1.22)
        ax.text(0.5, -0.26, hint, transform=ax.transAxes, ha="center", fontsize=8,
                color=GREY)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=8.5)
    fig.suptitle("同一个模型，四把尺给出的四个故事", fontsize=12.5, color=DARK, y=1.05)
    save(fig, "图07_四把尺.png")


# ---------------- 图 8 泄漏阶梯 ----------------
def fig8():
    lad = E13["ladder"]
    fig, ax = plt.subplots(figsize=(9.2, 4.2))
    xs = [10, 20, 30]
    ax.plot(xs, [E13["ladder"][0]["recall_vina"][str(x)] for x in xs], marker="o",
            ms=6, lw=2, color=GREY, label="Vina 对接分数（⑩）")
    cols = [RED, "#F08A8A", AMBER]
    for r, c in zip(lad, cols):
        ax.plot(xs, [r["recall"][str(x)] for x in xs], marker="o", ms=6, lw=2, color=c,
                label=f"QSAR · {r['tag']}")
    ax.plot(xs, [x * 30 / 92 for x in xs], linestyle=":", color="#BBBBBB", lw=1.4,
            label="随机期望")
    ax.set_xticks(xs)
    ax.set_xlabel("预算（挑前 k 个送去验证）", fontsize=9.5)
    ax.set_ylabel("命中的真活性分子数（目标 = 真实 Top-30）", fontsize=9.5)
    ax.set_ylim(0, 31)
    ax.legend(fontsize=8.4, frameon=False, loc="lower right")
    ax.set_title("同 92 个分子：把训练集里的“近亲”一层层删掉，排序还剩多少",
                 fontsize=11.5, color=DARK)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", alpha=.25)
    ax.tick_params(labelsize=9)
    save(fig, "图08_泄漏阶梯.png")


if __name__ == "__main__":
    fig1()
    fig2()
    fig3()
    fig4()
    fig5()
    fig6()
    fig7()
    fig8()
