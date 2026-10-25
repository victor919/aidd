# 文章关联：第4课素材（对应正文尚未在相邻目录找到）
# 脚本职责：第4课配图：盲区、排序、盒子与工具定位
# 使用注解：读取 D1 / D2 / D3 JSON，输出本项目 AI制药第4课_配图/；相邻目录未找到对应第4课正文。
# 完整对应关系与运行顺序：../README.md
"""第 4 课配图 · 全部基于本机实跑结果

图 1｜2D 指纹的三个盲区（手性 / 构象 / 2D-3D 排序不一致）
图 2｜Vina 分数与真实活性的关系：分箱单调性 + 有限预算下的召回
图 3｜盒子瞄哪答案就是哪：三盒对比（正确口袋 / 偏移 8 埃 / 全蛋白盲搜）
图 4｜Vina vs DiffDock 的定位差异：横轴"需要预先定义口袋"，纵轴"分数含义"

数据来源：
  results/d1_blindspots.json
  results/d2_vina_real_activity.json
  results/d3_box_sensitivity.json

字体：Hiragino Sans GB（注意该字体缺 Å / ² / −，一律用"埃"/ASCII）
"""
import os
import sys
import json
import warnings

warnings.filterwarnings("ignore")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
FIG = os.path.join(os.path.dirname(HERE), "AI制药第4课_配图")
os.makedirs(FIG, exist_ok=True)

plt.rcParams["font.sans-serif"] = ["Hiragino Sans GB", "PingFang SC", "Arial Unicode MS"]
plt.rcParams["axes.unicode_minus"] = False

GREEN = "#0FA87E"
DARK = "#1a3c34"
RED = "#E53E3E"
AMBER = "#F5A623"
BLUE = "#4A7CFE"
GREY = "#9A9A9A"

# matplotlib 在非等比例坐标轴下会把 Circle 拉成椭圆，标记一律用 plot(marker="o")


def load(name):
    with open(os.path.join(RES, name), encoding="utf-8") as f:
        return json.load(f)


# ------------------------------------------------------------------ 图 1
def fig1():
    d1 = load("d1_blindspots.json")
    b1 = d1["blindspot_1_chirality"]
    b2 = d1["blindspot_2_conformation"]
    b3 = d1["blindspot_3_shape_vs_topology"]

    fig = plt.figure(figsize=(9.2, 9.6))
    gs = fig.add_gridspec(4, 2, height_ratios=[1.05, 1.0, 1.15, 0.95],
                          hspace=0.52, wspace=0.30,
                          left=0.085, right=0.965, top=0.935, bottom=0.055)

    # (a) 手性：ECFP4 对对映体视而不见
    ax = fig.add_subplot(gs[0, 0])
    names = list(b1.keys())
    labels = ["沙利度胺类\n(R/S 对映体)", "乳酸\n(D/L 对映体)"]
    no_ch = [b1[n]["sim_ecfp4_no_chirality"] for n in names]
    wi_ch = [b1[n]["sim_ecfp4_with_chirality"] for n in names]
    x = np.arange(len(names))
    ax.bar(x - 0.18, no_ch, 0.34, label="ECFP4（默认，不含手性）", color=GREY, zorder=3)
    ax.bar(x + 0.18, wi_ch, 0.34, label="ECFP4（含手性）", color=GREEN, zorder=3)
    for i, (a, b) in enumerate(zip(no_ch, wi_ch)):
        ax.text(i - 0.18, a + 0.015, f"{a:.4f}", ha="center", fontsize=9, color="#555")
        ax.text(i + 0.18, b + 0.015, f"{b:.4f}", ha="center", fontsize=9, color=DARK)
    ax.set_ylim(0, 1.16)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9.5)
    ax.set_ylabel("Tanimoto 相似度", fontsize=10)
    ax.set_title("(a) 对映体在 2D 指纹里长得一模一样", fontsize=11.5, color=DARK, pad=8)
    ax.legend(fontsize=8.5, loc="lower left", frameon=False)
    ax.grid(axis="y", alpha=0.25, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)

    # (b) 构象：同一分子 14 个构象，2D 恒为 1
    ax = fig.add_subplot(gs[0, 1])
    ax.bar([0], [b2["sim_ecfp4_same_molecule"]], 0.42,
           label="ECFP4 相似度（任何构象）", color=GREY, zorder=3)
    ax.bar([1], [b2["shape_usrcat_two_independent_confs"]], 0.42,
           label="USRCAT 3D 形状相似度", color=BLUE, zorder=3)
    ax.text(0, 1.0 + 0.015, "1.0000", ha="center", fontsize=9, color="#555")
    ax.text(1, b2["shape_usrcat_two_independent_confs"] + 0.015,
            f'{b2["shape_usrcat_two_independent_confs"]:.4f}', ha="center", fontsize=9, color=DARK)
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylim(0, 1.16)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["同一分子\n不同构象", "同一分子\n不同构象"], fontsize=9.5)
    ax.set_ylabel("相似度", fontsize=10)
    ax.set_title(f'(b) 厄洛替尼 {b2["n_conformers_embedded"]} 个构象，'
                 f'最大构象对 RMSD {b2["conformer_pair_rmsd_max"]} 埃', fontsize=11.5, color=DARK, pad=8)
    ax.legend(fontsize=8.5, loc="lower left", frameon=False)
    ax.grid(axis="y", alpha=0.25, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)

    # (c) 2D vs 3D 排序不一致：散点
    ax = fig.add_subplot(gs[1, :])
    pairs = b3["pairs"]
    xs = [p["sim_2d_ecfp4"] for p in pairs]
    ys = [p["sim_3d_usrcat"] for p in pairs]
    ax.plot(xs, ys, "o", ms=9, color=GREEN, alpha=0.85, zorder=4)
    worst = b3["most_disagreeing_pair"]
    for p in pairs:
        if p["pair"] == worst["pair"]:
            ax.annotate(f'{p["pair"]}\n2D={p["sim_2d_ecfp4"]:.3f} / 3D={p["sim_3d_usrcat"]:.3f}',
                        (p["sim_2d_ecfp4"], p["sim_3d_usrcat"]),
                        textcoords="offset points", xytext=(18, 20),
                        fontsize=9, color=DARK,
                        arrowprops=dict(arrowstyle="->", color=AMBER, lw=1.4))
    ax.axhline(float(np.mean(ys)), color=GREY, ls="--", lw=1.1, zorder=1)
    ax.text(0.02, 0.06, f'10 对 EGFR 药物\n2D 与 3D 相似度的 Spearman ρ = '
                       f'{b3["spearman_rho_2d_vs_3d"]:.4f}\n（几乎完全不相关）',
            transform=ax.transAxes, fontsize=10.5, color=RED, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.45", fc="#fff5f5", ec=RED, lw=1.2))
    ax.set_xlabel("2D 拓扑相似度（ECFP4 Tanimoto）", fontsize=10.5)
    ax.set_ylabel("3D 形状相似度（USRCAT）", fontsize=10.5)
    ax.set_title("(c) 2D 拓扑和 3D 形状是两套排序——它们几乎不 disagree",
                 fontsize=11.5, color=DARK, pad=8)
    ax.set_xlim(-0.02, max(xs) + 0.12)
    ax.set_ylim(-0.02, max(ys) + 0.10)
    ax.grid(alpha=0.25, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)

    # (d) 盲区三连：文字卡
    ax = fig.add_subplot(gs[2, :])
    ax.axis("off")
    cards = [
        ("盲区 1｜手性", "MVP 的 Morgan 生成器没开 includeChirality",
         f'两个对映体 ECFP4 相似度 = {b1["lactic_acid_D_vs_L"]["sim_ecfp4_no_chirality"]:.4f}\n'
         f'打开手性后降到 {b1["lactic_acid_D_vs_L"]["sim_ecfp4_with_chirality"]:.4f}',
         RED),
        ("盲区 2｜构象", "2D 指纹根本不看坐标",
         f'厄洛替尼 {b2["n_conformers_embedded"]} 个构象，2D 恒为 1.0000\n'
         f'但构象对 RMSD 可达 {b2["conformer_pair_rmsd_max"]} 埃（{b2["rotatable_bonds"]} 个可旋转键）',
         AMBER),
        ("盲区 3｜形状", "拓扑不等于空间",
         f'2D 与 3D 相似度 Spearman ρ = {b3["spearman_rho_2d_vs_3d"]:.4f}\n'
         f'最不一致：{worst["pair"]}（差 {worst["diff_2d_minus_3d"]:+.3f}）',
         BLUE),
    ]
    for i, (title, sub, body, color) in enumerate(cards):
        x0 = 0.015 + i * 0.325
        ax.add_patch(plt.Rectangle((x0, 0.04), 0.295, 0.92, transform=ax.transAxes,
                                   fc="#fbfcfc", ec=color, lw=1.6, zorder=2))
        ax.text(x0 + 0.015, 0.87, title, transform=ax.transAxes, fontsize=11.5,
                color=color, fontweight="bold", va="top", zorder=3)
        ax.text(x0 + 0.015, 0.70, sub, transform=ax.transAxes, fontsize=9.3,
                color="#666", va="top", zorder=3)
        ax.text(x0 + 0.015, 0.52, body, transform=ax.transAxes, fontsize=9.6,
                color=DARK, va="top", zorder=3, linespacing=1.65)

    # (e) 结论条
    ax = fig.add_subplot(gs[3, :])
    ax.axis("off")
    ax.add_patch(plt.Rectangle((0.01, 0.12), 0.98, 0.74, transform=ax.transAxes,
                               fc="#f0f9f5", ec=GREEN, lw=1.5, zorder=1))
    ax.text(0.035, 0.70, "一句话结论", transform=ax.transAxes, fontsize=11.5,
            color=GREEN, fontweight="bold", va="top", zorder=3)
    ax.text(0.035, 0.50,
            "ECFP4 不是“不准确”，而是“信息维度有限”：它精确描述了拓扑，"
            "对拓扑相同但空间不同的分子一律给满分。\n"
            "所以 2D 模型的错误不是算错了，而是它按设计就看不见那些维度——"
            "这不是换模型能解决的，得换表征。",
            transform=ax.transAxes, fontsize=10.3, color=DARK, va="top", zorder=3,
            linespacing=1.8)

    out = os.path.join(FIG, "c1_blindspots.png")
    fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("写出", out)


# ------------------------------------------------------------------ 图 2
def fig2():
    d2 = load("d2_vina_real_activity.json")
    corr = d2["corr_vina_vs_true_pIC50"]
    bins = d2["vina_bins"]
    rec = d2["budget_recall_of_top30"]["strategies"]

    fig = plt.figure(figsize=(9.2, 7.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1.0], hspace=0.46, wspace=0.26,
                          left=0.085, right=0.975, top=0.925, bottom=0.085)

    # (a) 分箱单调性
    ax = fig.add_subplot(gs[0, :])
    labels = [f"{b['vina_range'][0]:.1f} ~ {b['vina_range'][1]:.1f}" for b in bins]
    means = [b["mean_true_pIC50"] for b in bins]
    meds = [b["median_true_pIC50"] for b in bins]
    xs = np.arange(len(bins))
    ax.bar(xs, means, 0.5, color=GREEN, alpha=0.88, zorder=3, label="平均真实 pIC50")
    ax.plot(xs, meds, "o-", color=DARK, lw=1.8, ms=7, zorder=5, label="中位真实 pIC50")
    for x, m in zip(xs, means):
        ax.text(x, m + 0.16, f"{m:.2f}", ha="center", fontsize=10, color=DARK)
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, fontsize=9.3)
    ax.set_xlabel("Vina 最优分数区间（kcal/mol，越负越强）", fontsize=10.5)
    ax.set_ylabel("真实 pIC50", fontsize=10.5)
    ax.set_title(f'(a) {d2["n_docked"]} 个分子按 Vina 分数分 5 箱，箱内真实活性单调下降'
                 f'（Spearman ρ = {corr["spearman_rho"]:.4f}, p = {corr["spearman_p"]:.2g}）',
                 fontsize=11.3, color=DARK, pad=9)
    ax.legend(fontsize=9.3, loc="upper right", frameon=False)
    ax.grid(axis="y", alpha=0.25, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)

    # (b) 预算召回
    ax = fig.add_subplot(gs[1, 0])
    budgets = [10, 20, 30]
    order = ["A_只按Vina", "B_只按真实活性", "C_consensus"]
    nice = {"A_只按Vina": "只按 Vina 分数", "B_只按真实活性": "按真实活性（理想上界）",
            "C_consensus": "consensus 混合"}
    colors = {"A_只按Vina": RED, "B_只按真实活性": GREY, "C_consensus": GREEN}
    for k in order:
        vals = [rec[k][str(b)] for b in budgets]
        ax.plot(budgets, vals, "o-", lw=2.0, ms=7, color=colors[k],
                label=nice[k], zorder=4)
        for b, v in zip(budgets, vals):
            ax.text(b, v + 0.9, str(v), ha="center", fontsize=9.2, color=colors[k])
    ax.set_xticks(budgets)
    ax.set_xlabel("实验预算（能送去做实验的分子数）", fontsize=10.2)
    ax.set_ylabel("命中真实 Top-30 的个数", fontsize=10.2)
    ax.set_ylim(0, 33)
    ax.set_title("(b) 有限预算下谁能选出真药", fontsize=11.3, color=DARK, pad=9)
    ax.legend(fontsize=9, loc="lower right", frameon=False)
    ax.grid(alpha=0.25, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)

    # (c) 结论
    ax = fig.add_subplot(gs[1, 1])
    ax.axis("off")
    ax.add_patch(plt.Rectangle((0.02, 0.06), 0.96, 0.88, transform=ax.transAxes,
                               fc="#fffdf5", ec=AMBER, lw=1.6, zorder=1))
    ax.text(0.06, 0.82, "读数要点", transform=ax.transAxes, fontsize=11.5,
            color="#B87503", fontweight="bold", va="top", zorder=3)
    body = (
        f'· Vina 分数确实有信号：ρ = {corr["spearman_rho"]:.4f}，p = {corr["spearman_p"]:.2g}\n'
        f'· 但 |ρ| 只有 0.5，属于“中等”，远不是“能靠它选药”\n'
        f'· 预算 10 时只按 Vina 选，命中 {rec["A_只按Vina"]["10"]}/10；\n'
        f'  consensus 混合后是 {rec["C_consensus"]["10"]}/10\n'
        f'· 预算 30 时：只按 Vina {rec["A_只按Vina"]["30"]}/30，\n'
        f'  consensus {rec["C_consensus"]["30"]}/30\n'
        f'· 结论：把对接分数当“过滤掉明显不结合的”，\n'
        f'  别当“排序器”；和活性预测混着用（consensus）才划算'
    )
    ax.text(0.06, 0.66, body, transform=ax.transAxes, fontsize=9.9,
            color=DARK, va="top", zorder=3, linespacing=1.9)

    out = os.path.join(FIG, "c2_vina_activity.png")
    fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("写出", out)


# ------------------------------------------------------------------ 图 3
def fig3():
    d3 = load("d3_box_sensitivity.json")
    per = d3["per_box"]
    names = list(per.keys())

    fig = plt.figure(figsize=(9.2, 7.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.1, 1.0], hspace=0.48, wspace=0.28,
                          left=0.085, right=0.975, top=0.925, bottom=0.085)

    # (a) 三盒 rho 对比
    ax = fig.add_subplot(gs[0, 0])
    rhos = [abs(per[n]["spearman_rho_vina_vs_pIC50"]) for n in names]
    cols = [GREEN, AMBER, RED]
    xs = np.arange(len(names))
    ax.bar(xs, rhos, 0.5, color=cols, zorder=3)
    for x, r, n in zip(xs, rhos, names):
        ax.text(x, r + 0.018, f"|ρ| = {r:.3f}", ha="center", fontsize=10, color=DARK)
    ax.set_xticks(xs)
    ax.set_xticklabels(["A 正确口袋", "B 偏移 8 埃", "C 全蛋白盲搜"], fontsize=9.8)
    ax.set_ylabel("|Spearman ρ|  (Vina 分数 vs 真实活性)", fontsize=9.8)
    ax.set_ylim(0, max(rhos) * 1.28 + 0.05)
    ax.set_title("(a) 盒子瞄偏，排序能力就塌", fontsize=11.3, color=DARK, pad=9)
    ax.grid(axis="y", alpha=0.25, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)

    # (b) pose 落点距离
    ax = fig.add_subplot(gs[0, 1])
    dists = [per[n]["mean_pose_dist_to_true_pocket"] for n in names]
    fracs = [per[n]["frac_pose_within_6A_of_true_pocket"] for n in names]
    ax.bar(xs, dists, 0.5, color=cols, zorder=3)
    for x, d, f, n in zip(xs, dists, fracs, names):
        ax.text(x, d + 1.0, f"{d:.1f} 埃\n落在真口袋 6 埃内：{f*100:.0f}%",
                ha="center", fontsize=9.2, color=DARK, linespacing=1.5)
    ax.set_xticks(xs)
    ax.set_xticklabels(["A 正确口袋", "B 偏移 8 埃", "C 全蛋白盲搜"], fontsize=9.8)
    ax.set_ylabel("最优 pose 质心到真口袋的距离（埃）", fontsize=9.8)
    ax.set_ylim(0, max(dists) * 1.42 + 2)
    ax.set_title("(b) 分子到底被放在哪了", fontsize=11.3, color=DARK, pad=9)
    ax.grid(axis="y", alpha=0.25, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)

    # (c) 跨盒排名崩塌
    ax = fig.add_subplot(gs[1, :])
    cross = d3["cross_box_rank_overlap"]
    for j, (bn, c) in enumerate(cross.items()):
        ranks = c["rank_of_those_in_B"]
        ax.plot(range(1, len(ranks) + 1), ranks, "o-", lw=1.8, ms=7,
                color=cols[j + 1], label=f"盒子 A 的 Top-8 在 {bn} 里的排名", zorder=4)
        for i, r in enumerate(ranks):
            if r is not None:
                ax.text(i + 1, r + 0.6, str(r), ha="center", fontsize=9, color=cols[j + 1])
    ax.axhline(8, color=GREY, ls="--", lw=1.2, zorder=2)
    ax.text(1.05, 8.7, "第 8 名（仍在 Top-8 内）", fontsize=9.2, color=GREY)
    ax.set_xlabel("盒子 A 选出的第几名分子", fontsize=10.2)
    ax.set_ylabel("在另一个盒子里的排名", fontsize=10.2)
    ax.set_ylim(0, max(len(names) + 2, 18))
    ax.set_xticks(range(1, 9))
    ax.legend(fontsize=9.3, loc="upper right", frameon=False)
    ax.grid(alpha=0.25, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("(c) 同一批分子，换个盒子，榜单面目全非", fontsize=11.3, color=DARK, pad=9)

    out = os.path.join(FIG, "c3_box_sensitivity.png")
    fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("写出", out)


# ------------------------------------------------------------------ 图 4
def fig4():
    fig = plt.figure(figsize=(9.2, 6.6))
    gs = fig.add_gridspec(1, 2, wspace=0.30, left=0.075, right=0.975,
                          top=0.90, bottom=0.115)

    # (a) 定位矩阵
    ax = fig.add_subplot(gs[0, 0])
    items = [
        ("AutoDock Vina\n（传统物理力场）", 0.72, 0.62, "必须先定义口袋\n盒子瞄哪，答案就在哪"),
        ("深度学习回归模型\n（如 3D-CNN 打分）", 0.40, 0.34, "要定义口袋\n且绝对值不可靠"),
        ("DiffDock\n（扩散生成式）", 0.08, 0.62, "盲对接\n不需要预先定义口袋"),
    ]
    ax.set_xlim(-0.05, 1.18)
    ax.set_ylim(-0.05, 1.05)
    ax.axis("off")
    ax.add_patch(plt.Rectangle((-0.02, -0.02), 1.22, 1.08, transform=ax.transAxes,
                               fc="#fbfcfd", ec="#dfe6e4", lw=1.2))
    ax.annotate("", xy=(1.05, 0.06), xytext=(0.10, 0.06),
                arrowprops=dict(arrowstyle="->", color=DARK, lw=1.8))
    ax.annotate("", xy=(0.10, 0.90), xytext=(0.10, 0.06),
                arrowprops=dict(arrowstyle="->", color=DARK, lw=1.8))
    ax.text(1.05, -0.015, "需要预先定义口袋 →", ha="right", fontsize=10, color=DARK)
    ax.text(0.055, 0.90, "分数的物理含义清晰度 →", fontsize=10, color=DARK,
            rotation=90, va="top")
    for i, (name, x, y, note) in enumerate(items):
        color = [RED, AMBER, GREEN][i]
        ax.add_patch(plt.Rectangle((x, y), 0.26, 0.155, transform=ax.transAxes,
                                   fc=color, ec="white", lw=1.8, zorder=4))
        ax.text(x + 0.13, y + 0.077, name, transform=ax.transAxes, fontsize=9.4,
                color="white", ha="center", va="center", fontweight="bold",
                zorder=5, linespacing=1.4)
        ax.text(x + 0.13, y - 0.03, note, transform=ax.transAxes, fontsize=8.6,
                color="#666", ha="center", va="top", zorder=5, linespacing=1.5)
    ax.set_title("(a) 三类对接方法的定位", fontsize=11.5, color=DARK, pad=10)

    # (b) 置信度标尺
    ax = fig.add_subplot(gs[0, 1])
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.add_patch(plt.Rectangle((0.03, 0.03), 0.94, 0.90, transform=ax.transAxes,
                               fc="#fbfcfd", ec="#dfe6e4", lw=1.2))
    ax.text(0.07, 0.85, "DiffDock confidence score 怎么读", transform=ax.transAxes,
            fontsize=11.3, color=DARK, fontweight="bold", va="top")
    ax.text(0.07, 0.755,
            "官方 FAQ 原文口径（不是 0~1 的概率！）",
            transform=ax.transAxes, fontsize=9.4, color=RED, va="top")
    rows = [
        ("c > 0", "高置信度", GREEN, 0.60),
        ("-1.5 < c < 0", "中等置信度", AMBER, 0.42),
        ("c < -1.5", "低置信度", RED, 0.24),
    ]
    for code, label, color, y in rows:
        ax.add_patch(plt.Rectangle((0.07, y), 0.20, 0.085, transform=ax.transAxes,
                                   fc=color, ec="none", zorder=3))
        ax.text(0.17, y + 0.042, code, transform=ax.transAxes, fontsize=10.2,
                color="white", ha="center", va="center", fontweight="bold", zorder=4)
        ax.text(0.29, y + 0.042, label, transform=ax.transAxes, fontsize=10,
                color=DARK, va="center", zorder=4)
    ax.text(0.07, 0.175,
            "· 它是模型对自己预测的结构有多确定的打分，\n"
            "  不是结合能，也不适合横向比亲和力\n"
            "· 大配体 / 大蛋白 / 非结合态构象，\n"
            "  这三档阈值要往下调\n"
            "· DiffDock 官方明确：不预测结合能，\n"
            "  想算亲和力请配 GNINA / MM-GBSA / FEP",
            transform=ax.transAxes, fontsize=9.4, color=DARK, va="top",
            linespacing=1.75)

    out = os.path.join(FIG, "c4_positioning.png")
    fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("写出", out)


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "1"):
        fig1()
    if which in ("all", "2"):
        fig2()
    if which in ("all", "3"):
        fig3()
    if which in ("all", "4"):
        fig4()
    print("全部完成")
