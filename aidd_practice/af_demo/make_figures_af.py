#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 文章关联：AI制药细节篇1_靶点与口袋_公众号文章.html
# 脚本职责：正文配图：a1–a5
# 使用注解：读取 E1 / E3 / E7 结果及结构；先完成实验再绘图，a7–a9 另由 E5 生成。
# 完整对应关系与运行顺序：../README.md
"""
⑨ 细节篇第 1 课配图（全部基于本机实跑结果）
图1 序列→三维折叠链路（概念图）
图2 AF 全长模型的 pLDDT 沿序列分布（真实 B-factor 列）+ 激酶域高亮
图3 AF vs 晶体：偏差分层 + 迭代剔除收敛曲线
图4 口袋残基的 pLDDT 排名分布（证伪"pLDDT 找口袋"）
图5 几何口袋检测：三组结构的 top 口袋体积 / 召回 / 精确率 / DCC
"""
import os
import sys
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.abspath(os.path.join(HERE, "..", ".."))
MVP = os.path.join(WORK, "aidd-mvp")
sys.path.insert(0, os.path.join(MVP, "scripts"))
sys.path.insert(0, os.path.join(MVP, "vendor"))

OUT = os.path.join(WORK, "AI制药细节篇1_配图")
os.makedirs(OUT, exist_ok=True)

plt.rcParams["font.sans-serif"] = ["Hiragino Sans GB", "Arial Unicode MS", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

GREEN = "#0F6B52"
BLUE = "#2B5F8F"
GREY = "#6B7280"
AMBER = "#C77B27"
RED = "#B03A2E"


def box(ax, x, y, w, h, text, fc, tc="white", fs=9.5, bold=False):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                boxstyle="round,pad=0.015,rounding_size=0.02",
                                fc=fc, ec="none", lw=0))
    ax.text(x, y, text, ha="center", va="center", color=tc,
            fontsize=fs, fontweight="700" if bold else "400", zorder=5)


def arrow(ax, x1, y1, x2, y2, c=GREY):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                 mutation_scale=11, lw=1.5, color=c,
                                 shrinkA=1, shrinkB=1, zorder=4))


# ============================================================ 图1
def fig1():
    fig, ax = plt.subplots(figsize=(9.2, 3.5), dpi=170)
    ax.set_xlim(0, 100); ax.set_ylim(0, 34)
    ax.axis("off")

    box(ax, 10, 24, 17, 8, "FASTA 序列\n（一维字母）", "#8FA8B8", fs=9)
    box(ax, 33, 24, 18, 8, "多序列比对 MSA\n（演化耦合信息）", "#7FA6C4", fs=9)
    box(ax, 57, 24, 19, 8, "Evoformer / Transformer\n（配对+单序列表征）", BLUE, fs=9)
    box(ax, 82, 24, 19, 8, "结构模块\n输出三维坐标", GREEN, fs=9)

    arrow(ax, 19, 24, 23, 24); arrow(ax, 42.5, 24, 47, 24); arrow(ax, 67, 24, 72, 24)

    box(ax, 82, 11.5, 24, 7.5, "每个残基附一个 pLDDT\n（0-100 局部置信度）", AMBER, fs=9)
    arrow(ax, 82, 19.8, 82, 15.5, AMBER)

    ax.text(50, 3.0, "图 1｜从一维序列到三维坐标：序列是输入，MSA 提供约束，结构模块给出坐标与置信度",
            ha="center", color=GREY, fontsize=8.5)
    fig.tight_layout()
    p = os.path.join(OUT, "a1_folding.png")
    fig.savefig(p, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return p


# ============================================================ 图2
def fig2():
    from e1_af_vs_crystal import parse_pdb, OFFSET
    a_atoms, a_res = parse_pdb(os.path.join(HERE, "struct", "egfr_af_v6.pdb"))
    ks = sorted(a_res)
    plddt = np.array([a_res[k]["ca"][0] * 0 + next(
        (a["bfac"] for a in a_atoms if a["resnum"] == k and a["atomname"] == "CA"), np.nan)
        for k in ks])
    x = np.array(ks)
    fig, ax = plt.subplots(figsize=(9.2, 3.6), dpi=170)
    ax.plot(x, plddt, lw=1.0, color=GREEN, zorder=3)
    ax.axhspan(0, 50, color=RED, alpha=0.06, zorder=1)
    ax.axhspan(50, 70, color=AMBER, alpha=0.07, zorder=1)
    ax.axhspan(70, 90, color=BLUE, alpha=0.06, zorder=1)
    ax.axhspan(90, 100, color=GREEN, alpha=0.09, zorder=1)
    ax.axvspan(696, 1019, color=AMBER, alpha=0.14, zorder=2)
    ax.text(857, 96, "激酶域 UniProt 696-1019", ha="center", color=AMBER,
            fontsize=9, fontweight="700")
    ax.text(8, 30, "低置信（不宜解释坐标）", color=RED, fontsize=8)
    ax.text(8, 58, "一般", color=GREY, fontsize=8)
    ax.text(8, 78, "可信", color=BLUE, fontsize=8)
    ax.text(8, 93, "很高", color=GREEN, fontsize=8)
    ax.set_xlabel("残基编号（UniProt P00533）", fontsize=9, color=GREY)
    ax.set_ylabel("pLDDT", fontsize=9, color=GREY)
    ax.set_ylim(0, 100); ax.set_xlim(0, 1215)
    ax.tick_params(labelsize=8, colors=GREY)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#D9DEE3")
    ax.set_title("AlphaFold 对 EGFR 全长序列的逐残基置信度（本机实跑，AlphaFold DB v6）",
               fontsize=10, color="#1A1A1A", pad=8)
    fig.tight_layout()
    p = os.path.join(OUT, "a2_plddt_profile.png")
    fig.savefig(p, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return p


# ============================================================ 图3
def fig3():
    e1 = json.load(open(os.path.join(HERE, "results", "e1_af_vs_crystal.json")))
    e7 = json.load(open(os.path.join(HERE, "results", "e7_rigor_check.json")))
    g = e1["e1_alpha_vs_crystal"]
    q1 = e7["Q1_rmsd_口径"]

    # 左图：把「残基集合」与「叠合方式」两个变量拆开，避免"裁剪=精度提升"的误读
    labels = ["A 全部 312 个\n全部拟合",
              "B 高置信 239 个\n沿用 A 变换",
              "C 高置信 239 个\n子集独立拟合",
              "D 全部 312 个\n沿用 C 变换"]
    vals = [q1["A_全部残基+全部拟合"], q1["B_高置信子集+沿用全局变换"],
            q1["C_高置信子集+子集独立拟合"], q1["D_全部残基+沿用高置信拟合"]]
    colors = [RED, AMBER, GREEN, BLUE]

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.5), dpi=170,
                             gridspec_kw={"width_ratios": [1.35, 1]})
    ax = axes[0]
    bars = ax.bar(range(4), vals, color=colors, width=0.62, zorder=3)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.32, f"{v:.2f}",
                ha="center", fontsize=8.5, fontweight="700", color="#1A1A1A")
    ax.set_xticks(range(4)); ax.set_xticklabels(labels, fontsize=7.6, linespacing=1.35)
    ax.set_ylabel("Cα RMSD (埃)", fontsize=9, color=GREY)
    ax.set_ylim(0, max(vals) * 1.30)
    ax.set_title("同一模型，换残基集合 / 换叠合方式，RMSD 差 15 倍",
                 fontsize=9.6, color="#1A1A1A", pad=8)
    ax.tick_params(labelsize=8, colors=GREY)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#D9DEE3")
    ax.grid(axis="y", ls=":", color="#D9DEE3", lw=0.7, zorder=0)
    # 标注两个"差值"分别来自哪个变量
    ax.annotate("", xy=(1, q1["B_高置信子集+沿用全局变换"] + 0.4),
                xytext=(0, q1["A_全部残基+全部拟合"] - 0.4),
                arrowprops=dict(arrowstyle="-", color=AMBER, lw=0.9, ls="--"))
    ax.text(0.5, 7.6, "只换残基集合\n-6.82", ha="center", fontsize=7.4, color=AMBER)
    ax.annotate("", xy=(2, q1["C_高置信子集+子集独立拟合"] + 0.4),
                xytext=(1, q1["B_高置信子集+沿用全局变换"] - 0.4),
                arrowprops=dict(arrowstyle="-", color=GREEN, lw=0.9, ls="--"))
    ax.text(1.5, 3.2, "再换叠合方式\n-4.00", ha="center", fontsize=7.4, color=GREEN)

    ax = axes[1]
    ref = g["rmsd_iterative_trim"]
    nres = [s["n_residues"] for s in e7["Q2_迭代剔除"]]
    ax.plot(range(1, len(ref) + 1), ref, marker="o", ms=5,
            color=BLUE, lw=1.8, zorder=3)
    for i, (v, n) in enumerate(zip(ref, nres)):
        ax.annotate(f"{v:.2f}", (i + 1, v), textcoords="offset points",
                    xytext=(0, 8), ha="center", fontsize=8, color=BLUE)
    ax.set_xticks(range(1, len(ref) + 1))
    ax.set_xticklabels([f"{i+1}\n({n})" for i, n in enumerate(nres)], fontsize=7.6)
    ax.set_xlabel("迭代轮次（下方括号：该轮剩余残基数）", fontsize=8.5, color=GREY)
    ax.set_ylabel("Cα RMSD (埃)", fontsize=9, color=GREY)
    ax.set_ylim(0, max(ref) * 1.35)
    ax.set_title("每轮剔除偏差最大的 20% 后重新拟合", fontsize=9.6, color="#1A1A1A", pad=8)
    ax.tick_params(labelsize=8, colors=GREY)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#D9DEE3")
    ax.grid(axis="y", ls=":", color="#D9DEE3", lw=0.7, zorder=0)

    fig.tight_layout()
    p = os.path.join(OUT, "a3_rmsd_layers.png")
    fig.savefig(p, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return p


# ============================================================ 图4
def fig4():
    e1 = json.load(open(os.path.join(HERE, "results", "e1_af_vs_crystal.json")))
    e2 = e1["e2_plddt_as_pocket_detector"]
    ranks = np.array(e2["pocket_residue_pLDDT_ranks"])
    n = e2["n_residues_alpha_model"]
    base = e2["random_baseline_precision"] * 100

    fig, ax = plt.subplots(figsize=(9.2, 3.2), dpi=170)
    y = np.ones(len(ranks))
    ax.scatter(ranks, y, s=120, color=AMBER, zorder=4, edgecolor="white", lw=1)
    for i, r in enumerate(sorted(ranks)):
        ax.annotate(str(r), (r, 1), textcoords="offset points", xytext=(0, -18),
                    ha="center", fontsize=7.5, color=AMBER)
    ax.set_xlim(0, n)
    ax.set_ylim(0.4, 1.6)
    ax.set_yticks([])
    ax.axvspan(0, 200, color=RED, alpha=0.10, zorder=1)
    ax.text(100, 1.5, "pLDDT 前 200 名\n（口袋残基数 = 0）", ha="center",
            color=RED, fontsize=8.5, fontweight="700")
    ax.axvline(n * base / 100, ls="--", lw=1.2, color=GREY)
    ax.text(n * base / 100 + 12, 1.35, f"随机基线 %d%%" % round(base),
            fontsize=8, color=GREY)
    ax.set_xlabel("残基在全模型 pLDDT 降序排名中的位置（共 1210 个残基）",
            fontsize=8.5, color=GREY)
    ax.set_title("20 个真实口袋残基在 pLDDT 排序中的位置：全部落在中后段",
                 fontsize=10, color="#1A1A1A", pad=8)
    ax.tick_params(labelsize=8, colors=GREY)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color("#D9DEE3")
    fig.tight_layout()
    p = os.path.join(OUT, "a4_plddt_rank.png")
    fig.savefig(p, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return p


# ============================================================ 图5
def fig5():
    e3 = json.load(open(os.path.join(HERE, "results", "e3_pocket_detection.json")))
    keys = ["A_crystal_1M17", "B_alphafold_full", "C_alphafold_kinase_domain"]
    names = ["1M17\n晶体", "AF 全长", "AF 裁到\n激酶域"]
    vols, recalls, precs, dccs = [], [], [], []
    for k in keys:
        r = e3[k]
        vols.append(r["top1_volume_A3"])
        row = r["ligand_overlap"]["top5"][0] if r.get("ligand_overlap") else None
        recalls.append(row["recall"] * 100 if row else 0.0)
        precs.append(row["precision"] * 100 if row else 0.0)
        dccs.append(r["top1_dcc_A"])

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.3), dpi=170)
    ax = axes[0]
    bars = ax.bar([0, 1, 2], vols, color=[GREEN, RED, BLUE], width=0.62, zorder=3)
    for b, v in zip(bars, vols):
        ax.text(b.get_x() + b.get_width() / 2, v + 500, f"{int(v)}",
                ha="center", fontsize=8.5, fontweight="700", color="#1A1A1A")
    ax.set_xticks([0, 1, 2]); ax.set_xticklabels(names, fontsize=8.5)
    ax.set_ylabel("最大口袋体积 (立方埃)", fontsize=8.5, color=GREY)
    ax.set_title("几何法「找到」的口袋有多大", fontsize=10, color="#1A1A1A", pad=8)
    ax.set_ylim(0, max(vols) * 1.22)
    ax.tick_params(labelsize=8, colors=GREY)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#D9DEE3")
    ax.grid(axis="y", ls=":", color="#D9DEE3", lw=0.7, zorder=0)

    ax = axes[1]
    x = np.arange(3)
    w = 0.55
    bars = ax.bar(x, dccs, w, color=[GREEN, RED, BLUE], zorder=3)
    for b, v in zip(bars, dccs):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.25, f"{v:.2f}",
                ha="center", fontsize=8.5, fontweight="700", color="#1A1A1A")
    ax.set_xticks(x); ax.set_xticklabels(names, fontsize=8.5)
    ax.set_ylabel("DCC: 预测口袋质心到真实位点的距离 (埃)", fontsize=8, color=GREY)
    ax.set_title("预测口袋离真实位点还有多远", fontsize=10, color="#1A1A1A", pad=8)
    ax.set_ylim(0, max(dccs) * 1.3)
    ax.legend(fontsize=7.5, frameon=False, loc="upper right")
    ax.tick_params(labelsize=8, colors=GREY)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#D9DEE3")
    ax.grid(axis="y", ls=":", color="#D9DEE3", lw=0.7, zorder=0)

    fig.tight_layout()
    p = os.path.join(OUT, "a5_pocket_detection.png")
    fig.savefig(p, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return p


if __name__ == "__main__":
    paths = []
    for f in (fig1, fig2, fig3, fig4, fig5):
        try:
            paths.append(f())
            print("OK", paths[-1])
        except Exception as e:
            import traceback
            print("FAIL", f.__name__, type(e).__name__, e)
            traceback.print_exc()
    print("共生成", len(paths), "张图")
