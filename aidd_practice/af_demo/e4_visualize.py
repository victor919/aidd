#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 文章关联：AI制药细节篇1_靶点与口袋_公众号文章.html
# 脚本职责：实验 E4：结构、PyMOL / ChimeraX 与交互网页
# 使用注解：复用 E1 / E3 算法；输出 viz/ 与上一级 EGFR_AFPocket_3D.html，提供补充可视化。
# 完整对应关系与运行顺序：../README.md
"""
E4｜把 AlphaFold 预测结构 + 检测到的口袋，真正输出成"能看"的东西

产物分三类，覆盖三种使用习惯：
  1) 结构文件（viz/*.pdb）        —— 拖进 PyMOL / ChimeraX / Mol* 就能看
  2) 一键脚本（viz/egfr_pocket.pml / .cxc）—— 命令行党直接跑出图
  3) 交互网页（EGFR_AFPocket_3D.html）    —— 单文件内联坐标，浏览器打开可旋转

内容：
  - AF-P00533-F1 全长（1210 残基），B-factor 列 = pLDDT
  - 1M17 晶体（删水/配体）+ 共晶配体 AQ4（厄洛替尼）
  - AF 全长按高置信区(pLDDT>=70) Kabsch 叠合到晶体坐标系
  - 真实口袋 20 残基（配体任一重原子 5 埃内）
  - 几何法检测到的口袋体素表面点云（三组对照）
"""
import os
import sys
import json
import gzip
import base64

import numpy as np
from scipy.spatial import cKDTree
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.abspath(os.path.join(HERE, "..", ".."))
MVP = os.path.join(WORK, "aidd-mvp")
sys.path.insert(0, os.path.join(MVP, "scripts"))
sys.path.insert(0, os.path.join(MVP, "vendor"))
sys.path.insert(0, HERE)

from _common import raw  # noqa: E402
from e1_af_vs_crystal import parse_pdb, fit_transform, OFFSET, CUTOFF, LIGAND_RESN  # noqa: E402
from e3_pocket_detection import build_grid, detect, SOLVENTS  # noqa: E402

VIZ = os.path.join(HERE, "viz")
FIG = os.path.join(WORK, "AI制药细节篇1_配图")
os.makedirs(VIZ, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

TARGET_POCKET_POINTS = 2000
AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


def seq_from_atoms(orig_lines):
    """从原始 PDB 文本按 CA 出现顺序抽序列（用于核对编号映射）。"""
    seq, seen = [], set()
    for line in orig_lines:
        if not line.startswith("ATOM"):
            continue
        if line[12:16].strip() != "CA":
            continue
        rn = int(line[22:26])
        if rn in seen:
            continue
        seen.add(rn)
        seq.append((rn, line[17:20].strip()))
    return seq


def write_pdb(path, lines):
    with open(path, "w") as f:
        f.write("".join(lines))
    return os.path.getsize(path)


def points_to_pdb(points, resname="PKT", atomname="C", element="C"):
    """把体素点云写成 HETATM 伪原子，可视化时 show spheres 即可。"""
    out = []
    for i, p in enumerate(points):
        out.append(
            "HETATM%5d %-4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n"
            % (i + 1, atomname, resname, (i // 9999) + 1, p[0], p[1], p[2], element)
        )
    out.append("END\n")
    return out


def pocket_surface_points(mask, lo, spacing, target=TARGET_POCKET_POINTS):
    """只取口袋体素的『外壳』，再等间距采样到 target 个点。

    直接画全部体素会有上万个球，浏览器卡死；外壳 + 采样后视觉上是一层半透明壳，
    反而更能看出口袋的形状和位置。
    """
    inner = ndimage.binary_erosion(mask, structure=ndimage.generate_binary_structure(3, 1))
    shell = mask & ~inner
    idx = np.argwhere(shell)
    if len(idx) == 0:
        return np.argwhere(mask).astype(float) * spacing + lo
    if len(idx) > target:
        step = max(1, int(np.ceil(len(idx) / target)))
        idx = idx[::step]
    return idx.astype(float) * spacing + lo


def main():
    pdb_crystal = raw("1M17.pdb")
    pdb_af = os.path.join(HERE, "struct", "egfr_af_v6.pdb")

    # ---------- 1. 解析 ----------
    with open(pdb_crystal) as f:
        cry_all = [l for l in f if l.startswith(("ATOM", "HETATM"))]
    with open(pdb_af) as f:
        af_all = [l for l in f if l.startswith(("ATOM", "HETATM"))]

    c_atoms, c_res = parse_pdb(pdb_crystal)
    a_atoms, a_res = parse_pdb(pdb_af)

    lig = [a for a in c_atoms if a["resn"] == LIGAND_RESN]
    lig_xyz = np.array([a["xyz"] for a in lig])

    lt = cKDTree(lig_xyz)
    pocket_resnums = set()
    for rnum, r in c_res.items():
        if r["hetero"] or r["ca"] is None:
            continue
        for xyz in r["atoms"].values():
            dd, _ = lt.query(np.asarray(xyz).reshape(1, -1), k=1)
            if float(dd[0]) <= CUTOFF:
                pocket_resnums.add(rnum)
                break
    pocket_resnums = sorted(pocket_resnums)

    # ---------- 2. AF → 晶体坐标系（只用高置信区求变换） ----------
    af_by_pdb = {rn + OFFSET: r for rn, r in a_res.items()}
    af_plddt = {a["resnum"]: a["bfac"] for a in a_atoms if a["atomname"] == "CA"}

    common = [n for n in sorted(c_res)
              if n in af_by_pdb and not c_res[n]["hetero"]
              and c_res[n]["ca"] is not None and af_by_pdb[n]["ca"] is not None]
    hi = [i for i, n in enumerate(common) if af_plddt.get(n - OFFSET, 0) >= 70]
    P = np.array([af_by_pdb[n]["ca"] for n in common])
    Q = np.array([c_res[n]["ca"] for n in common])
    R, t = fit_transform(P[hi], Q[hi])          # 高置信区定标
    rmsd_hi = float(np.sqrt((((P[hi] @ R + t) - Q[hi]) ** 2).sum(1).mean()))
    rmsd_all = float(np.sqrt((((P @ R + t) - Q) ** 2).sum(1).mean()))

    # ---------- 3. 输出结构文件 ----------
    sizes = {}

    # 3.1 AF 全长，B 列保留 pLDDT（Mol*/PyMOL 里按 b-factor 着色即可）
    sizes["01_AF_P00533_full_plddt.pdb"] = write_pdb(
        os.path.join(VIZ, "01_AF_P00533_full_plddt.pdb"), af_all + ["END\n"])

    # 3.2 AF 全长 → 叠合到晶体坐标系
    sup = []
    n = 0
    for line in af_all:
        if not line.startswith("ATOM"):
            continue
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        nx = xyz @ R + t
        sup.append(line[:30] + "%8.3f%8.3f%8.3f" % (nx[0], nx[1], nx[2]) + line[54:])
        n += 1
    sizes["02_AF_superposed_on_1M17.pdb"] = write_pdb(
        os.path.join(VIZ, "02_AF_superposed_on_1M17.pdb"), sup + ["END\n"])

    # 3.3 晶体蛋白（删水 + 删配体）与配体本体
    prot = [l for l in cry_all
            if l[17:20].strip() not in SOLVENTS and l[17:20].strip() != LIGAND_RESN]
    ligl = [l for l in cry_all if l[17:20].strip() == LIGAND_RESN]
    sizes["03_1M17_protein.pdb"] = write_pdb(
        os.path.join(VIZ, "03_1M17_protein.pdb"), prot + ["END\n"])
    sizes["04_1M17_ligand_AQ4.pdb"] = write_pdb(
        os.path.join(VIZ, "04_1M17_ligand_AQ4.pdb"), ligl + ["END\n"])

    # 3.4 真实口袋残基（抽出来单独一个文件，方便高亮）
    pk_lines, pk_ca = [], []
    for l in prot:
        rn = int(l[22:26])
        if rn in pocket_resnums:
            pk_lines.append(l)
            if l[12:16].strip() == "CA":
                pk_ca.append(np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]))
    sizes["05_pocket_true_20res.pdb"] = write_pdb(
        os.path.join(VIZ, "05_pocket_true_20res.pdb"), pk_lines + ["END\n"])
    pk_ca = np.array(pk_ca)
    true_center = pk_ca.mean(axis=0)

    # 3.5 几何检测口袋的表面点云（三组）
    af_kinase_resnums = {n - OFFSET for n in a_res if 672 <= (n + OFFSET) <= 995}
    jobs = [
        ("06_pocket_geo_crystal.pdb", pdb_crystal, None),
        ("07_pocket_geo_AF_full.pdb", pdb_af, None),
        ("08_pocket_geo_AF_kinase.pdb", pdb_af, af_kinase_resnums),
    ]
    geo = {}
    for fname, path, filt in jobs:
        atoms, residues = parse_pdb(path)
        heavy = np.array([
            a["xyz"] for a in atoms
            if not a["atomname"].startswith("H")
            and a["resn"] not in SOLVENTS and a["resn"] != LIGAND_RESN
            and (filt is None or a["resnum"] in filt)
        ])
        grid = build_grid(heavy)
        pockets = detect(grid)
        top = pockets[0]
        pts = pocket_surface_points(top["mask"], grid["lo"], grid["spacing"])
        # 与晶体坐标系对齐（只有 AF 组需要）
        if path == pdb_af:
            pts = pts @ R + t
        sizes[fname] = write_pdb(os.path.join(VIZ, fname), points_to_pdb(pts))
        cen = pts.mean(axis=0)
        geo[fname] = dict(
            n_voxel=top["volume_A3"], n_drawn=len(pts),
            dcc_to_true_A=round(float(np.linalg.norm(cen - true_center)), 2),
        )

    # ---------- 4. PyMOL / ChimeraX 一键脚本 ----------
    pml = """# PyMOL 一键出图：load egfr_pocket.pml（在该文件所在目录打开 PyMOL 后 @egfr_pocket.pml）
load 03_1M17_protein.pdb, crystal
load 04_1M17_ligand_AQ4.pdb, ligand
load 02_AF_superposed_on_1M17.pdb, afmodel
load 07_pocket_geo_AF_full.pdb, pocket_af
load 06_pocket_geo_crystal.pdb, pocket_cry

hide everything
# 晶体：灰色卡通；AF：按 B 因子(pLDDT)彩虹着色
show cartoon, crystal;  color grey40, crystal
show cartoon, afmodel;  spectrum b, blue_white_red, afmodel, 50, 100
# 配体：绿色球棍
show sticks, ligand;    show spheres, ligand;  color green, ligand
# 检测到的口袋：半透明球壳
show spheres, pocket_af;   set sphere_scale, 1.0, pocket_af;  color orange, pocket_af;  set transparency, 0.55, pocket_af
show spheres, pocket_cry;  set sphere_scale, 1.0, pocket_cry; color marine,  pocket_cry; set transparency, 0.55, pocket_cry
zoom ligand, 12
ray 1600, 1200; png egfr_pocket.png
"""
    cxc = """# ChimeraX 一键出图：open egfr_pocket.cxc
open 03_1M17_protein.pdb
open 04_1M17_ligand_AQ4.pdb
open 02_AF_superposed_on_1M17.pdb
open 07_pocket_geo_AF_full.pdb
open 06_pocket_geo_crystal.pdb
hide
cartoon
color #6b7280 #1
color bybfactor #2 palette alphafold
show #3
style #3 stick
color #3 lime
style #4 sphere
color #4 orange
transparency #4 55
style #5 sphere
color #5 cornflowerblue
transparency #5 55
view #3
save egfr_pocket.png width 1600 height 1200 supersample 3
"""
    with open(os.path.join(VIZ, "egfr_pocket.pml"), "w") as f:
        f.write(pml)
    with open(os.path.join(VIZ, "egfr_pocket.cxc"), "w") as f:
        f.write(cxc)

    # ---------- 5. 静态多角度 PNG ----------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    plt.rcParams["font.family"] = "Hiragino Sans GB"
    plt.rcParams["axes.unicode_minus"] = False

    # 叠合后的 AF 全长 Cα
    af_ca = []
    af_b = []
    for l in af_all:
        if l.startswith("ATOM") and l[12:16].strip() == "CA":
            xyz = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
            af_ca.append(xyz @ R + t)
            af_b.append(float(l[60:66]))
    af_ca = np.array(af_ca)
    af_b = np.array(af_b)

    kry_ca = np.array([r["ca"] for n, r in sorted(c_res.items())
                       if not r["hetero"] and r["ca"] is not None])

    pkt_af = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                       for l in open(os.path.join(VIZ, "07_pocket_geo_AF_full.pdb"))
                       if l.startswith("HETATM")])
    pkt_cry = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                        for l in open(os.path.join(VIZ, "06_pocket_geo_crystal.pdb"))
                        if l.startswith("HETATM")])

    fig = plt.figure(figsize=(9.6, 10.4))
    views = [
        ("① AF 全长按 pLDDT 着色（灰白=无序区）", 20, 45),
        ("② 叠合到晶体坐标系，聚焦激酶域", 25, 130),
        ("③ 检测到的口袋（橙=AF / 蓝=晶体）+ 厄洛替尼（绿）", 15, 210),
    ]
    for k, (title, elev, azim) in enumerate(views, 1):
        ax = fig.add_subplot(3, 1, k, projection="3d")
        if k == 1:
            sc = ax.scatter(af_ca[:, 0], af_ca[:, 1], af_ca[:, 2],
                            c=af_b, cmap="RdYlBu", vmin=40, vmax=100,
                            s=7, alpha=0.85)
            cb = fig.colorbar(sc, ax=ax, shrink=0.72, pad=0.06)
            cb.set_label("pLDDT", fontsize=9)
            ax.plot(kry_ca[:, 0], kry_ca[:, 1], kry_ca[:, 2],
                    "-", color="#111111", lw=1.6, alpha=0.9, label="1M17 晶体（激酶域）")
            ax.legend(loc="upper left", fontsize=8)
        elif k == 2:
            ax.plot(af_ca[:, 0], af_ca[:, 1], af_ca[:, 2],
                    "-", color="#cccccc", lw=0.7, alpha=0.8, label="AF 全长")
            ax.plot(kry_ca[:, 0], kry_ca[:, 1], kry_ca[:, 2],
                    "-", color="#0FA87E", lw=2.0, label="1M17 晶体（激酶域）")
            ax.scatter(lig_xyz[:, 0], lig_xyz[:, 1], lig_xyz[:, 2],
                       c="#E53E3E", s=55, depthshade=False, label="厄洛替尼 AQ4")
            ax.legend(loc="upper left", fontsize=8)
        else:
            ax.scatter(af_ca[:, 0], af_ca[:, 1], af_ca[:, 2],
                       c="#e5e7eb", s=4, alpha=0.35)
            ax.scatter(pkt_af[:, 0], pkt_af[:, 1], pkt_af[:, 2],
                       c="#F5A623", s=18, alpha=0.30, label="AF 检测口袋")
            ax.scatter(pkt_cry[:, 0], pkt_cry[:, 1], pkt_cry[:, 2],
                       c="#4A7CFE", s=22, alpha=0.65, label="晶体检测口袋")
            ax.scatter(lig_xyz[:, 0], lig_xyz[:, 1], lig_xyz[:, 2],
                       c="#0FA87E", s=70, depthshade=False, label="厄洛替尼 AQ4（真值）")
            ax.legend(loc="upper left", fontsize=8)
        ax.set_title(title, fontsize=10.5, pad=6, color="#1A1A1A")
        ax.view_init(elev=elev, azim=azim)
        ax.set_axis_off()
    fig.suptitle("EGFR（UniProt P00533）· AlphaFold 预测结构与活性口袋",
                 fontsize=13, color="#0F6B52", y=0.975)
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    out_png = os.path.join(FIG, "a6_3d_view.png")
    fig.savefig(out_png, dpi=150)
    plt.close(fig)

    # ---------- 6. 单文件交互 HTML ----------
    def read_text(p):
        with open(p) as f:
            return f.read()

    # 视图①的"全长飘带"同样可以用叠合后的坐标（刚体变换不改变 pLDDT 着色效果），
    # 这样只内联一份全长结构，省掉近 1 MB。
    payload = {
        "af_sup": read_text(os.path.join(VIZ, "02_AF_superposed_on_1M17.pdb")),
        "cry": read_text(os.path.join(VIZ, "03_1M17_protein.pdb")),
        "lig": read_text(os.path.join(VIZ, "04_1M17_ligand_AQ4.pdb")),
        "pkt_af": read_text(os.path.join(VIZ, "07_pocket_geo_AF_full.pdb")),
        "pkt_cry": read_text(os.path.join(VIZ, "06_pocket_geo_crystal.pdb")),
    }
    # PDB 文本重复度极高，gzip 后约为 1/4，再 base64 内联进单文件 HTML
    b64, raw_sz, gz_sz = {}, 0, 0
    for k, v in payload.items():
        b = v.encode()
        raw_sz += len(b)
        c = gzip.compress(b, 9)
        gz_sz += len(c)
        b64[k] = base64.b64encode(c).decode()

    html = HTML_TPL
    for k, v in b64.items():
        html = html.replace("__B64_%s__" % k.upper(), v)

    out_html = os.path.join(WORK, "EGFR_AFPocket_3D.html")
    with open(out_html, "w") as f:
        f.write(html)

    summary = dict(
        rmsd_ca_pLDDT_ge70_A=round(rmsd_hi, 3),
        rmsd_ca_all_312_A=round(rmsd_all, 3),
        n_true_pocket_residues=len(pocket_resnums),
        pocket_resnums_pdb=pocket_resnums,
        geo_pockets=geo,
        files=sizes,
        static_png=os.path.relpath(out_png, WORK),
        interactive_html=os.path.relpath(out_html, WORK),
        html_payload_raw_KB=round(raw_sz / 1024.0, 1),
        html_payload_gzip_KB=round(gz_sz / 1024.0, 1),
    )
    with open(os.path.join(HERE, "results", "e4_visualize.json"), "w") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    print("=== 叠合 ===")
    print("  高置信区(pLDDT>=70) Cα RMSD : %.3f 埃" % rmsd_hi)
    print("  全部 312 残基 Cα RMSD       : %.3f 埃" % rmsd_all)
    print("  真实口袋残基数              : %d" % len(pocket_resnums))
    print("  口袋残基(PDB 编号)          : %s" % pocket_resnums)
    print("=== 几何检测口袋（表面点） ===")
    for k, v in geo.items():
        print("  %-28s 体素%6d  绘点%5d  DCC %6.2f 埃"
              % (k, v["n_voxel"], v["n_drawn"], v["dcc_to_true_A"]))
    print("=== 产物 ===")
    for k, v in sizes.items():
        print("  viz/%-30s %8.1f KB" % (k, v / 1024.0))
    print("  %s" % summary["static_png"])
    print("  %s  (%.1f KB，内联坐标 gzip 后 %.1f KB / 原始 %.1f KB)"
          % (summary["interactive_html"], os.path.getsize(out_html) / 1024.0,
             gz_sz / 1024.0, raw_sz / 1024.0))


HTML_TPL = """<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>EGFR · AlphaFold 结构与活性口袋 3D 交互演示</title>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>
  body{margin:0;background:#F0F4F3;color:#3F3F3F;
       font-family:-apple-system,BlinkMacSystemFont,"PingFang SC","Hiragino Sans GB","Microsoft YaHei",sans-serif;}
  .wrap{max-width:920px;margin:0 auto;padding:22px 16px 40px;}
  h1{font-size:20px;color:#0F6B52;font-weight:500;margin:0 0 6px;}
  .sub{font-size:13px;color:#9A9A9A;margin:0 0 16px;line-height:1.7;}
  .bar{display:flex;flex-wrap:wrap;gap:8px;margin:0 0 12px;}
  button{border:1px solid #d7dde0;background:#fff;color:#0F6B52;font-size:13px;
         padding:7px 13px;border-radius:3px;cursor:pointer;font-family:inherit;}
  button:hover{background:#eaf6f1;}
  button.on{background:#0FA87E;border-color:#0FA87E;color:#fff;}
  #view{width:100%;height:470px;background:#ffffff;border:1px solid #E0E3E6;
        border-radius:4px;position:relative;}
  .hint{font-size:12px;color:#9A9A9A;margin:10px 0 0;line-height:1.7;}
  .legend{margin:14px 0 0;font-size:12.5px;line-height:1.9;color:#5a5a5a;
          background:#F5F7F9;border-radius:3px;padding:12px 16px;}
  .dot{display:inline-block;width:9px;height:9px;border-radius:50%;margin-right:6px;
       vertical-align:middle;}
  .note{font-size:13px;line-height:1.8;margin:14px 0 0;}
  #fallback{display:none;padding:16px;background:#FBF8F1;border-radius:3px;
            font-size:13px;line-height:1.8;color:#7a6a4a;}
</style>
</head>
<body>
<div class="wrap">
  <h1>EGFR · AlphaFold 预测结构与活性口袋（3D 交互演示）</h1>
  <p class="sub">
    UniProt P00533 · AlphaFold DB AF-P00533-F1（v6，1210 残基）｜对照晶体：PDB 1M17（激酶域，共晶配体 AQ4 = 厄洛替尼）<br>
    鼠标左键拖动旋转，滚轮缩放，右键平移。所有坐标均为本机实跑生成，非示意图。
  </p>

  <div class="bar">
    <button id="b1" class="on">① AF 全长 · pLDDT 着色</button>
    <button id="b2">② AF 叠合到晶体</button>
    <button id="b3">③ 口袋检测对照</button>
    <button id="b4">④ 只看口袋 + 配体</button>
  </div>

  <div id="view"></div>
  <p class="hint">提示：pLDDT &lt; 50 的区域（灰白色飘带）是无序区，不应作为对接受体的可靠依据。</p>

  <div class="legend">
    <span class="dot" style="background:#2b6cb0"></span>pLDDT ≥ 90（置信度高）
    &nbsp;&nbsp;<span class="dot" style="background:#63b3ed"></span>70–90（骨架基本可信）
    &nbsp;&nbsp;<span class="dot" style="background:#f6ad55"></span>50–70（谨慎）
    &nbsp;&nbsp;<span class="dot" style="background:#e2e8f0"></span>&lt; 50（无序）<br>
    <span class="dot" style="background:#0FA87E"></span>厄洛替尼 AQ4（共晶配体，口袋真值位置）
    &nbsp;&nbsp;<span class="dot" style="background:#F5A623"></span>AF 结构几何检测口袋
    &nbsp;&nbsp;<span class="dot" style="background:#4A7CFE"></span>晶体结构几何检测口袋
  </div>

  <div id="fallback">
    3D 渲染库（3Dmol.js）没能从网络加载。离线情况下可以直接看静态图
    <code>AI制药细节篇1_配图/a6_3d_view.png</code>，或把 <code>aidd_practice/af_demo/viz/</code>
    里的 .pdb 文件拖进 PyMOL / ChimeraX 查看。
  </div>
</div>

<script>
// 坐标以 gzip + base64 内联，保证单文件、可离线（不依赖任何数据请求）
var GZ = {
  af_sup : "__B64_AF_SUP__",
  cry    : "__B64_CRY__",
  lig    : "__B64_LIG__",
  pkt_af : "__B64_PKT_AF__",
  pkt_cry: "__B64_PKT_CRY__"
};
var PDB = {};

function gunzip(b64){
  var bin = Uint8Array.from(atob(b64), function(c){return c.charCodeAt(0);});
  if (typeof DecompressionStream === "undefined") return Promise.reject("no DecompressionStream");
  var s = new Blob([bin]).stream().pipeThrough(new DecompressionStream("gzip"));
  return new Response(s).text();
}

var viewer = null;
function boot(){
  if (typeof $3Dmol === "undefined") {
    document.getElementById("view").style.display = "none";
    document.getElementById("fallback").style.display = "block";
    return;
  }
  var keys = Object.keys(GZ);
  Promise.all(keys.map(function(k){ return gunzip(GZ[k]); }))
    .then(function(texts){
      keys.forEach(function(k, i){ PDB[k] = texts[i]; });
      viewer = $3Dmol.createViewer("view", {backgroundColor:"white"});
      show(1);
    })
    .catch(function(e){
      document.getElementById("view").style.display = "none";
      var fb = document.getElementById("fallback");
      fb.style.display = "block";
      fb.innerHTML += "<br>（解压失败：" + String(e) + "）";
    });
}

// pLDDT → 颜色
function plddtColor(b){
  if (b >= 90) return "#2b6cb0";
  if (b >= 70) return "#63b3ed";
  if (b >= 50) return "#f6ad55";
  return "#e2e8f0";
}

function clearAll(){
  viewer.removeAllModels(); viewer.removeAllSurfaces(); viewer.removeAllLabels();
  viewer.removeAllShapes();
}

// 按 AlphaFold 官方色带给卡通上色：分四档手动指定 resi 列表
function addPlddtCartoon(pdbText){
  var m = viewer.addModel(pdbText, "pdb");
  var buckets = {"#2b6cb0":[], "#63b3ed":[], "#f6ad55":[], "#e2e8f0":[]};
  var seen = {};
  m.selectedAtoms({}).forEach(function(a){
    if (seen[a.resi]) return;      // 一个残基只取一次
    seen[a.resi] = 1;
    buckets[plddtColor(a.b)].push(a.resi);
  });
  Object.keys(buckets).forEach(function(c){
    if (!buckets[c].length) return;
    m.setStyle({resi: buckets[c]}, {cartoon:{color:c, opacity:0.95}});
  });
  return m;
}

function show(n){
  if(!viewer) return;
  clearAll();
  if(n===1){
    addPlddtCartoon(PDB.af_sup);
    viewer.zoomTo();
    viewer.setViewStyle({style:"outline"});
  } else if(n===2){
    var m1 = addPlddtCartoon(PDB.af_sup);
    var m2 = viewer.addModel(PDB.cry, "pdb");
    m2.setStyle({}, {cartoon:{color:"#0FA87E"}});
    var m3 = viewer.addModel(PDB.lig, "pdb");
    m3.setStyle({}, {stick:{color:"#E53E3E", radius:0.25}, sphere:{radius:0.5}});
    viewer.zoomTo({model:m3});
  } else if(n===3){
    var a = viewer.addModel(PDB.af_sup, "pdb");
    a.setStyle({}, {cartoon:{color:"#dfe6ea", opacity:0.85}});
    var c = viewer.addModel(PDB.cry, "pdb");
    c.setStyle({}, {cartoon:{color:"#9AA7B0", opacity:0.9}});
    var pa = viewer.addModel(PDB.pkt_af, "pdb");
    pa.setStyle({}, {sphere:{color:"#F5A623", radius:0.95, opacity:0.5}});
    var pc = viewer.addModel(PDB.pkt_cry, "pdb");
    pc.setStyle({}, {sphere:{color:"#4A7CFE", radius:0.95, opacity:0.5}});
    var l = viewer.addModel(PDB.lig, "pdb");
    l.setStyle({}, {stick:{color:"#0FA87E", radius:0.3}, sphere:{radius:0.55}});
    viewer.zoomTo({model:l});
  } else if(n===4){
    var pa2 = viewer.addModel(PDB.pkt_af, "pdb");
    pa2.setStyle({}, {sphere:{color:"#F5A623", radius:1.0, opacity:0.45}});
    var pc2 = viewer.addModel(PDB.pkt_cry, "pdb");
    pc2.setStyle({}, {sphere:{color:"#4A7CFE", radius:1.0, opacity:0.45}});
    var l2 = viewer.addModel(PDB.lig, "pdb");
    l2.setStyle({}, {stick:{color:"#0FA87E", radius:0.35}, sphere:{radius:0.7}});
    viewer.zoomTo({model:l2});
  }
  viewer.render();
}

["b1","b2","b3","b4"].forEach(function(id,i){
  document.getElementById(id).onclick = function(){
    ["b1","b2","b3","b4"].forEach(function(j){document.getElementById(j).classList.remove("on");});
    this.classList.add("on");
    show(i+1);
  };
});

window.onload = boot;
</script>
</body>
</html>
"""


if __name__ == "__main__":
    main()
