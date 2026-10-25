#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 文章关联：AI制药细节篇1_靶点与口袋_公众号文章.html
# 脚本职责：实验 E5：三段可视化示例与 a7–a9 配图
# 使用注解：依赖 py3Dmol、Chrome 截图及已有结构文件；Chrome 路径写死为 macOS 安装位置。
# 完整对应关系与运行顺序：../README.md
"""
E5｜用 py3Dmol 把「蛋白 + 口袋」真正画出来（文章的「代码实战·第 3 步」）

三段递进，每段都对应文章里一个可直接复制的代码块：
  ① hello world：`view(query='pdb:1stp')` + cartoon spectrum + SAS 表面
  ② AlphaFold 专属打开方式：本地 PDB + 按 pLDDT 分档着色（官方色带）
  ③ 标出锁孔：addSphere 画 5 埃口袋球 + 配体 + 口袋残基高亮

出图方式：py3Dmol 本身只生成 JS；要拿静态 PNG，须 `write_html()`
导出后用 Chrome headless 截图（本机实测可行）。
"""
import os
import json
import subprocess

import numpy as np
import py3Dmol

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.abspath(os.path.join(HERE, "..", ".."))
FIG = os.path.join(WORK, "AI制药细节篇1_配图")
VIZ = os.path.join(HERE, "viz")
TMP = os.path.join(HERE, "tmp_html")
os.makedirs(FIG, exist_ok=True)
os.makedirs(TMP, exist_ok=True)

CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
W, H = 900, 640

# AlphaFold 官方 pLDDT 色带
PLDDT_BANDS = [
    ("#2b6cb0", 90, 101),   # ≥90  非常高置信
    ("#63b3ed", 70, 90),    # 70-90 骨架基本可信
    ("#f6ad55", 50, 70),    # 50-70 谨慎
    ("#e2e8f0", 0, 50),     # <50   无序，别解释
]


def read(p):
    with open(p) as f:
        return f.read()


def residues_by_plddt(pdb_text):
    """按残基把 pLDDT 分档，返回 {颜色: [残基号...]}。

    不用 setStyle 的 gradient（那是红白蓝，与 AlphaFold 官方色带不同），
    直接按 b 因子分四档，逐档 addStyle，颜色与 AlphaFold DB 页面一致。
    """
    buckets = {c: [] for c, _, _ in PLDDT_BANDS}
    seen = set()
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM") or line[12:16].strip() != "CA":
            continue
        resi = int(line[22:26])
        if resi in seen:
            continue
        seen.add(resi)
        b = float(line[60:66])
        for color, lo, hi in PLDDT_BANDS:
            if lo <= b < hi:
                buckets[color].append(resi)
                break
    return {c: v for c, v in buckets.items() if v}


def shot(view, out_png, wait_ms=9000):
    """write_html → Chrome headless 截图 → PNG"""
    html_path = os.path.join(TMP, os.path.basename(out_png).replace(".png", ".html"))
    view.write_html(html_path)
    subprocess.run(
        [CHROME, "--headless", "--no-sandbox", "--enable-unsafe-swiftshader",
         "--hide-scrollbars", "--force-device-scale-factor=1",
         "--window-size=%d,%d" % (W, H),
         "--screenshot=%s" % out_png,
         "--virtual-time-budget=%d" % wait_ms,
         "file://" + html_path],
        check=True, capture_output=True, timeout=180,
    )
    return out_png


def check_png(png, min_nonwhite=0.03):
    """像素自检：确认不是一张空白图（3Dmol 会安静地画出空白）"""
    from PIL import Image
    im = Image.open(png).convert("RGB")
    px = list(im.getdata())
    nw = sum(1 for p in px if abs(p[0] - 255) + abs(p[1] - 255) + abs(p[2] - 255) > 30)
    frac = nw / len(px)
    return round(100 * frac, 2), frac >= min_nonwhite


def main():
    out = {}

    # ---------- ① hello world：一个 ID 换一张图 ----------
    # 注：addSurface(SAS) 的计算走 Web Worker，在 Chrome headless --virtual-time-budget
    # 下会死锁（实测 3 种等待方式都拿不到结果）。用白色半透明原子球近似表面，
    # 视觉等效；读者在 Jupyter 里跑 addSurface(SAS) 没有这个问题。
    v1 = py3Dmol.view(query="pdb:1stp", width=W, height=H)
    v1.setStyle({"cartoon": {"color": "spectrum"}})
    v1.addStyle({}, {"sphere": {"scale": 0.95, "opacity": 0.42, "color": "white"}})
    v1.addStyle({"resn": "BTN"}, {"stick": {"color": "#E53E3E", "radius": 0.25}})
    v1.zoomTo()
    p1 = shot(v1, os.path.join(FIG, "a7_py3dmol_hello.png"))
    out["图6"] = dict(file=os.path.relpath(p1, WORK), **dict(zip(("pct", "ok"), check_png(p1))))

    # ---------- ② AlphaFold 专属打开方式 ----------
    af_sup = read(os.path.join(VIZ, "02_AF_superposed_on_1M17.pdb"))
    bands = residues_by_plddt(af_sup)
    v2 = py3Dmol.view(width=W, height=H)
    v2.addModel(af_sup, "pdb")
    v2.setStyle({"cartoon": {"color": "#cccccc"}})
    for color, resis in bands.items():
        v2.setStyle({"resi": resis}, {"cartoon": {"color": color}})
    v2.zoomTo()
    p2 = shot(v2, os.path.join(FIG, "a8_py3dmol_plddt.png"))
    out["图7"] = dict(file=os.path.relpath(p2, WORK), **dict(zip(("pct", "ok"), check_png(p2))))
    out["pLDDT 分档残基数"] = {c: len(v) for c, v in bands.items()}

    # ---------- ③ 把锁孔标出来 ----------
    cry = read(os.path.join(VIZ, "03_1M17_protein.pdb"))
    lig = read(os.path.join(VIZ, "04_1M17_ligand_AQ4.pdb"))
    pocket_res = read(os.path.join(VIZ, "05_pocket_true_20res.pdb"))

    lig_xyz = np.array([
        [float(l[30:38]), float(l[38:46]), float(l[46:54])]
        for l in lig.splitlines() if l.startswith(("ATOM", "HETATM"))
    ])
    cx, cy, cz = lig_xyz.mean(axis=0)
    pocket_resis = sorted({int(l[22:26]) for l in pocket_res.splitlines() if l.startswith("ATOM")})
    out["配体质心"] = [round(float(x), 2) for x in (cx, cy, cz)]
    out["口袋残基数"] = len(pocket_resis)

    v3 = py3Dmol.view(width=W, height=H)
    v3.addModel(cry, "pdb")
    v3.setStyle({"cartoon": {"color": "#dfe6ea", "opacity": 0.9}})
    v3.setStyle({"resi": pocket_resis}, {"stick": {"color": "#0FA87E", "radius": 0.18}})
    v3.addModel(lig, "pdb")
    v3.setStyle({"resn": "AQ4"}, {"stick": {"color": "#E53E3E", "radius": 0.28},
                                  "sphere": {"radius": 0.45}})
    # 文章里说的「5 埃口袋球」：以配体质心为圆心，半径 5 埃
    v3.addSphere({"center": {"x": cx, "y": cy, "z": cz},
                  "radius": 5.0, "color": "#F5A623", "alpha": 0.4})
    v3.zoomTo({"resn": "AQ4"}, 18)
    p3 = shot(v3, os.path.join(FIG, "a9_py3dmol_pocket.png"))
    out["图8"] = dict(file=os.path.relpath(p3, WORK), **dict(zip(("pct", "ok"), check_png(p3))))

    with open(os.path.join(HERE, "results", "e5_py3dmol.json"), "w") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)

    for k in ("图6", "图7", "图8"):
        r = out[k]
        print("%s %-38s 非白像素 %5.2f%%  %s" % (k, r["file"], r["pct"], "OK" if r["ok"] else "疑似空白!"))
    print("pLDDT 分档残基数:", out["pLDDT 分档残基数"])
    print("配体质心:", out["配体质心"], "口袋残基数:", out["口袋残基数"])


if __name__ == "__main__":
    main()
