# -*- coding: utf-8 -*-
"""⑬ 候选卡片图：按四条固定规则选例（重复出现也不替换），本机实跑数据。"""
import os, json
from rdkit import Chem, RDLogger
from rdkit.Chem import Draw
from PIL import Image, ImageDraw, ImageFont
RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results")
OUT = "/Users/dujiepeng/WorkBuddy/Work/AI制药第13篇_配图/图12_候选卡片.png"
e26 = json.load(open(os.path.join(RES, "e26_disagreement.json")))
cards = e26["generated"]["cards"]

W, H = 1200, 900
img = Image.new("RGB", (W, H), "white")
d = ImageDraw.Draw(img)
F = "/System/Library/Fonts/Hiragino Sans GB.ttc"
def font(sz, bold=False):
    try: return ImageFont.truetype(F, sz, index=1 if bold else 0)
    except Exception: return ImageFont.load_default()
fb, fn, fs = font(25, True), font(19), font(16)

d.text((30, 16), "图 12｜按固定规则选例：重复出现也不事后替换", font=fb, fill="#1A1A1A")
d.text((30, 52), "四个数字都来自本机实跑；高分/分歧/相似性都不是实验活性。", font=fs, fill="#9A9A9A")

for i, c in enumerate(cards):
    col, row = i % 2, i // 2
    x0, y0 = 30 + col * 590, 92 + row * 390
    d.rounded_rectangle((x0, y0, x0 + 560, y0 + 365), radius=10, outline="#D5DBDF", width=2, fill="#FAFBFC")
    d.text((x0 + 16, y0 + 12), c["rule"], font=font(20, True), fill="#0F6B52")
    d.text((x0 + 16, y0 + 42),
           f"RF {c['RF']:.2f} / GNN {c['GNN']:.2f} / NN {c['nn']:.2f}", font=fs, fill="#3F3F3F")
    d.text((x0 + 16, y0 + 64),
           f"MW {c['MW']:.0f} / SA {c['SA']:.2f} / QED {c['QED']:.2f} / 拟合池{'已有' if c['in_fit'] else '未见'}",
           font=fs, fill="#C0392B" if c["in_fit"] else "#3F3F3F")
    m = Chem.MolFromSmiles(c["smiles"])
    if m:
        try:
            mi = Draw.MolToImage(m, size=(480, 230), kekulize=True)
            img.paste(mi, (x0 + 40, y0 + 98))
        except Exception as e:
            d.text((x0 + 40, y0 + 150), f"结构绘制失败：{e}", font=fs, fill="#9A9A9A")

img.save(OUT, quality=95)
print("wrote", OUT)
