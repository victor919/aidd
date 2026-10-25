# -*- coding: utf-8 -*-
# 文章关联：AI制药实践篇_公众号文章.html
# 脚本职责：实验：第四站及 KRAS 案例，共价弹头 SMARTS 过滤
# 使用注解：SMARTS 命中表示结构片段匹配，不能据此确认共价结合或真实活性。
# 完整对应关系与运行顺序：README.md
"""共价弹头 SMARTS 过滤 —— KRAS 第四站实操（本机实测）
分子结构全部来自 PubChem 官方记录（2026-10-04）
"""
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors
RDLogger.DisableLog('rdApp.*')

MOLS = [
    # 名称, SMILES, PubChem CID, 身份
    ("Sotorasib",  "C[C@H]1CN(CCN1C2=NC(=O)N(C3=NC(=C(C=C32)F)C4=C(C=CC=C4F)O)C5=C(C=CN=C5C(C)C)C)C(=O)C=C", 137278711, "KRAS G12C 共价抑制剂 / 上市"),
    ("Adagrasib",  "CN1CCC[C@H]1COC2=NC3=C(CCN(C3)C4=CC=CC5=C4C(=CC=C5)Cl)C(=N2)N6CCN([C@H](C6)CC#N)C(=O)C(=C)F", 138611145, "KRAS G12C 共价抑制剂 / 上市"),
    ("Divarasib",  "C[C@H]1CN(CCN1C2=NC(=NC3=C(C(=C(C=C32)Cl)C4=C(C(=CC(=N4)N)C)C(F)(F)F)F)OC[C@@H]5CCCN5C)C(=O)C=C", 146624881, "KRAS G12C 共价抑制剂 / III 期"),
    ("Osimertinib","CN1C=C(C2=CC=CC=C21)C3=NC(=NC=C3)NC4=C(C=C(C(=C4)NC(=O)C=C)N(C)CCN(C)C)OC", 71496458, "EGFR 三代药 / 丙烯酰胺弹头"),
    ("Nirmatrelvir","CC1([C@@H]2[C@H]1[C@H](N(C2)C(=O)[C@H](C(C)(C)C)NC(=O)C(F)(F)F)C(=O)N[C@@H](C[C@@H]3CCNC3=O)C#N)C", 155903259, "3CLpro 抑制剂 / 腈基弹头"),
    ("Ibrutinib",  "C=CC(=O)N1CCC[C@H](C1)N2C3=C(C(=N2)C4=CC=C(C=C4)OC5=CC=CC=C5)C(=NC=N3)N", 24821094, "BTK 抑制剂 / 丙烯酰胺弹头"),
    ("Imatinib",   "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5", 5291, "BCR-ABL / 非共价对照"),
    ("Ritonavir",  "CC(C)C1=NC(=CS1)CN(C)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CC2=CC=CC=C2)C[C@@H]([C@H](CC3=CC=CC=C3)NC(=O)OCC4=CN=CS4)O", 392622, "CYP3A4 抑制剂 / 非共价对照"),
]

# 三类常见亲电弹头 SMARTS
WARHEADS = [
    ("丙烯酰胺 acrylamide", "C=CC(=O)N"),          # KRAS / EGFR / BTK 三联
    ("氯乙酰胺 chloroacetamide", "ClCC(=O)N"),      # 活性更高、选择性更差
    ("腈基 nitrile", "C#N"),                        # 可逆共价（奈玛特韦用）
]
pats = [(n, Chem.MolFromSmarts(s)) for n, s in WARHEADS]

print("=" * 116)
print("表 1｜共价弹头 SMARTS 扫描：谁带「抓手」，抓的是什么（本机实测，结构取自 PubChem）")
print("=" * 116)
print("%-14s %-30s %10s %8s %8s   %s" % ("分子", "身份", "丙烯酰胺", "氯乙酰胺", "腈基", "重原子数"))
rows = []
for name, smi, cid, ident in MOLS:
    m = Chem.MolFromSmiles(smi)
    assert m is not None, name
    hits = {}
    for wn, p in pats:
        hits[wn] = m.HasSubstructMatch(p)
    hac = m.GetNumHeavyAtoms()
    rows.append((name, ident, hits, hac, m, cid))
    print("%-14s %-30s %10s %8s %8s   %d" % (
        name, ident[:28], "● 命中" if hits["丙烯酰胺 acrylamide"] else "—",
        "● 命中" if hits["氯乙酰胺 chloroacetamide"] else "—",
        "● 命中" if hits["腈基 nitrile"] else "—", hac))

print()
print("=" * 116)
print("表 2｜Rule of Three（片段筛选标准）vs 真实上市药：起点和终点根本不是同一个尺寸")
print("=" * 116)
print("%-14s %8s %8s %6s %6s %8s   %s" % ("分子", "MW", "cLogP", "HBD", "HBA", "Ro3 通过", "说明"))
ro3_pass = 0
for name, ident, hits, hac, m, cid in rows:
    mw = Descriptors.MolWt(m)
    logp = Descriptors.MolLogP(m)
    hbd = rdMolDescriptors.CalcNumHBD(m)
    hba = rdMolDescriptors.CalcNumHBA(m)
    ok = mw <= 300 and logp <= 3 and hbd <= 3 and hba <= 3
    ro3_pass += ok
    note = "符合片段起点" if ok else "已超出片段区间（优化后的成品）"
    print("%-14s %8.1f %8.2f %6d %6d %8s   %s" % (name, mw, logp, hbd, hba, "通过" if ok else "不通过", note))
print()
print("→ 8 个分子里只有 %d 个满足 Rule of Three。片段筛选找的是「起点」，不是「成品」。" % ro3_pass)

print()
print("=" * 116)
print("表 3｜含丙烯酰胺的分子的成药性：带弹头的代价（本机实测）")
print("=" * 116)
print("%-14s %8s %8s %8s %8s %8s" % ("分子", "MW", "cLogP", "TPSA", "RotB", "Fsp3"))
for name, ident, hits, hac, m, cid in rows:
    if not hits["丙烯酰胺 acrylamide"]:
        continue
    print("%-14s %8.1f %8.2f %8.1f %8d %8.2f" % (
        name, Descriptors.MolWt(m), Descriptors.MolLogP(m),
        Descriptors.TPSA(m), rdMolDescriptors.CalcNumRotatableBonds(m),
        rdMolDescriptors.CalcFractionCSP3(m)))
