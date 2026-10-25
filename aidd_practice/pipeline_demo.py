# -*- coding: utf-8 -*-
# 文章关联：AI制药实践篇_公众号文章.html
# 脚本职责：实验：苗头到先导的分子性质、相似性与过滤演示
# 使用注解：结果输出到终端；真实分子案例用于解释流水线，不构成药效验证。
# 完整对应关系与运行顺序：README.md
"""AI 制药实践篇 · 本机实测脚本
数据源：PubChem PUG REST（2026-10-04 拉取，SMILES 与分子式均为官方记录）
计算：RDKit 2026.03.5
"""
import math, json
from rdkit import Chem, RDLogger, DataStructs
from rdkit.Chem import Descriptors, Crippen, QED, Lipinski, rdMolDescriptors
from rdkit.Chem import rdFingerprintGenerator, rdFMCS
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
RDLogger.DisableLog('rdApp.*')

# 真实分子：SMILES 与分子式来自 PubChem（CID 见注释）
MOLS = [
    # name, smiles, CID, 靶点/身份, Ki or IC50 (nM), 活性类型
    ("PF-00835231", "CC(C)C[C@@H](C(=O)N[C@@H](C[C@@H]1CCNC1=O)C(=O)CO)NC(=O)C2=CC3=C(N2)C=CC=C3OC", 11561899, "SARS-CoV-2 3CLpro 抑制剂（苗头/早期先导）", 0.27, "Ki"),
    ("Nirmatrelvir", "CC1([C@@H]2[C@H]1[C@H](N(C2)C(=O)[C@H](C(C)(C)C)NC(=O)C(F)(F)F)C(=O)N[C@@H](C[C@@H]3CCNC3=O)C#N)C", 155903259, "SARS-CoV-2 Mpro 抑制剂（优化后临床候选）", 3.11, "Ki"),
    ("Boceprevir", "CC1([C@@H]2[C@H]1[C@H](N(C2)C(=O)[C@H](C(C)(C)C)NC(=O)NC(C)(C)C)C(=O)NC(CC3CCC3)C(=O)C(=O)N)C", 10324367, "HCV NS3/4A 抑制剂（双环片段灵感来源）", None, None),
    ("Ensitrelvir", "CN1C=C2C=C(C(=CC2=N1)Cl)NC3=NC(=O)N(C(=O)N3CC4=CC(=C(C=C4F)F)F)CC5=NN(C=N5)C", 162533924, "SARS-CoV-2 3CLpro 非共价抑制剂", 13.0, "IC50"),
    ("Imatinib", "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5", 5291, "BCR-ABL 酪氨酸激酶抑制剂", None, None),
    ("Ritonavir", "CC(C)C1=NC(=CS1)CN(C)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CC2=CC=CC=C2)C[C@@H]([C@H](CC3=CC=CC=C3)NC(=O)OCC4=CN=CS4)O", 392622, "CYP3A4 抑制剂（Paxlovid 中的增效剂）", None, None),
]

def props(mol):
    return dict(
        MW=Descriptors.MolWt(mol),
        cLogP=Crippen.MolLogP(mol),
        TPSA=Descriptors.TPSA(mol),
        HBD=rdMolDescriptors.CalcNumHBD(mol),
        HBA=rdMolDescriptors.CalcNumHBA(mol),
        RotB=rdMolDescriptors.CalcNumRotatableBonds(mol),
        Fsp3=rdMolDescriptors.CalcFractionCSP3(mol),
        ArRing=rdMolDescriptors.CalcNumAromaticRings(mol),
        HAC=mol.GetNumHeavyAtoms(),
        QED=QED.qed(mol),
        Formula=rdMolDescriptors.CalcMolFormula(mol),
    )

rows = []
mols = {}
print("=" * 118)
print("表 1｜六个真实分子的成药性画像（RDKit 本机实测）")
print("=" * 118)
hdr = ("分子", "分子式", "MW", "cLogP", "TPSA", "HBD", "HBA", "RotB", "Fsp3", "HAC", "QED", "Lipinski违")
print("%-14s %-13s %6s %6s %6s %4s %4s %5s %6s %4s %6s %4s" % hdr)
for name, smi, cid, ident, act, act_type in MOLS:
    m = Chem.MolFromSmiles(smi)
    assert m is not None, name
    mols[name] = m
    p = props(m)
    # Lipinski 五规则违规数（MW<=500, logP<=5, HBD<=5, HBA<=10）
    viol = sum([p["MW"] > 500, p["cLogP"] > 5, p["HBD"] > 5, p["HBA"] > 10])
    p["Lipinski"] = viol
    # Veber: RotB<=10, TPSA<=140
    p["Veber"] = (p["RotB"] <= 10 and p["TPSA"] <= 140)
    p["name"], p["ident"], p["act"], p["act_type"] = name, ident, act, act_type
    rows.append(p)
    print("%-14s %-13s %6.1f %6.2f %6.1f %4d %4d %5d %6.2f %4d %6.3f %4d" % (
        name, p["Formula"], p["MW"], p["cLogP"], p["TPSA"], p["HBD"], p["HBA"], p["RotB"], p["Fsp3"], p["HAC"], p["QED"], viol))

print()
print("=" * 100)
print("表 2｜命中优化这一步到底优化了什么：PF-00835231 → Nirmatrelvir")
print("=" * 100)
a = next(r for r in rows if r["name"] == "PF-00835231")
b = next(r for r in rows if r["name"] == "Nirmatrelvir")
print("%-16s %12s %12s %10s" % ("指标", "PF-00835231", "Nirmatrelvir", "变化"))
for key, lab, fmt in [("Formula","分子式","%s"),("MW","分子量","%.1f"),("cLogP","cLogP","%.2f"),
                      ("TPSA","TPSA","%.1f"),("HBD","氢键供体","%d"),("HBA","氢键受体","%d"),
                      ("RotB","可旋转键","%d"),("Fsp3","Fsp3","%.2f"),("HAC","重原子数","%d"),("QED","QED","%.3f")]:
    va, vb = a[key], b[key]
    if isinstance(va, (int, float)):
        delta = "%+.2f" % (vb - va) if fmt != "%d" else "%+d" % (vb - va)
        print("%-16s %12s %12s %10s" % (lab, fmt % va, fmt % vb, delta))
    else:
        print("%-16s %12s %12s %10s" % (lab, va, vb, "-"))

print()
print("=" * 100)
print("表 3｜配体效率：活性高不等于分子好（pAct 由 Ki/IC50 换算，nM）")
print("=" * 100)
print("%-16s %8s %10s %6s %8s %8s" % ("分子", "Ki/IC50", "pAct", "HAC", "LE", "LipE"))
for r in rows:
    if r["act"] is None:
        continue
    pAct = 9 - math.log10(r["act"])          # nM -> pAct
    le = 1.37 * pAct / r["HAC"]
    lipE = pAct - r["cLogP"]
    print("%-16s %8.2f %10.2f %6d %8.3f %8.2f" % (r["name"], r["act"], pAct, r["HAC"], le, lipE))

print()
print("=" * 100)
print("表 4｜结构相似度（Morgan 指纹 r=2 / 2048 位，Tanimoto）")
print("=" * 100)
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
fps = {n: gen.GetFingerprint(m) for n, m in mols.items()}
order = ["PF-00835231", "Nirmatrelvir", "Boceprevir", "Ensitrelvir", "Imatinib", "Ritonavir"]
print("%-16s" % "" + "".join("%14s" % n[:12] for n in order))
for n1 in order:
    line = "%-16s" % n1
    for n2 in order:
        line += "%14.3f" % DataStructs.TanimotoSimilarity(fps[n1], fps[n2])
    print(line)

print()
mcs = rdFMCS.FindMCS([mols["Nirmatrelvir"], mols["Boceprevir"]], completeRingsOnly=True, ringMatchesRingOnly=True)
print("Nirmatrelvir 与 Boceprevir 的最大公共子结构：%d 个原子 / %d 个键，SMARTS = %s"
      % (mcs.numAtoms, mcs.numBonds, mcs.smartsString))
mcs2 = rdFMCS.FindMCS([mols["Nirmatrelvir"], mols["PF-00835231"]], completeRingsOnly=True, ringMatchesRingOnly=True)
print("Nirmatrelvir 与 PF-00835231 的最大公共子结构：%d 个原子 / %d 个键，SMARTS = %s"
      % (mcs2.numAtoms, mcs2.numBonds, mcs2.smartsString))

print()
print("=" * 100)
print("表 5｜结构警报过滤（RDKit FilterCatalog：PAINS + BRENK）——真药也会被『误杀』")
print("=" * 100)
params = FilterCatalogParams()
params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
catalog = FilterCatalog(params)
for n, m in mols.items():
    entry = catalog.GetFirstMatch(m)
    hits = catalog.GetMatches(m)
    names = [e.GetDescription() for e in hits]
    print("%-16s %-8s %s" % (n, "命中" if names else "通过", "；".join(sorted(set(names))[:3]) or "—"))

print()
print("=" * 100)
print("表 6｜苗头化合物的三道体检（Lipinski / Veber / QED / 结构警报）")
print("=" * 100)
print("%-16s %10s %10s %8s %12s" % ("分子", "Lipinski", "Veber", "QED", "结构警报"))
for r in rows:
    nm = r["name"]
    hits = catalog.GetMatches(mols[nm])
    print("%-16s %10s %10s %8.3f %12s" % (nm, "%d 违规" % r["Lipinski"] if r["Lipinski"] else "通过",
          "通过" if r["Veber"] else "超标", r["QED"], "有" if hits else "无"))
