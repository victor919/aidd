"""10 · ADMET 与成药性：第二个漏斗

输入：data/processed/dataset_with_scaffold.csv、predictions_scaffold.csv、docking.csv
输出：data/processed/admet_panel.csv、candidate_ranking.csv、10_admet.json

关键教学点：
    1) 活性强 ≠ 能成药，必须用第二个漏斗再过滤；
    2) 没有「唯一最佳分子」，只有不同取舍下的 Pareto 前沿；
    3) 公开数据库里存在物理上不可能的离群值（pIC50 > 11 意味着 IC50 < 0.01 nM），
       必须先体检再建模，否则 Pareto 前沿会被一个脏数据独占。
"""
import os
import sys
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen, QED, rdMolDescriptors
from rdkit.Chem import rdFingerprintGenerator
from rdkit import DataStructs

from _common import cfg, ensure_dirs, processed, log, banner, save_json

RDLogger.DisableLog("rdApp.*")
ensure_dirs()
banner("10 · ADMET 与成药性：从 Hit 到 Candidate 的第二个漏斗")

A = cfg["admet"]

# ---------------- 合成可及性 SA（RDKit Contrib） ----------------
try:
    from rdkit.Chem import RDConfig
    sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
    import sascorer
    _HAS_SA = True
except Exception:
    _HAS_SA = False
    log("sascorer 不可用，SA 列将留空", "WARN")

try:
    from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
    _params = FilterCatalogParams()
    _params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    _pains = FilterCatalog(_params)
    _HAS_PAINS = True
except Exception:
    _HAS_PAINS = False
    log("PAINS 不可用，跳过警报筛选", "WARN")


def calc_props(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    d = {
        "MW": Descriptors.MolWt(mol),
        "logP": Crippen.MolLogP(mol),
        "TPSA": Descriptors.TPSA(mol),
        "HBD": rdMolDescriptors.CalcNumHBD(mol),
        "HBA": rdMolDescriptors.CalcNumHBA(mol),
        "RotB": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "AromRings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "QED": QED.qed(mol),
    }
    d["SA"] = sascorer.calculateScore(mol) if _HAS_SA else np.nan
    if _HAS_PAINS:
        entry = _pains.GetFirstMatch(mol)
        d["PAINS"] = entry.GetDescription() if entry is not None else ""
    else:
        d["PAINS"] = ""
    return d


# ---------------- 1. 数据体检：先看有没有物理上不可能的值 ----------------
banner("体检：公开活性数据不是干净的")
full = pd.read_csv(processed("dataset_with_scaffold.csv"))
y = full["pIC50_median"].values
worst = full.loc[y.argmax()]
log(f"pIC50 范围 = {y.min():.2f} ~ {y.max():.2f}（均值 {y.mean():.2f}）")
# pIC50 = 11 意味着 IC50 = 1 pM。虽然 05 已经替我们挡掉了 >11 的极端值，
# 但处于边界上的这些点仍然值得标注：生化 assay 很少能可靠测到这个量级，
# 通常是数值/单位填写错误的结果。真正的项目里应该回去翻原始文献。
SUSPECT = 10.5
suspicious = full[y >= SUSPECT].copy()
log(f"pIC50 >= {SUSPECT} 的分子数 = {len(suspicious)}"
    f"（即 IC50 <= {10**(9-SUSPECT):.0e} nM，生化 assay 很少能可靠测到这个量级）")
for _, r in suspicious.head(3).iterrows():
    ic50 = 10 ** (9 - r["pIC50_median"])
    log(f"  例：pIC50={r['pIC50_median']:.2f} → IC50={ic50:.1e} nM "
        f"（n_measurements={int(r.get('n_measurements', -1))}）")
impossible = len(suspicious)
log("处理策略：建模保留，但候选列表里明确标注/剔除，避免一个脏数据独占 Pareto 前沿。")


# ---------------- 2. 计算成药性描述符面板 ----------------
banner(f"计算描述符面板（{len(full)} 个分子）")
rows = []
for i, smi in enumerate(full["smiles"]):
    p = calc_props(smi)
    if p is None:
        p = {k: np.nan for k in ["MW", "logP", "TPSA", "HBD", "HBA", "RotB",
                                 "AromRings", "QED", "SA"]}
        p["PAINS"] = ""
    rows.append(p)
    if (i + 1) % 2000 == 0:
        log(f"  {i+1}/{len(full)}")
props = pd.DataFrame(rows)
panel = pd.concat([full[["smiles", "scaffold", "pIC50_median", "n_measurements"]],
                   props], axis=1)

# Lipinski / Veber
panel["Lipinski_pass"] = ((panel.MW <= A["mw_max"]) & (panel.logP <= A["logp_max"]) &
                          (panel.HBD <= 5) & (panel.HBA <= 10))
panel["Veber_pass"] = (panel.RotB <= 10) & (panel.TPSA <= A["tpsa_max"])
panel["QED_pass"] = panel.QED >= A["qed_min"]
panel["activity_pass"] = panel.pIC50_median >= A["activity_min"]
panel["clean"] = panel.pIC50_median < SUSPECT          # 排除可疑高活性
panel.to_csv(processed("admet_panel.csv"), index=False)
log(f"Lipinski 通过率 {panel.Lipinski_pass.mean()*100:.1f}%  |  "
    f"Veber {panel.Veber_pass.mean()*100:.1f}%  |  QED>{A['qed_min']} {panel.QED_pass.mean()*100:.1f}%  |  "
    f"pIC50>={A['activity_min']} {panel.activity_pass.mean()*100:.1f}%")


# ---------------- 3. 三段漏斗：活性 → 成药性 → 候选 ----------------
banner("漏斗：活性 → 成药性 → 候选")
step1 = int(panel.activity_pass.sum())
step2 = int((panel.activity_pass & panel.Lipinski_pass & panel.Veber_pass).sum())
step3 = int((panel.activity_pass & panel.Lipinski_pass & panel.Veber_pass &
             panel.QED_pass & panel.clean).sum())
print(f"\n{'阶段':<34}{'剩余分子':>10}")
print("-" * 46)
print(f"{'全部分子':<34}{len(panel):>10}")
print(f"{'① pIC50 >= ' + str(A['activity_min']):<34}{step1:>10}")
print(f"{'② + Lipinski + Veber':<34}{step2:>10}")
print(f"{'③ + QED 门槛 + 剔除可疑高活性':<34}{step3:>10}")


# ---------------- 4. Pareto：没有唯一最好分子 ----------------
banner("Pareto 前沿：不要 sort_values('score') 取 Top1")
# 目标：最大化 pIC50、最大化 QED；两者经常冲突
sub = panel[panel.clean].copy()
sub = sub.dropna(subset=["pIC50_median", "QED"])


def pareto_front(df, cols, maxs):
    """返回不被任何其他点支配的索引"""
    vals = df[cols].values.astype(float)
    # 统一转成「越大越好」
    sign = np.array([1.0 if m else -1.0 for m in maxs])
    V = vals * sign
    keep = []
    for i in range(len(V)):
        dominated = False
        for j in range(len(V)):
            if i == j:
                continue
            if (V[j] >= V[i]).all() and (V[j] > V[i]).any():
                dominated = True
                break
        if not dominated:
            keep.append(i)
    return df.iloc[keep]


front = pareto_front(sub, ["pIC50_median", "QED"], [True, True])
log(f"Pareto 前沿上共有 {len(front)} 个分子（占 {len(front)/max(len(sub),1)*100:.2f}%）")
print(f"\n前沿上活性最高 / QED 最高 / 中间折中的三个：")
show_cols = ["smiles", "pIC50_median", "QED", "MW", "logP", "SA"]
f2 = front.sort_values("pIC50_median", ascending=False)
for tag, r in [("活性最强", f2.iloc[0]),
               ("QED 最好", f2.sort_values("QED", ascending=False).iloc[0]),
               ("折中", f2.iloc[len(f2) // 2])]:
    print(f"  {tag}: pIC50={r.pIC50_median:.2f}  QED={r.QED:.2f}  "
          f"MW={r.MW:.0f}  logP={r.logP:.2f}  SA={r.SA if pd.notna(r.SA) else float('nan'):.2f}")

print("\n  结论：活性最强的分子往往 QED 不好；QED 漂亮的又常常不够强。"
      "候选清单应该是「一组带视角的选项」，而不是一个「最佳分子」。")


# ---------------- 5. 候选清单：四种视角各出一列 ----------------
banner("候选清单：输出四组而不是一个 Winner")
cands = panel[panel.clean & panel.Lipinski_pass & panel.Veber_pass].copy()
picks = {}


def top_n(df, col, n=5, asc=False):
    return df.sort_values(col, ascending=asc).head(n)


groups = {
    "Top Activity": top_n(cands, "pIC50_median", 5),
    "Top Balanced": top_n(cands.assign(_b=cands.pIC50_median * 0.6 + cands.QED * 4), "_b", 5),
    "Top Druglike": top_n(cands, "QED", 5),
    "Top Low-risk": top_n(cands[cands.PAINS == ""], "MW", 5, asc=True),
}
for name, g in groups.items():
    for k, (_, r) in enumerate(g.iterrows()):
        picks.setdefault(name, []).append({
            "smiles": r.smiles, "pIC50": round(float(r.pIC50_median), 2),
            "QED": round(float(r.QED), 3), "MW": round(float(r.MW), 1),
            "logP": round(float(r.logP), 2), "TPSA": round(float(r.TPSA), 1),
            "SA": None if pd.isna(r.SA) else round(float(r.SA), 2),
            "PAINS": r.PAINS})

out_rows = []
for name, items in picks.items():
    for it in items:
        out_rows.append({"group": name, **it})
pd.DataFrame(out_rows).to_csv(processed("candidate_ranking.csv"), index=False)
log(f"候选清单已写出：{len(out_rows)} 条 / 4 组视角")

save_json({
    "n_total": int(len(panel)),
    "suspect_pic50_threshold": SUSPECT, "suspect_pic50_count": int(impossible),
    "funnel": {"activity": step1, "plus_lipinski_veber": step2, "final_candidates": step3},
    "pass_rates": {
        "Lipinski": round(float(panel.Lipinski_pass.mean()), 4),
        "Veber": round(float(panel.Veber_pass.mean()), 4),
        "QED": round(float(panel.QED_pass.mean()), 4),
        "activity": round(float(panel.activity_pass.mean()), 4),
    },
    "pareto_front_size": int(len(front)),
    "pareto_extremes": {
        "max_activity": {"pIC50": round(float(f2.iloc[0].pIC50_median), 2),
                         "QED": round(float(f2.iloc[0].QED), 3)},
        "max_qed": {"pIC50": round(float(f2.sort_values('QED', ascending=False).iloc[0].pIC50_median), 2),
                    "QED": round(float(f2.sort_values('QED', ascending=False).iloc[0].QED), 3)},
    },
    "candidate_groups": {k: len(v) for k, v in picks.items()},
}, processed("10_admet.json"))
banner(f"10 done · 漏斗终 {step3} 个候选，Pareto 前沿 {len(front)} 个")
