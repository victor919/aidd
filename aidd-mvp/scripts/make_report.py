"""把所有中间 json 汇总成一份可读的报告 reports/summary.md"""
import os
import json
import datetime

from _common import cfg, P, ensure_dirs, processed, path, log, banner

ensure_dirs()
banner("汇总报告")


def load(name):
    p = processed(name)
    if not os.path.exists(p):
        return None
    with open(p, encoding="utf-8") as f:
        return json.load(f)


d = {
    "target": load("01_target_rank.json"),
    "structure": load("02_structure_qc.json"),
    "pocket": load("pocket.json"),
    "target_id": load("04_target_identity.json"),
    "activity": load("04_activity_qc.json"),
    "curation": load("05_curation.json"),
    "baseline": load("07_baseline.json"),
    "chemprop": load("08_chemprop.json"),
    "docking": load("09_docking.json"),
    "admet": load("10_admet.json"),
    "closed_loop": load("11_closed_loop.json"),
}

missing = [k for k, v in d.items() if v is None]
if missing:
    log(f"以下环节尚未运行，报告中跳过：{missing}", "WARN")

L = []
w = L.append

w("# aidd-mvp · EGFR 端到端流水线运行汇总\n")
w(f"> 生成时间：{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}　"
  f"> 所有数字均为本机实跑结果\n")

# ---- 1 靶点 ----
if d["target"]:
    w("## ① 靶点排序\n")
    w(f"- Top1：**{d['target']['ranking'][0]['target']}**")
    if d["target"].get("top1_stable_under_weight_change") is not None:
        w(f"- 权重扰动下是否稳定：{d['target']['top1_stable_under_weight_change']}")
    w("")

# ---- 2/3 结构与口袋 ----
if d["pocket"]:
    w("## ②③ 蛋白结构与口袋\n")
    w(f"- PDB `{d['pocket'].get('pdb_id')}`，配体残基 `{d['pocket'].get('ligand')}`，"
      f"半径 {d['pocket'].get('radius')} 埃")
    res = d["pocket"].get("pocket_residues", [])
    w(f"- 口袋残基数：**{len(res)}**")
    key = [r for r in res if r.get("uniProt_num") in (718, 719, 745, 790, 793, 797, 844, 854, 855)]
    if key:
        w("- 命中关键位点：" + "、".join(
            f"{r['resname']}{r['uniProt_num']}" for r in sorted(key, key=lambda x: x["uniProt_num"])))
    w(f"- 编号偏移：1M17 = UniProt − **24**（结构 C773 = 临床 C797）\n")

# ---- 4/5 数据 ----
if d["activity"] and d["curation"]:
    w("## ④⑤ 活性数据获取与清洗\n")
    a, c = d["activity"], d["curation"]
    w(f"- ChEMBL `{d['target_id'].get('target_chembl_id', 'CHEMBL203')}`："
      f"{a.get('n_records')} 条活性 → {c.get('n_molecules')} 个分子")
    w(f"- 多测量分子 {c.get('n_multi_measured')} 个，其中位标准差 "
      f"{c.get('median_std_of_multi', 0):.3f} log")
    w(f"- pIC50 范围 {c['pIC50_range'][0]:.2f} ~ {c['pIC50_range'][1]:.2f}\n")

# ---- 6/7/8 模型 ----
if d["baseline"]:
    w("## ⑥⑦⑧ 模型对照（同一个 scaffold split）\n")
    w("| Split / 模型 | MAE | RMSE | R2 |")
    w("|---|---:|---:|---:|")
    b = d["baseline"]
    for k in ("random", "scaffold"):
        if k in b:
            v = b[k]
            w(f"| ECFP4+RF · {k} | {v['MAE']:.4f} | {v['RMSE']:.4f} | {v['R2']:.4f} |")
    if d["chemprop"]:
        cp = d["chemprop"]
        w(f"| Chemprop D-MPNN · scaffold | {cp['MAE']:.4f} | {cp['RMSE']:.4f} | {cp['R2']:.4f} |")
    if d["chemprop"]:
        w(f"\nΔMAE = {d['chemprop'].get('delta_MAE'):+.4f} → {d['chemprop'].get('verdict')}")
        w(f"（{d['chemprop'].get('epochs')} epoch，可训练参数 {d['chemprop'].get('n_params'):,}）\n")
    else:
        w("")

# ---- 9 对接 ----
if d["docking"]:
    w("## ⑨ 分子对接\n")
    rd = d["docking"].get("redock_check", {})
    w(f"- **自检（重对接 erlotinib）**：晶体自打分 {rd.get('crystal_pose_self_score')} kcal/mol，"
      f"最优 pose {rd.get('docked_best_score')} kcal/mol，**RMSD {rd.get('best_pose_rmsd')} 埃**")
    w(f"- 虚拟筛选：{d['docking'].get('n_docked')}/{d['docking'].get('n_screened')} 成功")
    w(f"- 预测 pIC50 与对接分数 Pearson r = **{d['docking'].get('pearson_pred_vs_docking')}**"
      f"（接近 0 说明两条证据正交，合并才有意义）")
    w(f"- Consensus 权重：活性 {d['docking']['weights']['activity']} / "
      f"对接 {d['docking']['weights']['docking']}\n")

# ---- 10 ADMET ----
if d["admet"]:
    w("## ⑩ ADMET 与候选\n")
    f = d["admet"].get("funnel", {})
    w(f"- 漏斗：{d['admet'].get('n_total')} → 活性 {f.get('activity')} → "
      f"+成药性 {f.get('plus_lipinski_veber')} → 最终 **{f.get('final_candidates')}**")
    pr = d["admet"].get("pass_rates", {})
    w(f"- 通过率：Lipinski {pr.get('Lipinski',0)*100:.1f}% · Veber {pr.get('Veber',0)*100:.1f}% · "
      f"QED {pr.get('QED',0)*100:.1f}% · 活性 {pr.get('activity',0)*100:.1f}%")
    w(f"- Pareto 前沿：**{d['admet'].get('pareto_front_size')}** 个分子")
    pe = d["admet"].get("pareto_extremes", {})
    if pe:
        w(f"  - 活性最强：pIC50 {pe['max_activity']['pIC50']} / QED {pe['max_activity']['QED']}")
        w(f"  - QED 最好：pIC50 {pe['max_qed']['pIC50']} / QED {pe['max_qed']['QED']}")
    w("")

# ---- 11 闭环 ----
if d["closed_loop"]:
    w("## ⑪ 干湿闭环\n")
    w("| 策略 | 起始 MAE | 最终 MAE | 最终 R2 | 样本数 |")
    w("|---|---:|---:|---:|---:|")
    for k, v in d["closed_loop"].get("summary", {}).items():
        w(f"| {k} | {v['start_MAE']:.4f} | {v['final_MAE']:.4f} | {v['final_R2']:.4f} | {v['final_n']} |")
    ntr = d["closed_loop"].get("n_to_reach_same_MAE", {})
    if ntr:
        w("\n达到同等 MAE 所需样本数：" + "，".join(f"{k} {v}" for k, v in ntr.items()))
    w("")

out = path(P["reports"], "summary.md")
with open(out, "w", encoding="utf-8") as f:
    f.write("\n".join(L) + "\n")
log(f"报告已写出：{os.path.relpath(out, path('.'))}")
print("\n".join(L))
