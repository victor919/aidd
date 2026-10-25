"""01 · 靶点排序（模拟靶点发现逻辑）

输入：data/raw/target_evidence.csv（首次运行自动生成一份教学用证据表）
输出：data/processed/target_rank.csv + reports/01_target_rank.json

要点：靶点选择本质是多证据决策，不是"模型输出 0.97 就结束"。
"""
import os
import pandas as pd

from _common import cfg, ensure_dirs, raw, processed, log, banner, save_json, exists

ensure_dirs()
banner("01 · Target Ranking：多证据加权排序")

# ---------- 1. 证据表（真实项目里应来自 Open Targets / DepMap / 文献） ----------
ev_path = raw("target_evidence.csv")
if not exists(ev_path):
    rows = [
        # target, disease, genetics, pathway, druggability, safety
        ("EGFR",    0.95, 0.82, 0.90, 0.95, 0.30),
        ("KRAS",    0.88, 0.90, 0.86, 0.35, 0.45),
        ("Target_B",0.70, 0.91, 0.65, 0.50, 0.45),
        ("Target_C",0.80, 0.55, 0.82, 0.35, 0.25),
        ("Target_D",0.60, 0.40, 0.55, 0.80, 0.10),
    ]
    df = pd.DataFrame(rows, columns=["target", "disease", "genetics", "pathway",
                                     "druggability", "safety"])
    df.to_csv(ev_path, index=False)
    log(f"生成教学用证据表 -> {os.path.basename(ev_path)}（真实项目请替换为 Open Targets 导出）")
else:
    log(f"复用已有证据表 {os.path.basename(ev_path)}")

df = pd.read_csv(ev_path)
log(f"候选靶点 {len(df)} 个：{list(df['target'])}")

# ---------- 2. 加权评分 ----------
W = dict(disease=0.30, genetics=0.25, pathway=0.20, druggability=0.25, safety=-0.20)
df["score"] = (
    W["disease"] * df["disease"]
    + W["genetics"] * df["genetics"]
    + W["pathway"] * df["pathway"]
    + W["druggability"] * df["druggability"]
    + W["safety"] * df["safety"]        # safety 是风险，负权重
)
df = df.sort_values("score", ascending=False).reset_index(drop=True)
df["rank"] = df.index + 1

print("\n" + df[["rank", "target", "score"]].to_string(index=False,
      float_format=lambda x: f"{x:.3f}"))

out = processed("target_rank.csv")
df.to_csv(out, index=False)
log(f"wrote data/processed/target_rank.csv")

# ---------- 3. 敏感性分析：权重一变，排序还稳吗？ ----------
# 这是教材里不会写、但真实项目必做的一步
alt = {
    "默认（可成药性 0.25）":    dict(disease=.30, genetics=.25, pathway=.20, druggability=.25, safety=-.20),
    "只看生物学证据":            dict(disease=.40, genetics=.35, pathway=.25, druggability=.00, safety=-.00),
    "只看可成药性+安全":         dict(disease=.10, genetics=.10, pathway=.10, druggability=.50, safety=-.20),
}
print("\n权重敏感性分析（Top1 会不会变？）：")
sens = {}
for name, w in alt.items():
    s = (w["disease"] * df["disease"] + w["genetics"] * df["genetics"]
         + w["pathway"] * df["pathway"] + w["druggability"] * df["druggability"]
         + w["safety"] * df["safety"])
    top1 = df.loc[s.idxmax(), "target"]
    sens[name] = top1
    print(f"  {name:<22s} -> Top1 = {top1}")
stable = len(set(sens.values())) == 1
print(f"\n  => 排序对权重{'不敏感（稳健）' if stable else '敏感（需要更多证据）'}")

save_json({
    "script": "01_target_rank.py",
    "n_candidates": int(len(df)),
    "ranking": df[["rank", "target", "score"]].to_dict(orient="records"),
    "weights": W,
    "sensitivity": sens,
    "top1_stable_under_weight_change": bool(stable),
    "selected_target": cfg["target"]["gene"],
}, processed("01_target_rank.json"))

banner("01 done · 选中靶点 = %s（UniProt %s）"
       % (cfg["target"]["gene"], cfg["target"]["uniprot"]))
