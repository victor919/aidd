#!/usr/bin/env python
"""e27 · 借鉴参考稿的诊断视角，用⑬自己的数据补跑（全部本机实跑）。

1) 逐对 B-A 差值（同起点同种子 9 对，全部列出，不挑方向）
2) 终点身份回查：终点规范 SMILES 是否重现拟合池（结构身份，非只看 Tanimoto）
3) 动作偏置：按臂统计分子量/SA/QED 的逐轮变化与"保留当前"比例
4) 产出统计：槽位/尝试/有效/唯一/相对拟合池与完整快照的新分子比例
5) 属性筛选代价：对终点按 MW/SA/QED 逐层过滤的保留数（事后演示，不进搜索）
输出 results/e27_audit_extras.json
"""
import os, json, collections
import numpy as np
from rdkit import Chem

HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results")
w25 = json.load(open(os.path.join(RES, "e25_closedloop.json")))
runs = w25["runs"]
recs = w25["all_records"]           # 修复后的全量保留记录（27 次运行）
e23 = json.load(open(os.path.join(RES, "e23_pool.json")))
out = {}

# 拟合池与完整快照的规范 SMILES 集合
import pandas as pd
df = pd.read_csv(os.path.join(HERE, "..", "..", "aidd-mvp", "data", "processed", "dataset_with_scaffold.csv"))
fit_smis = set(df.smiles.values[e23["fit_indices"]] if "fit_indices" in e23 else [])
if not fit_smis:
    idx = json.load(open(os.path.join(RES, "e23_indices.json")))
    fit_smis = set(df.smiles.values[np.array(idx["fit"])])
all_smis = set(df.smiles.values)

# ---------- 1) 逐对 B-A ----------
by = collections.defaultdict(dict)
for r in runs:
    by[(r["seed"], r["start"])][r["arm"]] = r["rounds"]
pairs = []
for (s, st), arms in sorted(by.items()):
    if "A" not in arms or "B" not in arms: continue
    a, b = arms["A"], arms["B"]
    da_pred = b[-1]["pred_med"] - b[0]["pred_med"] - (a[-1]["pred_med"] - a[0]["pred_med"])
    da_nn = b[-1]["nn_med"] - b[0]["nn_med"] - (a[-1]["nn_med"] - a[0]["nn_med"])
    pairs.append(dict(seed=s, start=st, d_pred_BminusA=round(da_pred, 4),
                      d_nn_BminusA=round(da_nn, 4),
                      B_nn_lower_than_A=bool(b[-1]["nn_med"] < a[-1]["nn_med"])))
out["pairwise_B_minus_A"] = dict(
    n=len(pairs),
    d_pred_median=round(float(np.median([p["d_pred_BminusA"] for p in pairs])), 4),
    d_nn_median=round(float(np.median([p["d_nn_BminusA"] for p in pairs])), 4),
    B_lower_count=sum(p["B_nn_lower_than_A"] for p in pairs),
    rows=pairs)

# ---------- 2) 终点身份回查 ----------
endpoints = []
for r in runs:
    rr = r["rounds"]
    if not rr: continue
    # 该运行最后保留的候选（从 all_records 中按 seed/start/arm 与最大 round 取）
    sub = [x for x in recs if x["seed"] == r["seed"] and x["start"] == r["start"] and x["arm"] == r["arm"]]
    if not sub: continue
    mx = max(x["round"] for x in sub)
    last = [x for x in sub if x["round"] == mx]
    endpoints.append(dict(seed=r["seed"], start=r["start"], arm=r["arm"],
                          n_last=len(last),
                          smis=[x["smiles"] for x in last],
                          in_fit=[x["smiles"] in fit_smis for x in last],
                          nn=[x["nn"] for x in last]))
byarm = {}
for arm in "ABC":
    es = [e for e in endpoints if e["arm"] == arm]
    in_fit_any = sum(any(e["in_fit"]) for e in es)
    nn1_exact = sum(any(x in fit_smis for x in e["smis"]) and all(n >= .9999 for n in e["nn"]) for e in es)
    byarm[arm] = dict(n_runs=len(es), runs_with_fit_endpoint=in_fit_any,
                      runs_with_nn1_endpoint=nn1_exact,
                      endpoint_unique=len(set(s for e in es for s in e["smis"])))
out["endpoint_identity"] = dict(fit_smis=len(fit_smis), by_arm=byarm,
                                note="身份回查用规范 SMILES 集合精确匹配；NN≈1 只作辅助，不单凭相似度下结论")

# ---------- 3) 动作偏置 ----------
bias = {}
for arm in "ABC":
    rs = [r for r in runs if r["arm"] == arm]
    drops = [r["rounds"][-1]["MW_med"] - r["rounds"][0]["MW_med"] for r in rs if len(r["rounds"]) >= 2]
    sas = [r["rounds"][-1]["SA_med"] - r["rounds"][0]["SA_med"] for r in rs if len(r["rounds"]) >= 2]
    bias[arm] = dict(n=len(drops),
                     MW_change_median=round(float(np.median(drops)), 1),
                     MW_change_range=[round(min(drops), 1), round(max(drops), 1)],
                     SA_change_median=round(float(np.median(sas)), 3),
                     QED_change_median=round(float(np.median([
                         r["rounds"][-1]["QED_med"] - r["rounds"][0]["QED_med"]
                         for r in rs if len(r["rounds"]) >= 2])), 3))
# 重访/回退比例：子记录中相邻轮次 SMILES 重复出现（用 all_records）
revisit = {}
for arm in "ABC":
    sub = [x for x in recs if x["arm"] == arm]
    byrun = collections.defaultdict(list)
    for x in sub: byrun[(x["seed"], x["start"])].append(x)
    tot = rep = 0
    for k, xs in byrun.items():
        seen = set()
        for x in sorted(xs, key=lambda z: z["round"]):
            tot += 1
            if x["smiles"] in seen: rep += 1
            seen.add(x["smiles"])
    revisit[arm] = dict(slots=tot, repeats=rep, frac=round(rep / tot, 4) if tot else 0)
out["action_bias"] = dict(by_arm=bias, revisit=revisit,
                          note="随机臂无选择压力，片段替换近似随机游走；分子量/SA 变化是动作本身的偏置，不是活性结论")

# ---------- 4) 产出统计（每臂独立） ----------
prod = {}
for arm in "ABC":
    sub = [x for x in recs if x["arm"] == arm]
    smis = set(x["smiles"] for x in sub)
    prod[arm] = dict(slots=len(sub), unique=len(smis),
                     novel_vs_fit=round(sum(s not in fit_smis for s in smis) / max(len(smis), 1), 4),
                     novel_vs_snapshot=round(sum(s not in all_smis for s in smis) / max(len(smis), 1), 4))
out["production"] = dict(by_arm=prod, queries_total=w25["queries_total"],
                         cache_hits=w25["cache_hits"], invalid=w25["invalid_proposals"],
                         note="槽位含重复；分母不同，不能连成通关率")

# ---------- 5) 属性筛选代价（终点去重后，事后演示） ----------
def canon(s):
    m = Chem.MolFromSmiles(s)
    return Chem.MolToSmiles(m) if m else s
filt = {}
for arm in "ABC":
    es = [e for e in endpoints if e["arm"] == arm]
    smis = sorted(set(s for e in es for s in e["smis"]))
    attrs = [x for x in recs if x["arm"] == arm and x["smiles"] in smis]
    seen = {}
    for x in attrs: seen[x["smiles"]] = x   # 取每条记录的最后属性（同 smiles 属性一致）
    vals = list(seen.values())
    mw = [x for x in vals if 150 <= x["MW"] <= 600]
    mwsa = [x for x in mw if x["SA"] <= 4.5]
    mwsaqed = [x for x in mwsa if x["QED"] >= 0.4]
    filt[arm] = dict(endpoint_unique=len(smis), mw=len(mw), mw_sa=len(mwsa), mw_sa_qed=len(mwsaqed))
out["attr_filter_demo"] = dict(criteria="MW 150-600; SA<=4.5; QED>=0.4（搜索中未使用，事后演示）",
                               by_arm=filt,
                               note="剩得多不等于活性更好；阈值不参与证明搜索优劣")

json.dump(out, open(os.path.join(RES, "e27_audit_extras.json"), "w", encoding="utf-8"),
          ensure_ascii=False, indent=2)
print("wrote e27_audit_extras.json")
for k in out:
    if k in ("endpoint_identity", "production", "attr_filter_demo"):
        print(k, json.dumps(out[k], ensure_ascii=False)[:600])
    elif k == "pairwise_B_minus_A":
        print(k, out[k]["n"], "d_pred_med", out[k]["d_pred_median"], "d_nn_med", out[k]["d_nn_median"],
              "B lower", out[k]["B_lower_count"])
    else:
        print(k, json.dumps(out[k], ensure_ascii=False)[:600])
