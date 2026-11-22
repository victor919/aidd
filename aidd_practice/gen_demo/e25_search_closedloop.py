#!/usr/bin/env python
"""e25 · 阶段2：三臂闭环搜索（局部提议 + 状态更新）

A 随机保留 / B 按 RF 预测保留 / C 满足 NN>=0.4 后按预测保留（可选）。
三臂共用：同一起点、同一提议函数、同一每轮候选数、同一基本化学检查、同一查询预算。
唯一差别 = 保留规则。预算 = oracle 查询次数（新候选才计费，重复命中缓存不计费）。

输出：results/e25_closedloop.json（逐轮指标 + 逐候选审计记录）
"""
import os, json, random, time
import numpy as np
from gen_common import World, stats, diversity

HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results")
T0 = time.time()
def log(*a): print(f"[{time.time()-T0:6.1f}s]", *a, flush=True)

w = World(); w.oracle_fit()
pool_meta = json.load(open(os.path.join(RES, "e23_pool.json")))
starts = pool_meta["start_indices"]
SEEDS = [42, 43, 44]
QUERY_BUDGET = 400          # 每臂每次运行的评分器查询上限
BEAM = 5                    # 每轮保留数
PER_BEAM = 10               # 每个保留分子每轮提议数 → 每轮 50 次查询
AD_FLOOR = 0.4              # C 臂约束（⑪ 域外分箱下沿，仅作本实验规则）

res = {"design": dict(query_budget=QUERY_BUDGET, beam=BEAM, per_beam=PER_BEAM,
                      rounds=QUERY_BUDGET // (BEAM * PER_BEAM), seeds=SEEDS,
                      n_starts=len(starts), ad_floor=AD_FLOOR,
                      budget_unit="oracle 查询次数（新唯一有效候选才计费；重复走缓存）",
                      note="A/B/C 同起点、同提议函数、同每轮候选数、同基本检查；只改保留规则"),
       "runs": []}

for seed in SEEDS:
    for si, st in enumerate(starts):
        for arm in ["A", "B", "C"]:
            rng = random.Random(seed * 10000 + si * 100)  # 三臂同一随机流：提议相同，唯一差别是保留规则
            beam = [w.smis[st]]           # 同一起点
            rows, keep_records = [], []
            query0 = w.query_count  # 本臂查询起点（每臂独立预算记账）
            # 起点本身不消耗查询（来自拟合池，已有标签），但记录其身份
            for rnd in range(1, res["design"]["rounds"] + 1):
                cands = []
                for parent in beam:
                    for _ in range(PER_BEAM):
                        p = w.propose(parent, rng)
                        if p: cands.append((parent, p))
                # 去重（同轮内）；缓存命中不计费
                seen = set(); uniq = []
                for parent, p in cands:
                    if p in seen: continue
                    seen.add(p); uniq.append((parent, p))
                new = [(p, w.evaluate(p)) for parent, p in uniq]
                if not new: break
                scored = [(s, d) for s, d in new]
                if arm == "A":
                    chosen = rng.sample(scored, min(BEAM, len(scored)))
                elif arm == "B":
                    chosen = sorted(scored, key=lambda x: -x[1]["pred"])[:BEAM]
                else:
                    allowed = [x for x in scored if x[1]["nn"] >= AD_FLOOR]
                    chosen = sorted(allowed, key=lambda x: -x[1]["pred"])[:BEAM]
                    if not chosen:  # 无可保留：记停滞，不强行放宽
                        chosen = []
                for s, d in chosen:
                    keep_records.append(dict(round=rnd, arm=arm, seed=seed, start=si,
                                             smiles=s, pred=d["pred"], nn=d["nn"],
                                             SA=d["SA"], QED=d["QED"], MW=d["MW"]))
                m = stats([s for s, _ in chosen], w)
                m.update(round=rnd, queries=w.query_count - query0)
                rows.append(m)
                beam = [s for s, _ in chosen]
                if not beam: break
            res["runs"].append(dict(seed=seed, start=si, arm=arm, rounds=rows))
            last = rows[-1] if rows else {}
            log(f"seed{seed} start{si} {arm}: rounds={len(rows)} queries={w.query_count-query0} "
                f"last pred={last.get('pred_med','-')} nn={last.get('nn_med','-')}")

res["queries_total"] = w.query_count
res["cache_hits"] = w.cache_hits
res["invalid_proposals"] = w.invalid
res["audit_records"] = keep_records
json.dump(res, open(os.path.join(RES, "e25_closedloop.json"), "w", encoding="utf-8"),
          ensure_ascii=False, indent=2)
log("wrote e25_closedloop.json; records", len(keep_records), "queries", w.query_count,
    "cache", w.cache_hits, "invalid", w.invalid)
