#!/usr/bin/env python
"""e24 · 阶段1：同池"随机保留 vs 按分数排序"的离线对照

诚实标注：这是**离线排序**，不是闭环爬山——池子是固定的，没有状态更新。
用途是先把"分数选择会不会把候选带离参考分布"这个问题在最简单的设置下问一遍。
预算：两臂都查询同一批候选（相同查询次数），A 忽略分数随机留，B 按分数留。
"""
import os, json, random, time
import numpy as np
from gen_common import World, stats, diversity

HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results")
T0 = time.time()
def log(*a): print(f"[{time.time()-T0:6.1f}s]", *a, flush=True)

w = World(); w.oracle_fit()
idx = json.load(open(os.path.join(RES, "e23_indices.json")))
starts = json.load(open(os.path.join(RES, "e23_pool.json")))["start_indices"]

# 固定候选池：每个起点独立提议 400 次（三臂共享），种子 42/43/44 各自生成同一批提议
POOL_PER_START = 400
KEPT = 100
res = {"design": dict(pool_per_start=POOL_PER_START, kept_per_start=KEPT,
                      seeds=[42, 43, 44], n_starts=len(starts),
                      note="离线排序：同池、同查询次数；A 随机保留、B 按 RF 预测保留；无状态更新，不称闭环爬山"),
       "runs": []}

for seed in [42, 43, 44]:
    for si, st in enumerate(starts):
        rng = random.Random(seed * 100 + si)
        cand = []
        for _ in range(POOL_PER_START):
            p = w.propose(w.smis[st], rng)
            if p: cand.append(p)
        cand = list(dict.fromkeys(cand))
        for s in cand: w.evaluate(s)
        rng2 = random.Random(seed * 1000 + si)
        keep_a = rng2.sample(cand, min(KEPT, len(cand)))
        keep_b = sorted(cand, key=lambda s: -w.pool[s]["pred"])[:KEPT]
        ra, rb = stats(keep_a, w), stats(keep_b, w)
        ra["diversity"] = diversity(keep_a, w); rb["diversity"] = diversity(keep_b, w)
        res["runs"].append(dict(seed=seed, start=si, n_pool=len(cand), A=ra, B=rb))
        log(f"seed{seed} start{si} pool={len(cand)} | A pred {ra['pred_med']:.3f} nn {ra['nn_med']:.3f} "
            f"| B pred {rb['pred_med']:.3f} nn {rb['nn_med']:.3f}")

res["queries_total"] = w.query_count
res["cache_hits"] = w.cache_hits
res["invalid_proposals"] = w.invalid
json.dump(res, open(os.path.join(RES, "e24_offline.json"), "w", encoding="utf-8"),
          ensure_ascii=False, indent=2)
log("wrote e24_offline.json; queries", w.query_count, "cache", w.cache_hits)
