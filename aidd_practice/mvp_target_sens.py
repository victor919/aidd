# 文章关联：AI制药实践篇B_实操篇_公众号文章.html
# 脚本职责：补充实验：靶点排序权重敏感性
# 使用注解：读取正式工程 target_rank.csv，输出 mvp_target_sensitivity.json。
# 完整对应关系与运行顺序：README.md
"""01 的补充：靶点排序的权重敏感性（完整排名，不只看 Top1）

同一张证据表，换三套权重，看名次怎么抖。
输入：aidd-mvp/data/processed/target_rank.csv
"""
import os, json
import pandas as pd

ROOT = os.path.expanduser("~/WorkBuddy/Work/aidd-mvp")
OUT = os.path.expanduser("~/WorkBuddy/Work/aidd_practice")
df = pd.read_csv(os.path.join(ROOT, "data", "processed", "target_rank.csv"))
print(df.to_string(index=False))

schemas = {
    "baseline（均匀兼顾）": dict(disease=.30, genetics=.25, pathway=.20, druggability=.25, safety=-.20),
    "只看生物学证据": dict(disease=.40, genetics=.35, pathway=.25, druggability=.00, safety=.00),
    "可成药性优先": dict(disease=.15, genetics=.10, pathway=.10, druggability=.55, safety=-.10),
    "安全优先": dict(disease=.20, genetics=.15, pathway=.10, druggability=.15, safety=-.40),
}

res = {}
print("\n" + "=" * 72)
for name, w in schemas.items():
    s = (w["disease"] * df["disease"] + w["genetics"] * df["genetics"] +
         w["pathway"] * df["pathway"] + w["druggability"] * df["druggability"] +
         w["safety"] * df["safety"])
    order = df.assign(score=s.round(4)).sort_values("score", ascending=False)
    rank = list(order["target"])
    res[name] = {"weights": w,
                 "ranking": [{"target": t, "score": float(sc)} for t, sc in zip(order["target"], order["score"])]}
    print(f"\n[{name}]")
    for i, (t, sc) in enumerate(zip(order["target"], order["score"]), 1):
        print(f"   {i}. {t:<9} {sc:.3f}")

# 名次稳定性：每个靶点在不同方案下的名次范围
rows = []
for t in df["target"]:
    ranks = [ [x["target"] for x in res[k]["ranking"]].index(t) + 1 for k in schemas ]
    scores = [ [x["score"] for x in res[k]["ranking"]][[x["target"] for x in res[k]["ranking"]].index(t)] for k in schemas ]
    rows.append(dict(target=t, best_rank=min(ranks), worst_rank=max(ranks),
                     score_min=round(min(scores),3), score_max=round(max(scores),3)))
stab = pd.DataFrame(rows).sort_values("best_rank")
print("\n" + "=" * 72)
print("名次稳定性（四套权重下）：")
print(stab.to_string(index=False))

json.dump(res, open(os.path.join(OUT, "mvp_target_sensitivity.json"), "w", encoding="utf-8"),
          indent=2, ensure_ascii=False)
print(f"\nwrote {OUT}/mvp_target_sensitivity.json")
