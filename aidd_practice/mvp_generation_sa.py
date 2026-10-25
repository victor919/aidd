# 文章关联：AI制药实践篇B_实操篇_公众号文章.html
# 脚本职责：补充实验：BRICS 玩具生成与 SAScore 对照
# 使用注解：写真实分子属性缓存、mvp_generation_sa.json 及 t4；与 fast 版共享输出名但生成方法不同。
# 完整对应关系与运行顺序：README.md
"""验证「Day-1 就要把 SAScore 写进筛选」这句话。

做法（玩具级生成器，不是真的生成模型，但足以说明问题）：
  ① 把 ChEMBL 里 10468 个真分子做 BRICS 分解，收集片段；
  ② 随机挑片段重新拼回去，得到 3000 个「AI 生成」的分子；
  ③ 对比 真分子 vs 生成分子 的 SAScore / QED / MW 分布。

如果说 SAScore 对真实化合物库几乎无区分度，那它到底是在拦谁？
答案就在这组对照里。
"""
import os, json, time, random, warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, QED, BRICS
from rdkit import RDConfig

RDLogger.DisableLog("rdApp.*")

ROOT = os.path.expanduser("~/WorkBuddy/Work/aidd-mvp")
OUT = os.path.expanduser("~/WorkBuddy/Work/aidd_practice")
import sys
sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer  # noqa: E402

random.seed(42); np.random.seed(42)
T0 = time.time()
def log(m): print(f"[{time.time()-T0:6.1f}s] {m}", flush=True)

df = pd.read_csv(os.path.join(ROOT, "data", "processed", "dataset.csv"))
log(f"读入 {len(df)} 个真实分子")

# ---------- 1. BRICS 拆解真分子，攒片段池 ----------
# 注意：BRICSDecompose 默认返回 SMILES 字符串元组，不是 mol 对象
pieces, n_ok = set(), 0
for smi in df["smiles"]:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        continue
    try:
        frags = BRICS.BRICSDecompose(mol)
        for f in frags:
            if isinstance(f, str) and "*" in f:
                pieces.add(f)
        n_ok += 1
    except Exception:
        pass
pieces = list(pieces)
log(f"成功拆解 {n_ok} 个分子 → BRICS 片段池 {len(pieces)} 个（去重）")

def props(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    return dict(smiles=smi, SA=round(sascorer.calculateScore(mol), 2),
                QED=round(QED.qed(mol), 3), MW=round(Descriptors.MolWt(mol), 1))

# ---------- 2. 真分子属性 ----------
CACHE = os.path.join(OUT, "mvp_real_props.csv")
if os.path.exists(CACHE):
    real = pd.read_csv(CACHE)
    log(f"复用缓存的真分子属性 {len(real)} 个")
else:
    real = pd.DataFrame([p for p in (props(s) for s in df["smiles"]) if p])
    real.to_csv(CACHE, index=False)
    log(f"真分子属性计算完成 {len(real)} 个（已缓存）")

# ---------- 3. 玩具生成器：随机拼 2~4 个片段 ----------
TARGET_N = 3000
gen_rows, tries = [], 0
while len(gen_rows) < TARGET_N and tries < TARGET_N * 40:
    tries += 1
    break_ = False
    k = random.choice([2, 2, 2, 2, 3])
    picks = [random.choice(pieces) for _ in range(k)]
    # BRICS 拼接：把带连接点(*)的片段随机拼成一个完整分子
    frag_mols = [mm for mm in (Chem.MolFromSmiles(p) for p in picks) if mm is not None]
    if len(frag_mols) < 2:
        continue
    try:
        n_from_this = 0
        for built in BRICS.BRICSBuild(frag_mols):
            s = Chem.MolToSmiles(built)
            if s and "." not in s:
                p = props(s)
                if p:
                    gen_rows.append(p)
                    n_from_this += 1
            if len(gen_rows) >= TARGET_N:
                break_ = True
                break
            if n_from_this >= 3:   # 每组片段最多取 3 个，避免无限枚举
                break
    except Exception:
        continue
    if break_:
        break
gen = pd.DataFrame(gen_rows).drop_duplicates("smiles")
log(f"生成分子 {len(gen)} 个（尝试 {tries} 次）")

def summarize(d, name):
    out = dict(group=name, n=len(d),
               SA_mean=round(float(d["SA"].mean()), 2), SA_median=round(float(d["SA"].median()), 2),
               SA_max=round(float(d["SA"].max()), 2),
               pass_SA45=round(float((d["SA"] <= 4.5).mean()), 4),
               pass_SA60=round(float((d["SA"] <= 6.0).mean()), 4),
               QED_mean=round(float(d["QED"].mean()), 3),
               MW_mean=round(float(d["MW"].mean()), 1),
               pass_lipinski=round(float(((d["MW"] <= 500)).mean()), 4))
    return out

S = [summarize(real, "真实分子（ChEMBL，都被合成过）"), summarize(gen, "玩具生成器拼出来的分子")]
res = dict(piece_pool=len(pieces), summary=S,
           generated_examples=gen.nlargest(3, "SA")[["smiles", "SA", "QED", "MW"]].to_dict("records"))
print("\n" + "=" * 68)
for s in S:
    print(f"  【{s['group']}】 n={s['n']}")
    print(f"     SAScore 均值 {s['SA_mean']}  中位 {s['SA_median']}  最难 {s['SA_max']}")
    print(f"     SA<=4.5 通过率 {s['pass_SA45']:.1%}   SA<=6.0 通过率 {s['pass_SA60']:.1%}")
    print(f"     QED 均值 {s['QED_mean']}   平均分子量 {s['MW_mean']}")
print(f"\n  最难合成的生成分子样例：")
for r in res["generated_examples"]:
    print(f"     SA={r['SA']}  QED={r['QED']}  MW={r['MW']}  {r['smiles'][:70]}")

json.dump(res, open(os.path.join(OUT, "mvp_generation_sa.json"), "w", encoding="utf-8"),
          indent=2, ensure_ascii=False)
log(f"wrote {OUT}/mvp_generation_sa.json")

# ---------- 4. 图 ----------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Hiragino Sans GB"
plt.rcParams["axes.unicode_minus"] = False
GREEN, RED = "#0FA87E", "#C0392B"
fig, axes = plt.subplots(2, 1, figsize=(9.2, 6.6))
bins = np.linspace(1, 9, 45)
axes[0].hist(real["SA"], bins=bins, alpha=.72, color=GREEN, label=f"真实分子 n={len(real)}")
axes[0].hist(gen["SA"], bins=bins, alpha=.72, color=RED, label=f"生成分子 n={len(gen)}")
axes[0].axvline(4.5, color="#333", ls="--", lw=1.5)
axes[0].text(4.6, axes[0].get_ylim()[1] * .8, "SA = 4.5 常用阈值", fontsize=9)
axes[0].set_xlabel("SAScore（越大越难合成）"); axes[0].set_ylabel("分子数")
axes[0].set_title(f"SAScore 在拦谁：真实分子 {S[0]['pass_SA45']:.0%} 通过，生成分子只有 {S[1]['pass_SA45']:.0%} 通过",
                  fontsize=12)
axes[0].legend(fontsize=10)

axes[1].hist(real["QED"], bins=40, alpha=.72, color=GREEN, label="真实分子")
axes[1].hist(gen["QED"], bins=40, alpha=.72, color=RED, label="生成分子")
axes[1].set_xlabel("QED（成药性，越大越好）"); axes[1].set_ylabel("分子数")
axes[1].set_title(f"同一批生成分子的成药性：QED 均值 {S[1]['QED_mean']} vs 真实 {S[0]['QED_mean']}", fontsize=12)
axes[1].legend(fontsize=10)
fig.tight_layout()
FIG = os.path.expanduser("~/WorkBuddy/Work/AI制药实操篇_配图")
fig.savefig(f"{FIG}/t4_generation_sa.png", dpi=170)
log("t4_generation_sa.png")
