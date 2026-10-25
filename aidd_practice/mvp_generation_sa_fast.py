# 文章关联：AI制药实践篇B_实操篇_公众号文章.html
# 脚本职责：补充实验：SMILES 拼接玩具生成对照
# 使用注解：先由 mvp_generation_sa.py 准备 mvp_real_props.csv；会覆盖同名 JSON 和 t4，不可混称 BRICS 结果。
# 完整对应关系与运行顺序：README.md
"""快速版生成对照：SMILES 字符串重组 toy 生成器。

比 BRICSBuild 粗暴得多（直接把两个真分子的 SMILES 各剪一刀接起来），
但化学上成立的分子更多、更怪，正适合验证「无约束生成会把 SAScore 推高」。
真分子属性复用 mvp_real_props.csv 缓存。
"""
import os, json, time, random, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit import RDConfig
import sys
sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer
from rdkit.Chem import Descriptors, QED

RDLogger.DisableLog("rdApp.*")
random.seed(42); np.random.seed(42)
T0 = time.time()
def log(m): print(f"[{time.time()-T0:6.1f}s] {m}", flush=True)

ROOT = os.path.expanduser("~/WorkBuddy/Work/aidd-mvp")
OUT = os.path.expanduser("~/WorkBuddy/Work/aidd_practice")
df = pd.read_csv(os.path.join(ROOT, "data", "processed", "dataset.csv"))
smis = df["smiles"].tolist()
log(f"读入 {len(smis)} 个真实分子")

def props(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None or mol.GetNumAtoms() < 8:
        return None
    try:
        return dict(smiles=smi, SA=round(sascorer.calculateScore(mol), 2),
                    QED=round(QED.qed(mol), 3), MW=round(Descriptors.MolWt(mol), 1),
                    rings=mol.GetRingInfo().NumRings())
    except Exception:
        return None

CACHE = os.path.join(OUT, "mvp_real_props.csv")
real = pd.read_csv(CACHE)
log(f"复用真分子属性 {len(real)} 个")

# ---------- 字符串重组 ----------
TARGET_N = 3000
gen_rows, tries = [], 0
seen = set(real["smiles"])
while len(gen_rows) < TARGET_N and tries < TARGET_N * 60:
    tries += 1
    s1, s2 = random.choice(smis), random.choice(smis)
    if len(s1) < 12 or len(s2) < 12:
        continue
    c1 = random.randrange(6, len(s1) - 4)
    c2 = random.randrange(6, len(s2) - 4)
    smi = s1[:c1] + s2[c2:]
    if smi in seen:
        continue
    p = props(smi)
    if p:
        gen_rows.append(p)
        seen.add(smi)
gen = pd.DataFrame(gen_rows)
log(f"生成 {len(gen)} 个分子（尝试 {tries} 次）")

def summarize(d, name):
    return dict(group=name, n=len(d),
                SA_mean=round(float(d["SA"].mean()), 2), SA_median=round(float(d["SA"].median()), 2),
                SA_max=round(float(d["SA"].max()), 2),
                pass_SA45=round(float((d["SA"] <= 4.5).mean()), 4),
                pass_SA60=round(float((d["SA"] <= 6.0).mean()), 4),
                QED_mean=round(float(d["QED"].mean()), 3),
                MW_mean=round(float(d["MW"].mean()), 1))

S = [summarize(real, "真实分子"), summarize(gen, "重组分子")]
res = dict(method="SMILES 字符串重组（toy）", summary=S,
           worst_examples=gen.nlargest(3, "SA")[["smiles", "SA", "QED", "MW"]].to_dict("records"))
print("=" * 64)
for s in S:
    print(f"  【{s['group']}】 n={s['n']}  SA 均值 {s['SA_mean']}  中位 {s['SA_median']}")
    print(f"      SA<=4.5 通过率 {s['pass_SA45']:.1%}   SA<=6.0 通过率 {s['pass_SA60']:.1%}   QED 均值 {s['QED_mean']}")
print("  最怪样例：")
for r in res["worst_examples"]:
    print(f"    SA={r['SA']}  {r['smiles'][:66]}")
json.dump(res, open(os.path.join(OUT, "mvp_generation_sa.json"), "w", encoding="utf-8"),
          indent=2, ensure_ascii=False)
log(f"wrote {OUT}/mvp_generation_sa.json")

# ---------- 图 ----------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Hiragino Sans GB"
plt.rcParams["axes.unicode_minus"] = False
GREEN, RED = "#0FA87E", "#C0392B"
FIG = os.path.expanduser("~/WorkBuddy/Work/AI制药实操篇_配图")
fig, axes = plt.subplots(2, 1, figsize=(9.2, 6.6))
bins = np.linspace(1, 9, 45)
axes[0].hist(real["SA"], bins=bins, alpha=.72, color=GREEN, label=f"真实分子 n={len(real)}")
axes[0].hist(gen["SA"], bins=bins, alpha=.72, color=RED, label=f"无约束重组分子 n={len(gen)}")
axes[0].axvline(4.5, color="#333", ls="--", lw=1.5)
axes[0].text(4.62, axes[0].get_ylim()[1] * .82, "SA = 4.5 常用阈值", fontsize=9)
axes[0].set_xlabel("SAScore（越大越难合成）"); axes[0].set_ylabel("分子数")
axes[0].set_title(f"SAScore 在拦谁：真实分子 {S[0]['pass_SA45']:.0%} 通过阈值，重组分子只有 {S[1]['pass_SA45']:.0%} 通过",
                  fontsize=12)
axes[0].legend(fontsize=10)
axes[1].hist(real["QED"], bins=40, alpha=.72, color=GREEN, label="真实分子")
axes[1].hist(gen["QED"], bins=40, alpha=.72, color=RED, label="无约束重组分子")
axes[1].set_xlabel("QED（成药性，越大越好）"); axes[1].set_ylabel("分子数")
axes[1].set_title(f"同一批重组分子的成药性：QED 均值 {S[1]['QED_mean']} vs 真实 {S[0]['QED_mean']}", fontsize=12)
axes[1].legend(fontsize=10)
fig.tight_layout()
fig.savefig(f"{FIG}/t4_generation_sa.png", dpi=170)
log("t4_generation_sa.png saved")
