# -*- coding: utf-8 -*-
# 文章关联：AI制药实践篇_公众号文章.html
# 脚本职责：实验：第七站，DMTA 增量更新与 LoRA 对照
# 使用注解：TDC Lipophilicity 骨架切分；配图 make_figures5.py 使用手工录入的实测摘要，不自动读取本次输出。
# 完整对应关系与运行顺序：README.md
"""干湿闭环 · 增量微调对照实验（本机实测）

场景（模拟真实 DMTA 循环）：
  模型先在「历史公开数据」上预训练；随后每一轮湿实验回传一批「自家新系列」分子。
  对比三种更新方式：(A) 不更新  (B) 全参微调  (C) LoRA 低秩适配器（rank r）
关键控制：三种方式使用完全相同的学习率、迭代轮数与数据，唯一差别是「哪些参数能动」。

数据：TDC Lipophilicity_AstraZeneca（4200 个真实分子，实验测定 logD7.4）
切分：按 Bemis-Murcko 骨架做 scaf fold 切分 → 新系列与历史系列骨架不相交
计算：PyTorch 2.13 + RDKit 2026.03.5
"""
import warnings, numpy as np, pandas as pd
warnings.filterwarnings("ignore")
import torch, torch.nn as nn
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from tdc.single_pred import ADME
RDLogger.DisableLog('rdApp.*')

torch.manual_seed(0)
rng = np.random.default_rng(0)

# ---------- 1. 数据 ----------
data = ADME(name="Lipophilicity_AstraZeneca").get_data()
print("数据集 Lipophilicity_AstraZeneca：%d 个真实分子，标签 = 实验测定 logD7.4" % len(data))
print("标签范围 %.2f ~ %.2f，均值 %.2f，标准差 %.2f" % (data.Y.min(), data.Y.max(), data.Y.mean(), data.Y.std()))

gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
X, y, scaffolds = [], [], []
for smi, val in zip(data.Drug, data.Y):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        continue
    try:
        scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=m, includeChirality=False)
    except Exception:
        scaf = ""
    X.append(np.array(gen.GetFingerprint(m), dtype=np.float32))
    y.append(float(val)); scaffolds.append(scaf)
X = np.array(X); y = np.array(y); scaffolds = np.array(scaffolds)
print("可解析分子 %d 个，指纹 %d 维，不同 Murcko 骨架 %d 个" % (len(y), X.shape[1], len(set(scaffolds))))

# ---------- 2. 骨架切分：历史系列 vs 自家新系列 ----------
groups = {}
for i, s in enumerate(scaffolds):
    groups.setdefault(s, []).append(i)
groups = sorted(groups.items(), key=lambda kv: -len(kv[1]))
new_idx, old_idx = [], []
for s, idxs in groups:
    (new_idx if len(new_idx) < 700 else old_idx).extend(idxs)
new_idx = np.array(sorted(new_idx)); old_idx = np.array(sorted(old_idx))
rng.shuffle(new_idx); rng.shuffle(old_idx)
old_test, old_train = old_idx[:400], old_idx[400:]
new_test, new_pool = new_idx[:300], new_idx[300:]
print("历史训练 %d | 历史留出 %d | 新系列候选池 %d | 新系列留出 %d"
      % (len(old_train), len(old_test), len(new_pool), len(new_test)))

mu, sd = y[old_train].mean(), y[old_train].std()

# ---------- 3. 模型：带 LoRA 插槽的 MLP ----------
class MLP(nn.Module):
    def __init__(self, d_in, h=256, r=8, alpha=16.0):
        super().__init__()
        self.fc1 = nn.Linear(d_in, h); self.fc2 = nn.Linear(h, h); self.fc3 = nn.Linear(h, 1)
        self.scaling = alpha / r
        self.lora_A1 = nn.Parameter(torch.zeros(r, d_in))
        self.lora_B1 = nn.Parameter(torch.zeros(h, r))
        self.lora_A2 = nn.Parameter(torch.zeros(r, h))
        self.lora_B2 = nn.Parameter(torch.zeros(h, r))
        self.use_lora = False

    def forward(self, x):
        h1 = self.fc1(x)
        h2 = self.fc2(h1)
        if self.use_lora:
            h1 = h1 + (x @ self.lora_A1.T @ self.lora_B1.T) * self.scaling
            h2 = self.fc2(h1) + (h1 @ self.lora_A2.T @ self.lora_B2.T) * self.scaling
        return self.fc3(torch.relu(h2)).squeeze(-1)

    def set_trainable(self, mode):
        names = {"all": ["fc1", "fc2", "fc3"],
                 "lora": ["lora_A1", "lora_B1", "lora_A2", "lora_B2"],
                 "none": []}[mode]
        for n, p in self.named_parameters():
            p.requires_grad = any(n.startswith(k) for k in names)
        if mode == "lora":
            nn.init.normal_(self.lora_A1, std=0.01); nn.init.normal_(self.lora_A2, std=0.01)
            nn.init.zeros_(self.lora_B1); nn.init.zeros_(self.lora_B2)
            self.use_lora = True
        else:
            self.use_lora = False
        return sum(p.numel() for p in self.parameters() if p.requires_grad)

    def n_total(self):
        return sum(p.numel() for p in self.parameters())

Xt = torch.tensor(X, dtype=torch.float32)

def rmse(model, idx):
    model.eval()
    with torch.no_grad():
        p = model(Xt[idx]).numpy() * sd + mu
    return float(np.sqrt(np.mean((p - y[idx]) ** 2)))

def fit(model, rows, epochs=40, lr=1e-3):
    if not any(p.requires_grad for p in model.parameters()):
        return model
    model.train()
    opt = torch.optim.Adam([p for p in model.parameters() if p.requires_grad], lr=lr)
    xs = torch.tensor(X[rows], dtype=torch.float32)
    ys = torch.tensor((y[rows] - mu) / sd, dtype=torch.float32)
    lf = nn.MSELoss()
    for _ in range(epochs):
        opt.zero_grad(); lf(model(xs), ys).backward(); opt.step()
    return model

def pretrain(seed):
    torch.manual_seed(seed)
    m = MLP(X.shape[1])
    m.set_trainable("all")
    fit(m, old_train, epochs=60, lr=1e-3)
    return m

# ---------- 4. 累积式 DMTA 轮次 ----------
ROUNDS = [8, 20, 50, 120, 250]      # 累计回传的新系列分子数
SEEDS = [0, 1, 2]
STRATS = ["不更新", "全参微调", "LoRA r=8"]
records = []
for seed in SEEDS:
    order = rng.permutation(new_pool)
    for strat in STRATS:
        m = pretrain(seed)
        base_new, base_old = rmse(m, new_test), rmse(m, old_test)
        for total in ROUNDS:
            labeled = np.array(sorted(order[:total]))
            mode = {"不更新": "none", "全参微调": "all", "LoRA r=8": "lora"}[strat]
            m_curr = pretrain(seed)
            if mode != "none":
                n_tr = m_curr.set_trainable(mode)
                # 从头开始，只在这一轮累计到的数据上训练（公平：同 lr、同轮数）
                fit(m_curr, labeled, epochs=40, lr=1e-3)
            else:
                n_tr = 0
            records.append(dict(seed=seed, strat=strat, total=total,
                                new=rmse(m_curr, new_test), old=rmse(m_curr, old_test),
                                trainable=n_tr, base_new=base_new, base_old=base_old))

H = pd.DataFrame(records)
TOTAL = pretrain(0).n_total()
agg = H.groupby(["strat", "total"])[["new", "old", "trainable"]].mean().reset_index()

print()
print("=" * 108)
print("【实测】每一轮湿实验累计回传 n 个自家分子后：新系列学得怎么样？旧知识还剩多少？（RMSE，越小越好）")
print("=" * 108)
print("%-12s %10s %16s %14s %18s" % ("更新方式", "累计回传", "新系列 RMSE(Δ)", "历史 RMSE(Δ)", "可训练参数"))
b0 = agg[(agg.strat == "不更新")].iloc[0]
for s in STRATS:
    for t in ROUNDS:
        r = agg[(agg.strat == s) & (agg.total == t)].iloc[0]
        if s == "不更新":
            print("%-12s %10d %16.3f %14.3f %18s" % (s, t, r.new, r.old, "0"))
        else:
            print("%-12s %10d %15.3f(%+.3f) %13.3f(%+.3f) %17d" %
                  (s, t, r.new, r.new - b0.new, r.old, r.old - b0.old, r.trainable))

print()
print("模型总参数 %d；LoRA 适配器 r=8 的可训练参数 %d（占 %.2f%%）。"
      % (TOTAL, agg[(agg.strat == "LoRA r=8")]["trainable"].iloc[0],
         100 * agg[(agg.strat == "LoRA r=8")]["trainable"].iloc[0] / TOTAL))
print("说明：LoRA 的 B 矩阵以 0 初始化，故训练起点与预训练模型完全一致；三种方式使用相同的 lr=1e-3、40 轮。")

# ---------- 5. 同一底座挂多个系列适配器 ----------
print()
print("=" * 108)
print("【实测】一个底座 + 多个系列适配器：能不能一份模型服务多条管线？")
print("=" * 108)
base = pretrain(0)
import copy, pickle
base_bytes = len(pickle.dumps(base.state_dict()))
size = base_bytes / 1024.0
print("底座权重大小：%.1f KB（冻结，所有系列共用）" % size)
rows = []
for k, sid in enumerate([0, 1, 2]):
    order = rng.permutation(new_pool)
    labeled = np.array(sorted(order[:50]))
    ada = copy.deepcopy(base)
    ada.set_trainable("lora")
    fit(ada, labeled, epochs=40, lr=1e-3)
    rows.append(dict(series="系列 %d" % (k + 1),
                     rmse_new=rmse(ada, new_test), rmse_old=rmse(ada, old_test),
                     kb=len(pickle.dumps({k: v for k, v in ada.state_dict().items() if "lora" in k})) / 1024.0,
                     same_base=bool(np.allclose(base.fc1.weight.detach().numpy(), ada.fc1.weight.detach().numpy()))))
R = pd.DataFrame(rows)
print("%-10s %16s %16s %14s %12s" % ("适配器", "本系列 RMSE", "历史 RMSE", "适配器大小", "底座是否未变"))
for _, r in R.iterrows():
    print("%-10s %16.3f %16.3f %13.1f KB %12s" % (r.series, r.rmse_new, r.rmse_old, r.kb, "是" if r.same_base else "否"))

# ---------- 6. 适配器也要调参：学习率敏感性 ----------
print()
print("=" * 108)
print("【实测】适配器自身的学习率很关键：lr 太大，所谓「低秩保护」并不存在")
print("=" * 108)
print("%-16s %10s %12s %18s %18s" % ("更新方式", "累计回传", "学习率", "新系列 RMSE(Δ)", "历史 RMSE(Δ)"))
for lr in [1e-3, 3e-4, 1e-4]:
    for total in [50, 250]:
        ns, os_ = [], []
        for seed in SEEDS:
            order = rng.permutation(new_pool)
            labeled = np.array(sorted(order[:total]))
            m1 = pretrain(seed)
            n_tr = m1.set_trainable("lora")
            fit(m1, labeled, epochs=40, lr=lr)
            ns.append(rmse(m1, new_test)); os_.append(rmse(m1, old_test))
        print("%-16s %10d %12s %17.3f(%+.3f) %17.3f(%+.3f)"
              % ("LoRA r=8", total, "%.0e" % lr, np.mean(ns), np.mean(ns) - b0.new,
                 np.mean(os_), np.mean(os_) - b0.old))
for total in [50, 250]:
    ns, os_ = [], []
    for seed in SEEDS:
        order = rng.permutation(new_pool)
        labeled = np.array(sorted(order[:total]))
        m1 = pretrain(seed); m1.set_trainable("all")
        fit(m1, labeled, epochs=40, lr=1e-3)
        ns.append(rmse(m1, new_test)); os_.append(rmse(m1, old_test))
    print("%-16s %10d %12s %17.3f(%+.3f) %17.3f(%+.3f)"
          % ("全参微调", total, "1e-3", np.mean(ns), np.mean(ns) - b0.new,
             np.mean(os_), np.mean(os_) - b0.old))

# ---------- 7. 第四站：只按分数取头部，化学多样性会坍缩 ----------
print()
print("=" * 108)
print("【实测】苗头名单的多样性：只按预测分数「掐尖」，骨架会迅速坍缩")
print("=" * 108)
model0 = pretrain(0)
model0.eval()
with torch.no_grad():
    pred_all = model0(Xt).numpy()
pool = new_pool  # 整个新系列候选池当作「待筛选虚拟库」（骨架与训练集不相交）
scaf_pool = np.array([MurckoScaffold.MurckoScaffoldSmiles(mol=Chem.MolFromSmiles(s), includeChirality=False)
                      for s in data.Drug.values[pool]])
order_by_score = pool[np.argsort(-pred_all[pool])]
print("%-22s %12s %16s" % ("取法", "分子数", "不同 Murcko 骨架"))
for k in [20, 50, 100]:
    top = order_by_score[:k]
    n_scaf = len(set(scaf_pool[np.searchsorted(pool, top)] if (np.all(np.diff(pool) > 0)) else
                     set(scaf_pool[[list(pool).index(i) for i in top]]) ))
    # 贪心多样性重排：先按分数取候选，再逐个挑选与已选不同骨架的
    cand = order_by_score[:k * 5]
    chosen, seen = [], set()
    cmap = {idx: scaf_pool[list(pool).index(idx)] for idx in cand}
    for idx in cand:
        s = cmap[idx]
        if s in seen:
            continue
        seen.add(s); chosen.append(idx)
        if len(chosen) == k:
            break
    n_div = len(set([cmap[i] for i in chosen]))
    print("%-22s %12d %16s" % ("只按预测分数掐尖 top%d" % k, k, "%d 个" % n_scaf))
    print("%-22s %12d %16s" % ("同候选 + 多样性重排 top%d" % k, k, "%d 个（占 %.0f%%）" % (n_div, 100 * n_div / k)))
