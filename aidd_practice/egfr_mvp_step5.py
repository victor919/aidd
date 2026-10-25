# 文章关联：AI制药实践篇B_实操篇_公众号文章.html
# 脚本职责：早期扩展：GCN 与 ECFP / RF 表示对照
# 使用注解：依赖 step3；不是本篇正式流水线的第五阶段，文章当前章节不能按文件编号直接对应。
# 完整对应关系与运行顺序：README.md
"""EGFR MVP Step5: 手写一个轻量 GCN，作为 ECFP+RF 的"图模型"对照（本机实跑）
说明：Chemprop(D-MPNN) 需要单独环境，本机未安装；这里用纯 PyTorch 手写一个
      GCN + mean pooling + MLP head，用于回答"换成图表示以后到底多学到了什么"。
"""
import os, json, time, warnings
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.model_selection import GroupShuffleSplit
from sklearn.metrics import mean_absolute_error, r2_score
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")
RDLogger.DisableLog("rdApp.*")

BASE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(BASE, "data", "egfr_mvp")
full = pd.read_csv(os.path.join(DATA, "dataset_with_scaffold.csv"))
full["scaffold"] = full["scaffold"].fillna("")
y = full["pIC50_median"].values

# ---------- 分子图特征 ----------
ATOM_LIST = list(range(1, 36))  # 原子序数 1..35
def atom_feat(a):
    return [
        ATOM_LIST.index(a.GetAtomicNum()) if a.GetAtomicNum() in ATOM_LIST else 0,
        int(a.GetDegree()),
        int(a.GetFormalCharge()) + 2,
        int(a.GetNumRadicalElectrons()),
        int(a.GetIsAromatic()),
        int(a.IsInRing()),
        int(a.GetHybridization()),
    ]
N_ATOM_F = 7
N_DEG, N_CHG, N_HYB = 11, 5, 8

def mol_to_graph(smiles):
    m = Chem.MolFromSmiles(smiles)
    n = m.GetNumAtoms()
    x = np.zeros((n, N_ATOM_F), dtype=np.float32)
    for i, a in enumerate(m.GetAtoms()):
        f = atom_feat(a)
        x[i, 0] = min(f[0], 34)
        x[i, 1] = min(f[1], N_DEG - 1)
        x[i, 2] = min(f[2], N_CHG - 1)
        x[i, 3] = min(f[3], 4)
        x[i, 4] = f[4]
        x[i, 5] = f[5]
        x[i, 6] = min(int(f[6]), N_HYB - 1)
    edges = []
    for b in m.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        edges.append((i, j)); edges.append((j, i))
    return x, edges

t0 = time.time()
graphs = [mol_to_graph(s) for s in full["std_smiles"]]
print(f"built {len(graphs)} molecular graphs in {time.time()-t0:.1f}s")
n_atoms = [g[0].shape[0] for g in graphs]
print(f"atoms per molecule: mean={np.mean(n_atoms):.1f} min={min(n_atoms)} max={max(n_atoms)}")

# ---------- 划分 ----------
splitter = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=42)
tr, te = next(splitter.split(full, groups=full["scaffold"]))
rng = np.random.RandomState(42)
perm = rng.permutation(len(full))
r_tr, r_te = perm[len(te):], perm[:len(te)]

# ---------- 模型 ----------
class GCN(nn.Module):
    def __init__(self, n_atom_feat, emb=64, hidden=128, n_layers=3, dropout=0.1):
        super().__init__()
        self.atom_emb = nn.Embedding(35, emb)
        self.proj = nn.Linear(n_atom_feat - 1, emb)
        self.convs = nn.ModuleList([nn.Linear(emb, emb) for _ in range(n_layers)])
        self.act = nn.ReLU()
        self.drop = nn.Dropout(dropout)
        self.head = nn.Sequential(nn.Linear(emb, hidden), nn.ReLU(),
                                  nn.Dropout(dropout), nn.Linear(hidden, 1))
        self.n_layers = n_layers

    def forward(self, x, adj_idx, batch):
        z = self.atom_emb(x[:, 0].long()) + self.proj(x[:, 1:])
        for conv in self.convs:
            m = torch.zeros_like(z)
            m.index_add_(0, adj_idx[0], z[adj_idx[1]])   # 邻域求和
            deg = torch.zeros(z.size(0), device=z.device).index_add_(0, adj_idx[0], torch.ones_like(adj_idx[1], dtype=z.dtype))
            m = m / deg.clamp(min=1).unsqueeze(1)
            z = self.act(conv(z + m))
            z = self.drop(z)
        g = torch.zeros(batch.max().item() + 1, z.size(1), device=z.device)
        g = g.index_add(0, batch, z)
        cnt = torch.zeros(g.size(0), device=z.device).index_add(0, batch, torch.ones(z.size(0), device=z.device))
        g = g / cnt.clamp(min=1).unsqueeze(1)
        return self.head(g).squeeze(-1)

def collate(idxs):
    xs, rows, cols, batch = [], [], [], []
    off = 0
    for k, i in enumerate(idxs):
        x, e = graphs[i]
        xs.append(torch.tensor(x))
        for (a, b) in e:
            rows.append(a + off); cols.append(b + off)
        batch.extend([k] * x.shape[0]); off += x.shape[0]
    return (torch.cat(xs, 0),
            torch.stack([torch.tensor(rows), torch.tensor(cols)]),
            torch.tensor(batch))

device = "mps" if torch.backends.mps.is_available() else "cpu"
print("device:", device)

def train_eval(train_idx, test_idx, tag, epochs=40, lr=1e-3):
    torch.manual_seed(42); np.random.seed(42)
    model = GCN(N_ATOM_F).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    lossf = nn.MSELoss()
    yt = torch.tensor(y[train_idx], dtype=torch.float32)
    ntr, nte = len(train_idx), len(test_idx)
    bs = 128
    t0 = time.time()
    for ep in range(epochs):
        model.train()
        order = np.random.RandomState(ep).permutation(ntr)
        tot = 0.0
        for s in range(0, ntr, bs):
            bidx = train_idx[order[s:s + bs]]
            x, adj, bat = collate(bidx)
            x, adj, bat = x.to(device), adj.to(device), bat.to(device)
            tgt = torch.tensor(y[bidx], dtype=torch.float32).to(device)
            opt.zero_grad()
            loss = lossf(model(x, adj, bat), tgt)
            loss.backward(); opt.step()
            tot += loss.item() * len(bidx)
        if (ep + 1) % 10 == 0 or ep == 0:
            print(f"      {tag} epoch {ep+1:>2d}  train MSE={tot/ntr:.4f}  ({time.time()-t0:.0f}s)", flush=True)
    model.eval()
    preds = []
    with torch.no_grad():
        for s in range(0, nte, 256):
            bidx = test_idx[s:s + 256]
            x, adj, bat = collate(bidx)
            p = model(x.to(device), adj.to(device), bat.to(device)).cpu().numpy()
            preds.append(p)
    p = np.concatenate(preds)
    mae = mean_absolute_error(y[test_idx], p)
    rmse = float(np.sqrt(np.mean((y[test_idx] - p) ** 2)))
    r2 = r2_score(y[test_idx], p)
    sp = spearmanr(y[test_idx], p).statistic
    print(f"  ==> {tag}: MAE={mae:.3f} RMSE={rmse:.3f} R2={r2:.3f} Spearman={sp:.3f} ({time.time()-t0:.0f}s)", flush=True)
    return dict(label=tag, mae=round(mae, 3), rmse=round(rmse, 3), r2=round(r2, 3), spearman=round(sp, 3))

print("\n=== GCN on random split ===")
r1 = train_eval(r_tr, r_te, "GCN / random split")
print("\n=== GCN on scaffold split ===")
r2_ = train_eval(tr, te, "GCN / scaffold split")

json.dump([r1, r2_], open(os.path.join(DATA, "gcn_results.json"), "w"), indent=2)
print("\nsaved -> gcn_results.json")
