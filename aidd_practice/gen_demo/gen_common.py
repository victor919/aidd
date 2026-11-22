# -*- coding: utf-8 -*-
"""⑬ 生成实验公共模块：局部提议、指纹/最近邻、属性、oracle 与记账。

铁律：A/B/C 三臂共用全部底层函数；唯一差别是"保留规则"。
"""
import os, json, random, sys, time
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import BRICS, Descriptors, QED, rdFingerprintGenerator
from rdkit import RDConfig

RDLogger.DisableLog("rdApp.*")
sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PROC = os.path.join(ROOT, "aidd-mvp", "data", "processed")
RES = os.path.join(HERE, "results")
GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


class World:
    def __init__(self, seed=42):
        self.df = pd.read_csv(os.path.join(PROC, "dataset_with_scaffold.csv"))
        self.y = self.df["pIC50_median"].values.astype(float)
        self.X = np.load(os.path.join(PROC, "X_ecfp.npy"))
        idx = json.load(open(os.path.join(RES, "e23_indices.json")))
        self.fit = np.array(idx["fit"]); self.test = np.array(idx["test"])
        self.smis = self.df["smiles"].values
        fp_cache = os.path.join(RES, "e23_fitfps.pkl")
        if os.path.exists(fp_cache):
            import pickle
            self.fit_fps = pickle.load(open(fp_cache, "rb"))
        else:
            self.fit_fps = [GEN.GetFingerprint(Chem.MolFromSmiles(self.smis[i])) for i in self.fit]
            import pickle
            pickle.dump(self.fit_fps, open(fp_cache, "wb"))
            print("cached fit fps", len(self.fit_fps), flush=True)
        frag = json.load(open(os.path.join(RES, "e23_fragments.json")))
        self.pieces = frag["pieces"]
        self.by_dummy = {int(k): v for k, v in frag["by_dummy"].items()}
        self.seed = seed
        self.pool = {}          # smiles -> dict(pred, nn, SA, QED, MW)
        self.query_count = 0
        self.cache_hits = 0
        self.invalid = 0

    def oracle_fit(self, seed=42):
        from sklearn.ensemble import RandomForestRegressor
        self.rf = RandomForestRegressor(n_estimators=500, random_state=seed, n_jobs=-1)
        self.rf.fit(self.X[self.fit], self.y[self.fit])
        # 训练用多核；单候选预测切回单进程，避免 joblib 每次启动的开销
        self.rf.n_jobs = 1

    def evaluate(self, smi):
        """一个新候选的评分：计入查询次数；重复候选走缓存不计费。"""
        if smi in self.pool:
            self.cache_hits += 1
            return self.pool[smi]
        m = Chem.MolFromSmiles(smi)
        fp = GEN.GetFingerprint(m)
        pred = float(self.rf.predict(np.asarray(fp, dtype=np.uint8).reshape(1, -1))[0])
        nn = float(max(DataStructs.BulkTanimotoSimilarity(fp, self.fit_fps)))
        d = dict(pred=round(pred, 4), nn=round(nn, 4),
                 SA=round(float(sascorer.calculateScore(m)), 3),
                 QED=round(float(QED.qed(m)), 3),
                 MW=round(float(Descriptors.MolWt(m)), 1),
                 atoms=m.GetNumHeavyAtoms())
        self.pool[smi] = d
        self.query_count += 1
        return d

    def propose(self, smi, rng):
        """局部动作：BRICS 片段替换（一次只换一个片段，保持接口数一致）。"""
        m = Chem.MolFromSmiles(smi)
        if m is None:
            self.invalid += 1; return None
        try:
            frags = list(BRICS.BRICSDecompose(m))
            if len(frags) < 2:
                self.invalid += 1; return None
            cand = [f for f in frags if 1 <= f.count("*") <= 2]
            if not cand:
                self.invalid += 1; return None
            out = rng.choice(cand)
            k = out.count("*")
            others = list(frags); others.remove(out)
            ins = None
            for _ in range(30):
                p = rng.choice(self.by_dummy.get(k, []))
                molp = Chem.MolFromSmiles(p)
                if molp is None: continue
                na = sum(1 for a in molp.GetAtoms() if a.GetSymbol() != "*")
                mo = Chem.MolFromSmiles(out)
                oa = sum(1 for a in mo.GetAtoms() if a.GetSymbol() != "*")
                if abs(na - oa) <= 12:
                    ins = p; break
            if ins is None:
                self.invalid += 1; return None
            fmols = [Chem.MolFromSmiles(x) for x in others + [ins]]
            fmols = [x for x in fmols if x is not None]
            built = None
            for b in BRICS.BRICSBuild(fmols, scrambleReagents=True):
                built = b; break
            if built is None:
                self.invalid += 1; return None
            Chem.SanitizeMol(built)
            s = Chem.MolToSmiles(built)
            m2 = Chem.MolFromSmiles(s)
            if (not s) or ("." in s) or m2 is None or m2.GetNumHeavyAtoms() > 60 \
                    or Descriptors.MolWt(m2) > 900:
                self.invalid += 1; return None
            return s
        except Exception:
            self.invalid += 1; return None


def stats(smis, world):
    if not smis:
        return dict(n=0)
    rows = [world.pool[s] for s in smis]
    nn = np.array([r["nn"] for r in rows]); pr = np.array([r["pred"] for r in rows])
    return dict(n=len(smis),
                pred_med=round(float(np.median(pr)), 4),
                nn_med=round(float(np.median(nn)), 4),
                frac_nn_lt_05=round(float((nn < .5).mean()), 4),
                frac_nn_lt_04=round(float((nn < .4).mean()), 4),
                SA_med=round(float(np.median([r["SA"] for r in rows])), 3),
                QED_med=round(float(np.median([r["QED"] for r in rows])), 3),
                MW_med=round(float(np.median([r["MW"] for r in rows])), 1),
                atoms_med=float(np.median([r["atoms"] for r in rows])))


def diversity(smis, world):
    if len(smis) < 2: return None
    fps = [GEN.GetFingerprint(Chem.MolFromSmiles(s)) for s in smis]
    vals = [DataStructs.TanimotoSimilarity(fps[i], fps[j])
            for i in range(len(fps)) for j in range(i + 1, len(fps))]
    return round(float(np.mean(vals)), 4)
