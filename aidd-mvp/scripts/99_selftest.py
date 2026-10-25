"""99 · 自检：把"踩过的坑"固化成断言，让它们不会再悄悄复发。

用法：
    python scripts/99_selftest.py            # 只读检查（秒级），跑在流水线之后
    python scripts/99_selftest.py --full     # 追加"可复现性"重算（约 1 分钟）

退出码：0 = 全通过；1 = 有 FAIL。

为什么需要它：
    这条流水线上每一个断言，背后都对应一个真实踩过的坑——
    编号差 24、受体里混进共晶配体、随机切分造成的"默写"、
    chemprop 归一化写反、docking 与活性正交……
    这些错误都不会让程序报错，只会让结论悄悄错掉。
    文档会过时，注释会被忽略，只有断言会被 CI 执行。
"""
import os
import sys
import json
import argparse

import numpy as np
import pandas as pd

from _common import cfg, raw, processed, log, banner

FULL = "--full" in sys.argv

RESULTS = []


def claim(cid, desc, ok, detail=""):
    """记录一条断言。ok=True/False/None(跳过)"""
    RESULTS.append((cid, desc, ok, detail))
    mark = "PASS" if ok is True else ("SKIP" if ok is None else "FAIL")
    print(f"  [{mark}] {cid}  {desc}" + (f"  ->  {detail}" if detail else ""))
    return ok


def main():
    banner("99 · Selftest：数据契约 / 科学不变量 / 可复现性")

    def need(p):
        if not os.path.exists(p):
            print(f"\n缺少产物 {os.path.basename(p)}，请先跑完流水线（bash run_all.sh）")
            raise SystemExit(1)
        return p

    ds = processed("dataset.csv")
    need(ds)
    df = pd.read_csv(ds)

    print("\n── A. 数据契约 ─────────────────────────────────────────")
    claim("A1", "样本量落在预期区间（8000~20000）",
          8000 <= len(df) <= 20000, f"{len(df)} 个分子")
    claim("A2", "SMILES 唯一且非空",
          bool(df.smiles.is_unique and df.smiles.notna().all()),
          f"重复 {int(len(df) - df.smiles.nunique())} 个")
    y = df["pIC50_median"]
    # 只守上界：弱活性（pIC50 低至 1.6，即 IC50 = 25 mM）是真实数据，不该被砍；
    # 上界 11 对应 IC50 = 0.01 nM，超过它基本是单位录错（清洗前见过 17.30）。
    claim("A3", "pIC50 上界 ≤ 11（超出的必是单位录错）",
          bool((y <= 11).all()),
          f"实际 [{y.min():.2f}, {y.max():.2f}]（清洗前曾出现 17.30）")
    claim("A3b", "pIC50 上界截断确实生效（最大值恰为 11，而非自然分布的尾部）",
          bool(abs(y.max() - 11.0) < 1e-6),
          f"max = {y.max():.4f}")

    dws = processed("dataset_with_scaffold.csv")
    X_path = processed("X_ecfp.npy")
    splits = json.load(open(need(processed("split_indices.json"))))
    if os.path.exists(dws) and os.path.exists(X_path):
        d2 = pd.read_csv(dws)
        X = np.load(X_path, mmap_mode="r")
        claim("A4", "特征矩阵与标签表行数一致", len(d2) == len(X),
              f"{len(d2)} 行 / X{X.shape}")
    for name, sp in splits.items():
        tr, te = set(sp["train"]), set(sp["test"])
        claim(f"A5.{name}", f"{name} 切分：训练/测试无交集且覆盖全集",
              bool(not (tr & te) and (tr | te) == set(range(len(df)))),
              f"train {len(tr)} / test {len(te)}")

    print("\n── B. 切分泄漏（这条流水线最贵的一个坑） ──────────────")
    if os.path.exists(dws):
        d2 = pd.read_csv(dws)
        te, tr = splits["scaffold"]["test"], splits["scaffold"]["train"]
        ov = set(d2.iloc[tr]["scaffold"].fillna("")) & set(d2.iloc[te]["scaffold"].fillna(""))
        claim("B1", "骨架切分下训练集/测试集骨架零重叠",
              len(ov) == 0, f"重叠 {len(ov)} 个骨架")
        rr_te, rr_tr = splits["random"]["test"], splits["random"]["train"]
        ovr = set(d2.iloc[rr_tr]["scaffold"].fillna("")) & set(d2.iloc[rr_te]["scaffold"].fillna(""))
        claim("B2", "随机切分确实存在骨架重叠（对照组，说明对比有意义）",
              len(ovr) > 0, f"重叠 {len(ovr)} 个骨架（这正是'默写'的来源）")

    print("\n── C. 结构与口袋（差 24 的编号偏移） ──────────────────")
    pk = json.load(open(need(processed("pocket.json"))))
    num = pk["numbering"]
    claim("C1", "自动检测出的编号偏移 == config 里写的偏移",
          num["detected_offset"] == num["config_offset"],
          f"detected={num['detected_offset']}  config={num['config_offset']}")
    claim("C2", "序列比对匹配率 ≥ 99%",
          num["match_rate"] >= 0.99, f"{num['matched']} = {num['match_rate']:.1%}")
    canon = pk["canonical_sites_in_pocket"]
    claim("C3", "共价残基 C797 在口袋内（PDB 773 + 24）",
          bool(canon.get("C797 共价残基")),
          "命中" if canon.get("C797 共价残基") else "未命中")
    # 盒子中心必须等于共晶配体质心（伪代码硬编码坐标的坑）
    pdb_path = raw(f"{cfg['structure']['pdb_id']}.pdb")
    if os.path.exists(pdb_path):
        lig = cfg["structure"]["ligand_resname"]
        pts = []
        for line in open(pdb_path, encoding="utf-8"):
            if line.startswith("HETATM") and line[17:20].strip() == lig:
                pts.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        cen = np.array(pts).mean(axis=0)
        box = json.load(open(processed("03_pocket.json")))["box"] \
            if os.path.exists(processed("03_pocket.json")) else None
        box = box or json.load(open(os.path.join(os.path.dirname(pdb_path),
                                                 "..", "interim", "docking_box.json")))["center"] \
            if os.path.exists(os.path.join(os.path.dirname(pdb_path),
                                           "..", "interim", "docking_box.json")) else None
        if box is not None:
            d = float(np.linalg.norm(cen - np.array(box)))
            claim("C4", "对接口中心 == 共晶配体质心（误差 < 0.6 埃）",
                  d < 0.6, f"偏差 {d:.3f} 埃（硬编码坐标会偏差 31 埃）")

    print("\n── D. 对接自证（受体里混进配体，是最隐蔽的坑） ────────")
    dk_path = processed("09_docking.json")
    if os.path.exists(dk_path):
        dk = json.load(open(dk_path))
        rc = dk["redock_check"]
        claim("D1", "晶体构象自打分 < 0（受体已剔除共晶配体）",
              rc["crystal_pose_self_score"] < 0,
              f"{rc['crystal_pose_self_score']} kcal/mol（含配体时会是 +438）")
        claim("D2", "重对接 RMSD < 2.0 埃",
              rc["best_pose_rmsd"] < 2.0, f"{rc['best_pose_rmsd']} 埃")
    else:
        claim("D1", "对接产物存在", None, "未运行 09，跳过")

    print("\n── E. 模型与闭环（负结果也要守住） ────────────────────")
    bl = json.load(open(need(processed("07_baseline.json"))))
    claim("E1", "骨架切分 MAE 不劣化到不可用（< 1.0）",
          bl["scaffold"]["MAE"] < 1.0, f"scaffold MAE = {bl['scaffold']['MAE']}")
    claim("E2", "骨架切分 MAE ≥ 随机切分 MAE（切分变难是预期方向）",
          bl["scaffold"]["MAE"] >= bl["random"]["MAE"],
          f"random {bl['random']['MAE']} -> scaffold {bl['scaffold']['MAE']}")
    if os.path.exists(dk_path):
        dk = json.load(open(dk_path))
        r = abs(dk["pearson_pred_vs_docking"])
        claim("E3", "预测活性与对接分接近正交（|r| < 0.5，守住'两者不可互相替代'）",
              r < 0.5, f"|r| = {r:.3f}")
    hist = pd.read_csv(need(processed("iteration_history.csv")))
    last = hist.sort_values("n_labeled").groupby("strategy").tail(1).set_index("strategy")["MAE"]
    if {"uncertainty", "greedy"} <= set(last.index):
        claim("E4", "闭环终点 uncertainty 不劣于 greedy",
              float(last["uncertainty"]) <= float(last["greedy"]),
              f"uncertainty {last['uncertainty']:.3f} vs greedy {last['greedy']:.3f}")

    if FULL:
        print("\n── F. 可复现性（重算一遍，而不是读缓存） ────────────")
        from rdkit import Chem, DataStructs, RDLogger
        from rdkit.Chem import rdFingerprintGenerator
        from rdkit.Chem.Scaffolds import MurckoScaffold
        from sklearn.model_selection import GroupShuffleSplit
        RDLogger.DisableLog("rdApp.*")
        m = cfg["model"]
        gen = rdFingerprintGenerator.GetMorganGenerator(
            radius=int(m["fingerprint_radius"]), fpSize=int(m["fingerprint_size"]))
        n = 500
        want = np.load(X_path, mmap_mode="r")[:n]
        got = []
        for smi in df["smiles"].tolist()[:n]:
            arr = np.zeros((int(m["fingerprint_size"]),), dtype=np.uint8)
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                DataStructs.ConvertToNumpyArray(gen.GetFingerprint(mol), arr)
            got.append(arr)
        claim("F1", "同参数重算前 500 个分子 ECFP，与产物逐位一致",
              bool(np.array_equal(np.stack(got), np.array(want))),
              "逐位相同" if np.array_equal(np.stack(got), np.array(want)) else "不一致")

        d2 = pd.read_csv(dws)
        scaf = d2["scaffold"].fillna("").values
        gss = GroupShuffleSplit(n_splits=1, test_size=float(m["test_size"]),
                                random_state=int(cfg["random_seed"]))
        tr2, te2 = next(gss.split(d2, groups=scaf))
        claim("F2", "同种子重算骨架切分，索引与产物一致",
              bool(np.array_equal(np.sort(tr2), np.array(sorted(splits["scaffold"]["train"])))),
              f"train {len(tr2)} / test {len(te2)}")

    print("\n" + "=" * 72)
    npass = sum(1 for r in RESULTS if r[2] is True)
    nfail = sum(1 for r in RESULTS if r[2] is False)
    nskip = sum(1 for r in RESULTS if r[2] is None)
    print(f"{npass} PASS / {nfail} FAIL / {nskip} SKIP")
    if nfail:
        print("失败项：" + ", ".join(r[0] for r in RESULTS if r[2] is False))
        return 1
    print("PASS：流水线产物自洽，踩过的坑都还锁着。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
