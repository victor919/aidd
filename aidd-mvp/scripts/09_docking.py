"""09 · Structure-based screening：AutoDock Vina 分子对接

输入：data/raw/1M17.pdb、data/processed/predictions_scaffold.csv
输出：data/processed/docking.csv、09_docking.json、data/interim/*.pdbqt

本脚本先做一件被大多数人跳过的事——**验证自己的对接 Pipeline 是否可信**：
拿共晶配体 erlotinib（1M17 的 AQ4）重新对接，看能不能找回晶体 pose。
如果这一步 RMSD 很大，后面几万个分子的打分都不值得信。
"""
import os
import sys
import json
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

from _common import cfg, P, ensure_dirs, raw, interim, processed, path, log, banner, save_json, exists

RDLogger.DisableLog("rdApp.*")
ensure_dirs()
banner("09 · Structure-based：Vina 对接（先自证，再筛分子）")

# meeko / vina 在工程内 vendor 或已装好的环境里
from meeko.polymer import Polymer, ResidueChemTemplates
from meeko import MoleculePreparation
from meeko.writer import PDBQTWriterLegacy
from vina import Vina

D = cfg["docking"]
STRUCT = cfg["structure"]
SEED = int(cfg["random_seed"])
N_SCREEN = int(D["n_screen"])
EXHAUST = int(D["exhaustiveness"])
BOX_SIZE = [float(x) for x in D["box_size"]]
N_POSES = int(D["n_poses"])

AA20 = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER "
           "THR TRP TYR VAL".split())


# ---------------------------------------------------------------- 工具
def pdb_coords(fn, hetatm=True):
    """读取 ATOM/HETATM 的坐标与元素"""
    xyz, el = [], []
    for l in open(fn):
        if l.startswith(("ATOM", "HETATM")):
            xyz.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
            el.append(l[76:78].strip().upper())
    return np.array(xyz), el


def rmsd_by_order(xyz_a, xyz_b):
    """两数组原子顺序一致时的 RMSD（Vina 输出与输入 pdbqt 顺序一致）"""
    n = min(len(xyz_a), len(xyz_b))
    return float(np.sqrt(((xyz_a[:n] - xyz_b[:n]) ** 2).sum(axis=1).mean()))


def mol_to_pdbqt(mol, out_file, seed=SEED):
    """RDKit mol（需已带 3D 坐标 + H）→ meeko pdbqt"""
    prep = MoleculePreparation()
    setups = prep.prepare(mol)
    pdbqt = PDBQTWriterLegacy.write_string(setups[0])
    pdbqt = pdbqt[0] if isinstance(pdbqt, tuple) else pdbqt
    with open(out_file, "w") as f:
        f.write(pdbqt)
    return out_file


def embed3d(smiles, seed=SEED):
    """SMILES → 带 H 的 3D 单体分子；失败返回 None"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    ps = AllChem.ETKDGv3()
    ps.randomSeed = seed
    if AllChem.EmbedMolecule(mol, ps) != 0:
        # ETKDG 失败时退回随机坐标，仍尽量给到一次机会
        ps.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, ps) != 0:
            return None
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        pass
    return mol


# ------------------------------------------------- 1. 受体：必须先删掉配体
# 常见翻车点：受体里留着共晶配体 → 配体-配体自碰撞 → 晶体 pose 自打分 +438 kcal/mol。
receptor_pdbqt = interim("receptor.pdbqt")
if exists(receptor_pdbqt):
    log(f"复用已有受体 {os.path.basename(receptor_pdbqt)}")
else:
    log("准备受体（剔除共晶配体 AQ4 与水）…")
    templates = ResidueChemTemplates.create_from_defaults()
    pdb_text = open(raw(f"{STRUCT['pdb_id']}.pdb"), encoding="utf-8", errors="ignore").read()
    # 收集需要删除的残基：所有非标准氨基酸 + 水
    to_del = set()
    for line in pdb_text.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            resn = line[17:20].strip()
            if resn not in AA20:
                to_del.add(f"{line[21] if len(line) > 21 else 'A'}:{int(line[22:26])}")
    to_del = sorted(to_del)
    log(f"将删除 {len(to_del)} 个非标准残基：{to_del[:6]}{' …' if len(to_del) > 6 else ''}")
    polymer = Polymer.from_pdb_string(pdb_text, templates,
                                      mk_prep=MoleculePreparation(),
                                      residues_to_delete=to_del,
                                      default_altloc=str(D["default_altloc"]),
                                      allow_bad_res=bool(D["allow_bad_res"]))
    rec_txt = PDBQTWriterLegacy.write_string_from_polymer(polymer)
    rec_txt = rec_txt[0] if isinstance(rec_txt, tuple) else rec_txt
    with open(receptor_pdbqt, "w") as f:
        f.write(rec_txt)
    n_atom = len([l for l in rec_txt.splitlines() if l.startswith("ATOM")])
    log(f"受体就绪：{n_atom} 个原子，含配体={'AQ4' in rec_txt}")


# ------------------------------------------------- 2. 盒子：以共晶配体为中心
cry_xyz, _ = pdb_coords(raw(f"{STRUCT['pdb_id']}.pdb"))
# 只取配体部分的坐标
lig_xyz = []
for line in open(raw(f"{STRUCT['pdb_id']}.pdb")):
    if line.startswith("HETATM") and line[17:20].strip() == STRUCT["ligand_resname"]:
        lig_xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
lig_xyz = np.array(lig_xyz)
center = lig_xyz.mean(axis=0).round(3)
log(f"对接盒子：center={list(center)} size={BOX_SIZE}（共晶配体 {len(lig_xyz)} 原子的质心）")
save_json({"center": [float(x) for x in center], "size": BOX_SIZE,
           "n_ligand_atoms": int(len(lig_xyz))}, interim("docking_box.json"))


# ------------------------------------------------- 3. 自证：erlotinib 重对接
banner("自检：把共晶配体重新对接回去，看能不能找回晶体 pose")
# 从晶体坐标恢复配体（补键序），而不是用 RDKit 从 SMILES 生成的构象——
# 后者坐标原点附近，不在盒子里，Vina 会直接报 "ligand is outside the grid box"。
ERLOTINIB = "C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1"  # erlotinib，PubChem CID 176870
redock = {}
try:
    tpl = Chem.MolFromSmiles(ERLOTINIB)
    block = "".join(l for l in open(raw(f"{STRUCT['pdb_id']}.pdb"))
                    if l.startswith("HETATM") and l[17:20].strip() == STRUCT["ligand_resname"])
    pmol = Chem.MolFromPDBBlock(block, sanitize=False, removeHs=False)
    ref_mol = AllChem.AssignBondOrdersFromTemplate(tpl, pmol)
    Chem.SanitizeMol(ref_mol)
    ref_mol = Chem.RemoveHs(ref_mol)
    Chem.MolToPDBFile(ref_mol, interim("aq4_crystal.pdb"))
    ref_mol_h = Chem.AddHs(ref_mol, addCoords=True)
    mol_to_pdbqt(ref_mol_h, interim("aq4_crystal.pdbqt"))
    log("晶体配体 erlotinib 恢复成功（含完整晶体坐标）")

    v = Vina(sf_name="vina", cpu=int(D["cpu"]), seed=SEED)
    v.set_receptor(receptor_pdbqt)
    v.set_ligand_from_file(interim("aq4_crystal.pdbqt"))
    v.compute_vina_maps(center=list(center), box_size=BOX_SIZE)
    score_crystal = float(v.score()[0])
    v.dock(exhaustiveness=EXHAUST, n_poses=N_POSES)
    v.write_poses(interim("aq4_redock.pdbqt"), n_poses=N_POSES, overwrite=True)
    scores = [float(e[0]) for e in v.energies(n_poses=N_POSES)]

    xyz_ref, _ = pdb_coords(interim("aq4_crystal.pdbqt"))
    best = None
    for i, blk in enumerate(open(interim("aq4_redock.pdbqt")).read().split("ENDMDL")):
        if "ATOM" not in blk:
            continue
        xyz_i = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                          for l in blk.splitlines() if l.startswith(("ATOM", "HETATM"))])
        r = rmsd_by_order(xyz_i, xyz_ref)
        if best is None or r < best[1]:
            best = (i + 1, r)
    redock = {"crystal_pose_self_score": round(score_crystal, 3),
              "docked_best_score": round(scores[0], 3),
              "n_poses": int(N_POSES),
              "best_pose_rank": int(best[0]),
              "best_pose_rmsd": round(float(best[1]), 3),
              "all_scores": [round(s, 3) for s in scores]}
    log(f"晶体 pose 自打分 = {score_crystal:.3f} kcal/mol（应为明显负值；若为正几百，说明受体里混入了配体）")
    log(f"对接最优 pose 打分 = {scores[0]:.3f} kcal/mol")
    log(f"最优 RMSD = {best[1]:.2f} 埃（第 {best[0]} 个 pose）")
    log("判读：RMSD < 2 埃 = pose 复现成功，这条对接 Pipeline 可信；"
        "> 3 埃 = 盒子/受体/质子化状态里有东西需要回头修。")
except Exception as e:
    log(f"重对接自检跳过：{type(e).__name__}: {e}", "WARN")
    redock = {"error": f"{type(e).__name__}: {e}"}


# ------------------------------------------------- 4. 虚拟筛选
banner(f"虚拟筛选：对 QSAR 预测 Top {N_SCREEN} 的分子做对接")
pred = pd.read_csv(processed("predictions_scaffold.csv"))
pred = pred[pred["split"] == "scaffold"]
pred["molecule_id"] = [f"mol{i:04d}" for i in range(len(pred))]
top = pred.sort_values("pred_pIC50", ascending=False).head(N_SCREEN).reset_index(drop=True)
log(f"取自 predictions_scaffold.csv 的 Top {len(top)}，预测 pIC50 范围 "
    f"{top['pred_pIC50'].min():.2f} ~ {top['pred_pIC50'].max():.2f}")

rows, skipped = [], []
lig_dir = interim("ligands")
os.makedirs(lig_dir, exist_ok=True)

for i, r in top.iterrows():
    smi = str(r["smiles"])
    mol_id = str(r.get("molecule_id", f"mol{i:04d}"))
    mol = embed3d(smi)
    if mol is None:
        skipped.append({"molecule_id": mol_id, "reason": "3D 构象生成失败"})
        continue
    pdbqt = os.path.join(lig_dir, f"{mol_id}.pdbqt")
    try:
        mol_to_pdbqt(mol, pdbqt)
    except Exception as e:
        skipped.append({"molecule_id": mol_id, "reason": f"pdbqt 转换失败：{type(e).__name__}"})
        continue
    try:
        vv = Vina(sf_name="vina", cpu=int(D["cpu"]), seed=SEED)
        vv.set_receptor(receptor_pdbqt)
        vv.set_ligand_from_file(pdbqt)
        vv.compute_vina_maps(center=list(center), box_size=BOX_SIZE)
        vv.dock(exhaustiveness=EXHAUST, n_poses=1)
        s = float(vv.energies(n_poses=1)[0][0])
        vv.write_poses(os.path.join(lig_dir, f"{mol_id}_dock.pdbqt"),
                       n_poses=1, overwrite=True)
    except Exception as e:
        skipped.append({"molecule_id": mol_id, "reason": f"对接失败：{type(e).__name__}"})
        continue
    rows.append({"molecule_id": mol_id, "smiles": smi,
                 "pred_pIC50": float(r["pred_pIC50"]),
                 "docking_score": round(s, 3)})
    if (i + 1) % 5 == 0 or i == len(top) - 1:
        log(f"  {i+1}/{len(top)}  最近={s:.2f} kcal/mol")

dock_df = pd.DataFrame(rows)
if dock_df.empty:
    raise RuntimeError("没有任何分子对接成功，检查配体准备或盒子设置")
log(f"对接完成：成功 {len(dock_df)} / {len(top)}，跳过 {len(skipped)}")
dock_df.to_csv(processed("docking.csv"), index=False)


# ------------------------------------------------- 5. 不要把两个分数直接相加
banner("Consensus：2D 活性模型 + 3D 对接，量纲不同，先标准化再合并")
from sklearn.preprocessing import StandardScaler

X = dock_df[["pred_pIC50", "docking_score"]].values
Xs = StandardScaler().fit_transform(X)
# pIC50 越高越好；docking_score 越负越好 → 取负号后同向
w1 = float(D["w_activity"])
w2 = float(D["w_docking"])
dock_df["consensus"] = w1 * Xs[:, 0] - w2 * Xs[:, 1]
dock_df = dock_df.sort_values("consensus", ascending=False).reset_index(drop=True)
dock_df.insert(0, "rank", range(1, len(dock_df) + 1))
dock_df.to_csv(processed("docking.csv"), index=False)

print(f"\n{'rank':<6}{'pred_pIC50':>12}{'docking':>10}{'consensus':>12}")
print("-" * 42)
for _, r in dock_df.head(10).iterrows():
    print(f"{int(r['rank']):<6}{r['pred_pIC50']:>12.2f}{r['docking_score']:>10.2f}{r['consensus']:>12.3f}")
print("\n  注意：consensus 只是一个排序规则，不是「成药概率」。")

# 两个证据的相关性：几乎不相关才说明 3D 提供了正交信息
corr = float(np.corrcoef(dock_df["pred_pIC50"], dock_df["docking_score"])[0, 1])
log(f"pred_pIC50 与 docking_score 的 Pearson r = {corr:.3f}"
    f"（接近 0 = 两条证据正交，合并才有意义）")

save_json({"redock_check": redock,
           "n_screened": int(len(top)), "n_docked": int(len(dock_df)),
           "n_skipped": int(len(skipped)), "skipped": skipped,
           "box": {"center": [float(x) for x in center], "size": BOX_SIZE},
           "weights": {"activity": w1, "docking": w2},
           "pearson_pred_vs_docking": round(corr, 3),
           "top5": dock_df.head(5)[["molecule_id", "pred_pIC50",
                                    "docking_score", "consensus"]].to_dict("records")},
          processed("09_docking.json"))
banner(f"09 done · 对接成功 {len(dock_df)}/{len(top)}，重对接 RMSD="
       f"{redock.get('best_pose_rmsd', 'n/a')} 埃")
