"""00 · 环境预检：跑之前先回答"这台机器能不能跑"。

用法：python scripts/00_preflight.py
退出码：0 = 可以跑；1 = 缺少必需依赖，先装再跑。

设计原则：
  * 必需 / 可选分开。缺必需 -> FAIL；缺可选 -> WARN（对应脚本会自动跳过）。
  * 只检查，不写任何产物，不联网（联网检查放在 02/04 各自脚本里）。
"""
import os
import sys
import shutil
import importlib

from _common import ROOT, cfg, path, P, log, banner

# name -> (import 名, 是否必需, 用途)
DEPS = [
    ("numpy", "numpy", True, "矩阵与数值"),
    ("pandas", "pandas", True, "表格读写"),
    ("yaml", "yaml", True, "config.yaml 解析"),
    ("rdkit", "rdkit", True, "分子标准化 / 指纹 / 描述符"),
    ("sklearn", "sklearn", True, "RandomForest + GroupShuffleSplit"),
    ("requests", "requests", True, "RCSB / UniProt / ChEMBL 拉取"),
    ("joblib", "joblib", True, "模型落盘"),
    ("Bio", "biopython", False, "03 口袋提取（PDB 解析 / 邻居搜索）"),
    ("meeko", "meeko", False, "09 受体制备与配体 pdbqt 写出"),
    ("gemmi", "gemmi", False, "09 结构读写兜底"),
    ("torch", "torch", False, "08 chemprop D-MPNN"),
    ("chemprop", "chemprop", False, "08 chemprop D-MPNN"),
    ("matplotlib", "matplotlib", False, "出图（非流水线必需）"),
]

REQUIRED_CFG_KEYS = [
    "random_seed",
    "target.uniprot", "target.chembl_id",
    "structure.pdb_id", "structure.ligand_resname",
    "structure.pocket_radius", "structure.numbering_offset_to_uniprot",
    "model.fingerprint_radius", "model.fingerprint_size", "model.test_size",
    "paths.raw", "paths.interim", "paths.processed",
]


def version_of(mod):
    try:
        m = importlib.import_module(mod)
        return getattr(m, "__version__", "?")
    except Exception:
        return None


def get_cfg(dotted):
    cur = cfg
    for k in dotted.split("."):
        if not isinstance(cur, dict) or k not in cur:
            return None
        cur = cur[k]
    return cur


def main():
    banner("00 · Preflight：环境与配置预检")
    fails, warns = [], []

    print(f"Python  {sys.version.split()[0]}   ({sys.executable})")
    print(f"工程根  {ROOT}")
    print(f"vendor  {'存在（依赖隔离）' if os.path.isdir(os.path.join(ROOT, 'vendor')) else '不存在（用系统环境）'}")
    print()

    print(f"{'依赖':<12}{'版本':<16}{'必需':<8}{'用途'}")
    print("-" * 66)
    for imp, mod, required, usage in DEPS:
        v = version_of(imp)
        if v is None:
            mark = "FAIL" if required else "WARN"
            (fails if required else warns).append(mod)
            print(f"{mod:<12}{'-':<16}{mark:<8}{usage}")
        else:
            print(f"{mod:<12}{str(v)[:14]:<16}{'OK':<8}{usage}")

    print()
    print("外部可执行文件：")
    vina = shutil.which("vina")
    if not vina:  # PATH 里没有时，再找工程 vendor/ 与当前解释器的 bin/
        for cand in (os.path.join(ROOT, "vendor", "vina"),
                     os.path.join(os.path.dirname(sys.executable), "vina"),
                     os.path.join(sys.prefix, "bin", "vina")):
            if os.path.exists(cand) and os.access(cand, os.X_OK):
                vina = cand
                break
    if vina:
        print(f"  vina            OK   {vina}")
    else:
        warns.append("vina")
        print("  vina            WARN 未找到（09_docking 会自动跳过对接，其余步骤不受影响）")

    print()
    print("config.yaml 必填键：")
    for k in REQUIRED_CFG_KEYS:
        v = get_cfg(k)
        if v is None:
            fails.append(f"config.{k}")
            print(f"  {k:<42} FAIL 缺失")
        else:
            print(f"  {k:<42} OK   {v}")

    print()
    print("目录可写性：")
    for k in ("raw", "interim", "processed", "models", "reports"):
        d = path(P[k])
        try:
            os.makedirs(d, exist_ok=True)
            t = os.path.join(d, ".preflight.tmp")
            with open(t, "w") as f:
                f.write("x")
            os.remove(t)
            print(f"  {P[k]:<16} OK   可写")
        except Exception as e:
            fails.append(f"dir:{P[k]}")
            print(f"  {P[k]:<16} FAIL {e}")

    print()
    if fails:
        print(f"FAIL：{len(fails)} 项必需依赖/配置不满足 -> {fails}")
        print("先补齐再跑：pip install <缺失包>，或检查 config.yaml。")
        return 1
    print("PASS：必需依赖与配置齐全" + (f"（{len(warns)} 项可选缺失：{warns}）" if warns else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
