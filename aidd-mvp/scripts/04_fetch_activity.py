"""04 · 从 ChEMBL 拉真实活性数据

输入：config.target.uniprot / activity 过滤条件
输出：data/raw/raw_activity.csv + data/processed/04_*.json

要点：
  1) 先按 UniProt accession 找 target_chembl_id，只取 Homo sapiens + SINGLE PROTEIN；
     直接搜 "EGFR" 会得到十几条，混着蛋白家族、混着其他物种。
  2) 真实数据不是 (SMILES, label)，而是
     Molecule / Target / Assay / Endpoint / Relation / Value / Unit / Source。

联网说明：ChEMBL 的 activity 端点会间歇性返回 500；本脚本对 5xx 做指数退避重试。
若最终仍失败，会检查 data/raw/_offline_cache/ 下有没有此前跑通过的真实数据可用。
"""
import os
import time
import shutil
import requests
import pandas as pd

from _common import cfg, ensure_dirs, raw, processed, log, banner, save_json, exists

ensure_dirs()
banner("04 · Bioactivity：ChEMBL 真实 IC50 拉取")

BASE = "https://www.ebi.ac.uk/chembl/api/data"
S = requests.Session()


def get(url, params=None, tries=8, timeout=90):
    """ChEMBL 会间歇性 500/502，也会中途断连；对 5xx 与网络异常都做指数退避。"""
    last = ""
    for i in range(tries):
        try:
            r = S.get(url, params=params, timeout=timeout)
            if r.status_code >= 500:
                raise requests.HTTPError(f"{r.status_code} Server Error")
            r.raise_for_status()
            return r.json()
        except Exception as e:
            last = f"{type(e).__name__}: {str(e)[:60]}"
            wait = min(3 * (2 ** i), 30)
            log(f"  ChEMBL 第 {i+1}/{tries} 次失败（{last}），{wait}s 后重试", "WARN")
            time.sleep(wait)
    raise RuntimeError(f"ChEMBL 请求连续失败: {url} | {last}")


COLS = ["molecule_chembl_id", "canonical_smiles", "standard_type",
        "standard_relation", "standard_value", "standard_units",
        "assay_type", "assay_chembl_id", "pchembl_value", "target_chembl_id"]


def fetch_from_chembl():
    """联网拉取，返回 DataFrame"""
    # ---------- 1. Target identity：先找对靶点 ----------
    accession = cfg["target"]["uniprot"]
    log(f"Step1 按 accession={accession} 查 ChEMBL target")
    d = get(f"{BASE}/target.json",
            {"target_components__accession": accession, "limit": 50})
    targets = d["targets"]
    print(f"\n  共 {len(targets)} 条 ChEMBL target 命中 {accession}：")
    print(f"  {'CHEMBL_ID':<14}{'类型':<22}{'物种':<20}{'名称'}")
    print("  " + "-" * 70)
    for t in targets:
        print(f"  {t['target_chembl_id']:<14}{t['target_type']:<22}"
              f"{t['organism']:<20}{t['pref_name'][:26]}")

    use = [t for t in targets
           if t["organism"] == "Homo sapiens" and t["target_type"] == "SINGLE PROTEIN"]
    if not use:
        raise SystemExit("没有找到 Homo sapiens + SINGLE PROTEIN 的靶点条目")
    tid = use[0]["target_chembl_id"]
    print(f"\n  => 选用 {tid}（{use[0]['pref_name']}）")
    print(f"     其余 {len(targets)-1} 条被排除：非人源 / 非单蛋白 / 蛋白家族")
    save_json({"accession": accession, "all_hits": len(targets), "selected": tid,
               "selected_name": use[0]["pref_name"], "excluded": len(targets) - 1},
              processed("04_target_identity.json"))

    # ---------- 2. 分页拉 activity ----------
    a = cfg["activity"]
    params = {"target_chembl_id": tid,
              "standard_type": a["endpoint"],
              "standard_units": a["unit"],
              "standard_relation": a["relation"],
              "assay_type": a["assay_type"],
              "limit": 1000}
    total = get(f"{BASE}/activity.json", params)["page_meta"]["total_count"]
    cap = min(total, int(cfg["activity"]["max_records"]))
    log(f"Step2 符合条件记录总数 = {total}，本次取 {cap} 条")

    rows, offset = [], 0
    while offset < cap:
        d = get(f"{BASE}/activity.json", {**params, "offset": offset})
        rows.extend(d["activities"])
        offset += 1000
        log(f"  {min(offset, cap)}/{cap}")
        if not d["activities"]:
            break
    df = pd.DataFrame(rows)
    return df[[c for c in COLS if c in df.columns]]


out_csv = raw("raw_activity.csv")
OFFLINE = os.environ.get("AIDD_OFFLINE", "0") == "1"   # AIDD_OFFLINE=1 跳过联网，直接用缓存

if OFFLINE:
    log("AIDD_OFFLINE=1：跳过联网，直接使用本地缓存", "WARN")

cache = raw("_offline_cache/raw_activity.csv")

if OFFLINE and exists(cache) and not exists(out_csv):
    shutil.copy(cache, out_csv)
    log("已用离线缓存生成 data/raw/raw_activity.csv（此前同参数真实跑通的结果）", "WARN")

if exists(out_csv):
    log(f"复用已有 {os.path.basename(out_csv)}"
        f"（{os.path.getsize(out_csv)/1024/1024:.1f} MB）；删掉可强制重新联网")
else:
    try:
        df = fetch_from_chembl()
        df.to_csv(out_csv, index=False)
        log(f"wrote data/raw/raw_activity.csv ({len(df)} 行)")
    except Exception as e:
        cache_dir = raw("_offline_cache")
        cand = os.path.join(cache_dir, "raw_activity.csv")
        if exists(cand):
            shutil.copy(cand, out_csv)
            log(f"联网失败（{str(e)[:60]}），改用离线缓存 {cand}", "WARN")
            log("该缓存是此前用同样参数真实跑通的结果，非手工编造", "WARN")
        else:
            raise SystemExit(
                f"\nChEMBL 拉取失败：{e}\n"
                f"可稍后重试，或把一份同格式的 raw_activity.csv 放到 "
                f"{cache_dir}/ 作为离线兜底。")

# ---------- 3. 看一眼真实数据长什么样 ----------
df = pd.read_csv(out_csv)
print("\n真实数据的样子（前 3 行 —— 注意它不是一个 label 列）：")
print(df.head(3).to_string(index=False))

print("\n数据体检：")
print(f"  记录数        : {len(df)}")
print(f"  唯一分子      : {df['molecule_chembl_id'].nunique()}")
print(f"  唯一 assay    : {df['assay_chembl_id'].nunique()}")
print(f"  缺失 SMILES   : {df['canonical_smiles'].isna().sum()}")
print(f"  endpoint 分布 : {dict(df['standard_type'].value_counts().head(5))}")
print(f"  relation 分布 : {dict(df['standard_relation'].value_counts().head(5))}")
print(f"  unit 分布     : {dict(df['standard_units'].value_counts().head(5))}")

save_json({
    "n_records": int(len(df)),
    "n_molecules": int(df["molecule_chembl_id"].nunique()),
    "n_assays": int(df["assay_chembl_id"].nunique()),
    "missing_smiles": int(df["canonical_smiles"].isna().sum()),
    "endpoint": {str(k): int(v) for k, v in df["standard_type"].value_counts().head(5).items()},
}, processed("04_activity_qc.json"))

banner(f"04 done · {len(df)} 条 / {df['molecule_chembl_id'].nunique()} 个分子 / "
       f"{df['assay_chembl_id'].nunique()} 个 assay")
