# 文章关联：AI制药实践篇B_实操篇_公众号文章.html
# 脚本职责：早期实验 2：ChEMBL EGFR IC50 数据抓取
# 使用注解：输出 raw_activity.csv，带缓存和断点；接早期 step3，不能与正式工程数据目录混用。
# 完整对应关系与运行顺序：README.md
"""EGFR MVP Step2: 从 ChEMBL 拉取 CHEMBL203 的真实 IC50 活性数据（带重试/断点续传）"""
import os, json, time
import requests
import pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(BASE, "data", "egfr_mvp")
os.makedirs(DATA, exist_ok=True)
CACHE = os.path.join(DATA, "raw_activity.csv")
PART = os.path.join(DATA, "raw_activity_partial.json")

URL = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
FIELDS = ["molecule_chembl_id", "canonical_smiles", "standard_type", "standard_relation",
          "standard_value", "standard_units", "assay_type", "assay_chembl_id", "pchembl_value"]
PARAMS = {
    "target_chembl_id": "CHEMBL203",
    "standard_type": "IC50",
    "standard_units": "nM",
    "standard_relation": "=",
    "assay_type": "B",
    "limit": 1000,
}

def fetch_page(offset, tries=6):
    p = dict(PARAMS); p["offset"] = offset; p["only"] = ",".join(FIELDS)
    last = None
    for k in range(tries):
        try:
            s = requests.Session()
            s.trust_env = True
            r = s.get(URL, params=p, timeout=90)
            r.raise_for_status()
            return r.json()
        except Exception as e:
            last = e
            time.sleep(2 + 3 * k)
    raise RuntimeError(f"offset {offset} failed after {tries} tries: {last}")

if os.path.exists(CACHE):
    df = pd.read_csv(CACHE)
    print("cache hit:", CACHE, df.shape)
else:
    rows = []
    if os.path.exists(PART):
        rows = json.load(open(PART))
        print(f"resume from partial: {len(rows)} rows")
    total = None
    t0 = time.time()
    while True:
        j = fetch_page(len(rows))
        total = j["page_meta"]["total_count"]
        acts = j["activities"]
        rows.extend(acts)
        print(f"  {len(rows)}/{total}  ({time.time()-t0:.0f}s)", flush=True)
        if len(rows) % 5000 < 1000:
            json.dump(rows, open(PART, "w"))
        if len(acts) < 1000 or len(rows) >= total:
            break
        time.sleep(0.4)
    json.dump(rows, open(PART, "w"))
    df = pd.DataFrame(rows)[FIELDS]
    df.to_csv(CACHE, index=False)
    print("saved ->", CACHE, df.shape)

print("\n=== raw ===")
print("rows:", len(df))
print("unique molecules:", df["molecule_chembl_id"].nunique())
print("unique assays   :", df["assay_chembl_id"].nunique())
print("null smiles     :", df["canonical_smiles"].isna().sum())
print("\nfirst 3 rows:")
print(df.head(3).to_string())
print("\nstandard_value dtype:", df["standard_value"].dtype)
