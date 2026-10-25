"""aidd-mvp 公共工具：配置加载、路径、依赖引导、统一日志。

每个脚本开头都做三件事：
    from _common import ROOT, cfg, log, ensure_dirs
这样所有脚本都共享同一套 config.yaml 与目录约定，
任何一个脚本都能被单独 python scripts/0X_xxx.py 执行。
"""
import os
import sys
import time
import datetime

# ---- 依赖引导：biopython / meeko 装在工程内 vendor/，避免污染 conda 环境 ----
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VENDOR = os.path.join(ROOT, "vendor")
if os.path.isdir(VENDOR) and VENDOR not in sys.path:
    sys.path.insert(0, VENDOR)


def load_cfg():
    import yaml
    with open(os.path.join(ROOT, "config.yaml"), encoding="utf-8") as f:
        return yaml.safe_load(f)


cfg = load_cfg()
P = cfg["paths"]


def path(*parts):
    """相对于工程根目录取路径"""
    return os.path.join(ROOT, *parts)


def raw(name):
    return path(P["raw"], name)


def interim(name):
    return path(P["interim"], name)


def processed(name):
    return path(P["processed"], name)


def ensure_dirs():
    for k in ("raw", "interim", "processed", "models", "reports", "figures"):
        os.makedirs(path(P[k]), exist_ok=True)


_T0 = time.time()


def log(msg, level="INFO"):
    dt = time.time() - _T0
    print(f"[{level:<5}] {dt:6.1f}s | {msg}", flush=True)


def banner(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72, flush=True)


def exists(p):
    return os.path.exists(p) and os.path.getsize(p) > 0


def save_json(obj, p):
    import json
    os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)
    with open(p, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=2)
    log(f"wrote {os.path.relpath(p, ROOT)}")


def stamp():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
