#!/usr/bin/env bash
# aidd-mvp 一键跑通（含前置环境预检 + 后置自检）
#
# 用法：
#   bash run_all.sh                                # 联网拉 ChEMBL
#   AIDD_OFFLINE=1 bash run_all.sh                 # 用 data/raw/_offline_cache 缓存
#   PYTHON=/path/to/python bash run_all.sh         # 指定解释器（推荐 conda env aidd）
#
# 三个阶段：
#   1) 00_preflight.py   环境与配置预检，缺必需依赖直接退出（退出码 1）
#   2) 01~11             流水线主体（按文件名字典序，即依赖顺序）
#   3) make_report / check_artifacts / 99_selftest 出报告 + 产物自洽 + 不变量断言
set -euo pipefail
cd "$(dirname "$0")"

PY="${PYTHON:-python}"

echo "################################################################"
echo "# 0/3  Preflight：环境与配置预检"
echo "################################################################"
"$PY" scripts/00_preflight.py

echo ""
echo "################################################################"
echo "# 1/3  Pipeline：01 -> 11"
echo "################################################################"
# 只跑 01~11：00 是预检（已在上面跑过），99 是自检（放在最后）
for s in scripts/[01][0-9]_*.py; do
  echo ""
  echo "----------------------------------------------------------------"
  echo "# $s"
  echo "----------------------------------------------------------------"
  case "$s" in
    scripts/04_*) AIDD_OFFLINE="${AIDD_OFFLINE:-1}" "$PY" "$s" ;;   # 默认走缓存，省 5 分钟
    *)            "$PY" "$s" ;;
  esac
done

echo ""
echo "################################################################"
echo "# 2/3  Report：汇总报告"
echo "################################################################"
"$PY" scripts/make_report.py

echo ""
echo "################################################################"
echo "# 3/3  Verify：产物自洽 + 不变量自检（任一 FAIL 即非零退出）"
echo "################################################################"
"$PY" scripts/check_artifacts.py
"$PY" scripts/99_selftest.py

echo ""
echo "================================================================"
echo "MVP 跑通并通过全部自检。关键产物："
echo "  data/processed/target_rank.csv        靶点排序"
echo "  data/processed/pocket.json            5 埃口袋 + 编号偏移校验"
echo "  data/processed/dataset.csv            清洗后的活性数据集"
echo "  data/processed/07_baseline.json       ECFP4+RF baseline"
echo "  data/processed/08_chemprop.json       D-MPNN 对照"
echo "  data/processed/docking.csv            Vina 对接 + consensus"
echo "  data/processed/candidate_ranking.csv  四组视角候选清单"
echo "  data/processed/iteration_history.csv  主动学习轮次历史"
echo "  reports/summary.md / validation.md    汇总与验证报告"
echo "================================================================"
