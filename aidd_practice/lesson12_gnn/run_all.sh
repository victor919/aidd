#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
PYTHON_BIN="${PYTHON_BIN:-python}"
"$PYTHON_BIN" scripts/e18_prepare.py
"$PYTHON_BIN" scripts/e19_benchmark.py
"$PYTHON_BIN" scripts/e20_sensitivity.py
"$PYTHON_BIN" scripts/e21_figures.py
"$PYTHON_BIN" scripts/e23_attachment_audit.py
"$PYTHON_BIN" scripts/e24_supplement_figures.py
"$PYTHON_BIN" scripts/e25_case_readout.py
"$PYTHON_BIN" scripts/make_course_map.py
"$PYTHON_BIN" scripts/e22_article.py
"$PYTHON_BIN" scripts/verify_article.py
