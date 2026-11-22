#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
PYTHON_BIN="${PYTHON_BIN:-python}"
"$PYTHON_BIN" scripts/experiment.py
"$PYTHON_BIN" scripts/analyze.py
"$PYTHON_BIN" scripts/attachment_readout.py
"$PYTHON_BIN" scripts/figures.py
"$PYTHON_BIN" scripts/article.py
"$PYTHON_BIN" scripts/verify.py
