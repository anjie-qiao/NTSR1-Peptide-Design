#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_DIR"
export PYTHONUTF8=1 PYTHONUNBUFFERED=1 WANDB_MODE=disabled
exec "${PYTHON_BIN:-python}" -u pipeline/cached_runner.py "$@"
