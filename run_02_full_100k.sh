#!/usr/bin/env bash
set -euo pipefail

python pipeline/02_pepmlm_generate.py \
  --config config/02_pepmlm_generate.json
