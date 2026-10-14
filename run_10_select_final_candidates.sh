#!/usr/bin/env bash
set -euo pipefail

python pipeline/10_select_final_candidates.py \
  --config config/10_select_final_candidates.json
