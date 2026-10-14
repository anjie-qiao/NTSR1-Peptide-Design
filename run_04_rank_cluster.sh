#!/usr/bin/env bash
set -euo pipefail

python pipeline/04_rank_cluster.py \
  --config config/04_rank_cluster.json
