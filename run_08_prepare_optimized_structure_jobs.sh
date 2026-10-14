#!/usr/bin/env bash
set -euo pipefail

python pipeline/08_prepare_optimized_structure_jobs.py \
  --config config/08_prepare_optimized_structure_jobs.json
