#!/usr/bin/env bash
set -euo pipefail

python pipeline/05_prepare_structure_jobs.py \
  --config config/05_prepare_structure_jobs.json
