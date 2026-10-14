#!/usr/bin/env bash
set -euo pipefail

python pipeline/03_sequence_filter.py \
  --config config/03_sequence_filter.json
