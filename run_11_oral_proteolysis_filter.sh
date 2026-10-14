#!/usr/bin/env bash
set -euo pipefail

python pipeline/11_oral_proteolysis_filter.py \
  --config config/11_oral_proteolysis_filter.json
