#!/usr/bin/env bash
set -euo pipefail

python pipeline/09_audit_optimized_structures.py \
  --config config/09_audit_optimized_structures.json
