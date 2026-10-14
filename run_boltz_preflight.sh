#!/usr/bin/env bash
set -euo pipefail

command -v boltz >/dev/null || { echo 'ERROR: boltz is not installed in this environment' >&2; exit 1; }
python - <<'PY'
import torch
print('torch:', torch.__version__)
print('torch CUDA:', torch.version.cuda)
print('CUDA available:', torch.cuda.is_available())
if not torch.cuda.is_available():
    raise SystemExit('ERROR: CUDA is unavailable; refusing to run Boltz on CPU')
print('GPU:', torch.cuda.get_device_name(0))
PY
boltz predict --help >/dev/null
echo 'Boltz/GPU preflight PASS'
