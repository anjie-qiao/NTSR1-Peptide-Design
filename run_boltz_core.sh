#!/usr/bin/env bash
set -euo pipefail

bash work/05_prepare_structure_jobs/run_boltz_preflight.sh

boltz predict work/05_prepare_structure_jobs/boltz_inputs_core \
  --out_dir work/05_prepare_structure_jobs/boltz_predictions_core \
  --accelerator gpu \
  --devices 1 \
  --recycling_steps 3 \
  --sampling_steps 200 \
  --diffusion_samples 1 \
  --output_format pdb \
  --write_full_pae \
  --use_potentials \
  --no_kernels \
  --max_msa_seqs 512
