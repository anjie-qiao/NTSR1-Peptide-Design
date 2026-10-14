"""Lazy torch/Transformers loading helpers shared by GPU stages."""

from __future__ import annotations

from pathlib import Path

from pepmlm_pipeline_common import CANONICAL_AAS


def select_device_and_dtype(torch, requested_device: str, requested_dtype: str):
    if requested_device == "auto":
        device = "cuda" if torch.cuda.is_available() else "cpu"
    else:
        device = requested_device
    if device.startswith("cuda") and not torch.cuda.is_available():
        raise RuntimeError("CUDA was requested but torch.cuda.is_available() is False")

    if requested_dtype == "auto":
        if device.startswith("cuda"):
            dtype = torch.bfloat16 if torch.cuda.is_bf16_supported() else torch.float16
        else:
            dtype = torch.float32
    else:
        if not hasattr(torch, requested_dtype):
            raise ValueError(f"Unknown torch dtype: {requested_dtype}")
        dtype = getattr(torch, requested_dtype)
    return device, dtype


def load_local_model(model_path: Path, requested_device: str = "auto", requested_dtype: str = "auto"):
    try:
        import torch
        from transformers import AutoModelForMaskedLM, AutoTokenizer
    except ImportError as exc:
        raise RuntimeError("This stage requires torch and transformers") from exc

    if not model_path.is_dir():
        raise FileNotFoundError(f"Local PepMLM directory not found: {model_path}")
    device, dtype = select_device_and_dtype(torch, requested_device, requested_dtype)
    tokenizer = AutoTokenizer.from_pretrained(str(model_path), local_files_only=True)
    kwargs = {"local_files_only": True, "low_cpu_mem_usage": True}
    try:
        model = AutoModelForMaskedLM.from_pretrained(str(model_path), dtype=dtype, **kwargs)
    except TypeError:
        model = AutoModelForMaskedLM.from_pretrained(str(model_path), torch_dtype=dtype, **kwargs)
    model = model.to(device)
    model.eval()
    return torch, tokenizer, model, device, dtype


def canonical_token_ids(tokenizer) -> list[int]:
    token_ids = [int(tokenizer.convert_tokens_to_ids(aa)) for aa in CANONICAL_AAS]
    if len(set(token_ids)) != len(CANONICAL_AAS):
        raise ValueError("Tokenizer does not map the 20 canonical amino acids to unique tokens")
    if tokenizer.unk_token_id in token_ids:
        raise ValueError("At least one canonical amino acid mapped to the unknown token")
    return token_ids
