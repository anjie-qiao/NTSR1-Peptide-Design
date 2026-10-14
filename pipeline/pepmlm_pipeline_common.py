"""Shared dependency-light helpers for the PepMLM pipeline."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import tempfile
from pathlib import Path
from typing import Iterable, Sequence


CANONICAL_AAS = "ACDEFGHIKLMNPQRSTVWY"
HYDROPHOBIC_AAS = frozenset("AVILMFWY")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_json_hash(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return sha256_bytes(payload.encode("utf-8"))


def read_single_fasta(path: Path) -> tuple[str, str]:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not lines or not lines[0].startswith(">"):
        raise ValueError(f"Expected a single-record FASTA: {path}")
    if any(line.startswith(">") for line in lines[1:]):
        raise ValueError(f"Multiple FASTA records are not supported: {path}")
    sequence = "".join(lines[1:]).upper()
    invalid = sorted(set(sequence) - set(CANONICAL_AAS))
    if not sequence or invalid:
        raise ValueError(f"Invalid target sequence in {path}; invalid symbols={invalid}")
    return lines[0][1:], sequence


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(text)
        os.replace(temp_name, path)
    except Exception:
        try:
            os.unlink(temp_name)
        except FileNotFoundError:
            pass
        raise


def atomic_write_json(path: Path, value: object) -> None:
    atomic_write_text(path, json.dumps(value, indent=2, ensure_ascii=False) + "\n")


def atomic_write_csv(path: Path, rows: Iterable[dict], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temp_name, path)
    except Exception:
        try:
            os.unlink(temp_name)
        except FileNotFoundError:
            pass
        raise


def allocate_integer(total: int, weights: Sequence[float]) -> list[int]:
    if total < 0:
        raise ValueError("total must be non-negative")
    if not weights or any(weight < 0 for weight in weights) or sum(weights) <= 0:
        raise ValueError("weights must contain at least one positive value and no negatives")
    normalized = [weight / sum(weights) for weight in weights]
    raw = [total * weight for weight in normalized]
    counts = [math.floor(value) for value in raw]
    remaining = total - sum(counts)
    order = sorted(range(len(weights)), key=lambda idx: (raw[idx] - counts[idx], -idx), reverse=True)
    for idx in order[:remaining]:
        counts[idx] += 1
    return counts


def validate_sampling_profiles(profiles: Sequence[dict]) -> None:
    if not profiles:
        raise ValueError("At least one sampling profile is required")
    names = set()
    for profile in profiles:
        required = {"name", "weight", "top_k", "temperature"}
        missing = required - set(profile)
        if missing:
            raise ValueError(f"Sampling profile is missing keys: {sorted(missing)}")
        if profile["name"] in names:
            raise ValueError(f"Duplicate sampling profile name: {profile['name']}")
        names.add(profile["name"])
        if not (1 <= int(profile["top_k"]) <= len(CANONICAL_AAS)):
            raise ValueError(f"top_k must be between 1 and {len(CANONICAL_AAS)}")
        if float(profile["temperature"]) <= 0:
            raise ValueError("temperature must be positive")
        if float(profile["weight"]) < 0:
            raise ValueError("profile weight cannot be negative")


def longest_run(sequence: str, alphabet: set[str] | frozenset[str] | None = None) -> int:
    best = current = 0
    previous = None
    for aa in sequence:
        if alphabet is None:
            current = current + 1 if aa == previous else 1
        else:
            current = current + 1 if aa in alphabet else 0
        best = max(best, current)
        previous = aa
    return best


def normalized_sequence_entropy(sequence: str) -> float:
    counts = {aa: sequence.count(aa) for aa in set(sequence)}
    entropy = -sum((count / len(sequence)) * math.log(count / len(sequence)) for count in counts.values())
    return entropy / math.log(min(len(sequence), len(CANONICAL_AAS))) if len(sequence) > 1 else 0.0


def approximate_net_charge(sequence: str, ph: float = 7.4) -> float:
    # Henderson-Hasselbalch approximation for free peptide termini and ionizable side chains.
    positive = 1 / (1 + 10 ** (ph - 8.0))
    positive += sequence.count("K") / (1 + 10 ** (ph - 10.5))
    positive += sequence.count("R") / (1 + 10 ** (ph - 12.5))
    positive += sequence.count("H") / (1 + 10 ** (ph - 6.0))
    negative = 1 / (1 + 10 ** (3.1 - ph))
    negative += sequence.count("D") / (1 + 10 ** (3.9 - ph))
    negative += sequence.count("E") / (1 + 10 ** (4.1 - ph))
    negative += sequence.count("C") / (1 + 10 ** (8.3 - ph))
    negative += sequence.count("Y") / (1 + 10 ** (10.1 - ph))
    return positive - negative


def sequence_properties(sequence: str) -> dict[str, float | int]:
    return {
        "net_charge_ph7_4": approximate_net_charge(sequence),
        "hydrophobic_fraction": sum(aa in HYDROPHOBIC_AAS for aa in sequence) / len(sequence),
        "max_hydrophobic_run": longest_run(sequence, HYDROPHOBIC_AAS),
        "max_homopolymer_run": longest_run(sequence),
        "normalized_entropy": normalized_sequence_entropy(sequence),
        "cysteine_count": sequence.count("C"),
    }
