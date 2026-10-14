"""Record installed versions and validate a GPU, without loading model weights."""
import argparse
import importlib.metadata as metadata
import json
import platform


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--backend", choices=("torch", "jax"), required=True)
    args = parser.parse_args()
    packages = {}
    for name in ("torch", "transformers", "jax", "jaxlib", "colabdesign", "boltz",
                 "numpy", "matplotlib", "chex", "optax", "dm-haiku", "wandb"):
        try:
            packages[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            continue
    data = {"python": platform.python_version(), "platform": platform.platform(),
            "libc": platform.libc_ver(), "packages": packages}
    if args.backend == "torch":
        import torch
        data["torch_cuda"] = torch.version.cuda
        if not torch.cuda.is_available():
            raise RuntimeError("CUDA GPU unavailable in this environment")
        data["gpus"] = [{"name": torch.cuda.get_device_name(i), "vram_bytes": torch.cuda.get_device_properties(i).total_memory}
                        for i in range(torch.cuda.device_count())]
    else:
        import jax
        import colabdesign
        data["gpus"] = [str(device) for device in jax.devices() if device.platform == "gpu"]
        if not data["gpus"]:
            raise RuntimeError("JAX GPU unavailable")
        data["jax_tree_map_available"] = hasattr(jax, "tree_map")
        if not data["jax_tree_map_available"]:
            raise RuntimeError("ColabDesign v1.1.3 requires legacy jax.tree_map; use the pinned baseline")
    print(json.dumps(data, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
