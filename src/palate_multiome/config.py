"""Config loading. Every threshold in the pipeline comes from here."""
from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

DEFAULT_CONFIG = Path(__file__).resolve().parents[2] / "config" / "config.yaml"


def load_config(path: str | Path | None = None) -> dict[str, Any]:
    """Read the YAML config.

    Parameters
    ----------
    path
        Config file. Defaults to ``config/config.yaml`` at the repo root.
    """
    path = Path(path) if path else DEFAULT_CONFIG
    if not path.exists():
        raise FileNotFoundError(f"config not found: {path}")
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    _validate(cfg)
    return cfg


def _validate(cfg: dict[str, Any]) -> None:
    required = ["seed", "data", "samples", "qc", "rna", "markers", "priming"]
    missing = [k for k in required if k not in cfg]
    if missing:
        raise KeyError(f"config missing required sections: {missing}")
    stages = {s["stage"] for s in cfg["samples"]}
    unknown = stages - set(cfg["stage_order"])
    if unknown:
        raise ValueError(f"samples reference stages absent from stage_order: {unknown}")


def sample_ids(cfg: dict[str, Any]) -> list[str]:
    return [s["id"] for s in cfg["samples"]]


def stage_map(cfg: dict[str, Any]) -> dict[str, str]:
    return {s["id"]: s["stage"] for s in cfg["samples"]}
