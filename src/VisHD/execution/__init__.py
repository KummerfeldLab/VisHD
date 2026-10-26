"""Execution backends for VisHD."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from VisHD.config import ExecutionConfig
from VisHD.execution.local_backend import LocalExecutor


def make_executor(config: ExecutionConfig, log_dir: str | Path) -> Any:
    """Create the selected executor without importing optional tools unnecessarily."""

    if config.backend == "local":
        return LocalExecutor(workers=config.workers)

    from VisHD.execution.slurm import make_slurm_executor

    return make_slurm_executor(config, log_dir)


__all__ = ["LocalExecutor", "make_executor"]
