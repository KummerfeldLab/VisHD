"""Optional SLURM support through submitit."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from VisHD.config import ExecutionConfig


def make_slurm_executor(config: ExecutionConfig, log_dir: str | Path) -> Any:
    """Create a submitit executor and provide a useful optional-dependency error."""

    if config.slurm is None:
        raise ValueError("SLURM settings are required for the SLURM backend")

    try:
        import submitit
    except ImportError as error:
        raise ImportError(
            "SLURM execution requires the optional dependency: "
            "install VisHD with the 'slurm' extra"
        ) from error

    executor = submitit.AutoExecutor(folder=Path(log_dir))
    parameters: dict[str, Any] = {
        "slurm_account": config.slurm.account,
        "mem_gb": config.slurm.memory_gb,
        "cpus_per_task": config.slurm.cpus,
        "slurm_time": config.slurm.time_limit,
    }
    if config.slurm.partition is not None:
        parameters["slurm_partition"] = config.slurm.partition
    executor.update_parameters(**parameters)
    return executor
