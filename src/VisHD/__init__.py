"""VisHD: spatial hotspot analysis for Visium HD and related platforms."""

from VisHD.config import (
    ConfigError,
    ExecutionConfig,
    PathConfig,
    SlurmConfig,
    VisHDConfig,
    load_config,
)


def run_pipeline(*args, **kwargs):
    """Run the permutation pipeline without importing heavy dependencies up front."""

    from VisHD.permutation_pipeline import run_pipeline as _run_pipeline

    return _run_pipeline(*args, **kwargs)


__all__ = [
    "ConfigError",
    "ExecutionConfig",
    "PathConfig",
    "SlurmConfig",
    "VisHDConfig",
    "load_config",
    "run_pipeline",
]
