"""User-facing configuration for VisHD.

Configuration is explicit: callers either construct these dataclasses in Python or
pass a TOML file to :func:`load_config`. VisHD does not search a home directory for
configuration and does not store machine-specific paths in package modules.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import tomllib
from typing import Any, Literal, Mapping


class ConfigError(ValueError):
    """Raised when a VisHD configuration is incomplete or invalid."""


@dataclass(frozen=True, slots=True)
class PathConfig:
    """Locations of user data and generated results."""

    spaceranger_root: Path
    output_root: Path
    positions_root: Path | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "spaceranger_root", Path(self.spaceranger_root))
        object.__setattr__(self, "output_root", Path(self.output_root))
        if self.positions_root is not None:
            object.__setattr__(self, "positions_root", Path(self.positions_root))


@dataclass(frozen=True, slots=True)
class SlurmConfig:
    """Settings used only when ``execution.backend`` is ``"slurm"``."""

    account: str
    memory_gb: int = 64
    cpus: int = 4
    time_limit: str = "06:00:00"
    partition: str | None = None

    def __post_init__(self) -> None:
        if not self.account.strip():
            raise ConfigError("execution.slurm.account cannot be empty")
        if self.memory_gb < 1:
            raise ConfigError("execution.slurm.memory_gb must be at least 1")
        if self.cpus < 1:
            raise ConfigError("execution.slurm.cpus must be at least 1")


@dataclass(frozen=True, slots=True)
class ExecutionConfig:
    """How VisHD should run work."""

    backend: Literal["local", "slurm"] = "local"
    workers: int = 1
    use_gpu: bool = False
    slurm: SlurmConfig | None = None

    def __post_init__(self) -> None:
        if self.backend not in {"local", "slurm"}:
            raise ConfigError("execution.backend must be 'local' or 'slurm'")
        if self.workers < 1:
            raise ConfigError("execution.workers must be at least 1")
        if self.backend == "slurm" and self.slurm is None:
            raise ConfigError(
                "execution.slurm is required when execution.backend is 'slurm'"
            )


@dataclass(frozen=True, slots=True)
class VisHDConfig:
    """Complete configuration for one VisHD run."""

    paths: PathConfig
    execution: ExecutionConfig = field(default_factory=ExecutionConfig)


def _only_keys(section: str, values: Mapping[str, Any], allowed: set[str]) -> None:
    unknown = set(values) - allowed
    if unknown:
        names = ", ".join(sorted(unknown))
        raise ConfigError(f"unknown key(s) in [{section}]: {names}")


def _resolve_path(value: Any, *, name: str, base: Path) -> Path:
    if not isinstance(value, (str, Path)):
        raise ConfigError(f"{name} must be a path string")
    path = Path(value).expanduser()
    return path if path.is_absolute() else (base / path).resolve()


def load_config(path: str | Path) -> VisHDConfig:
    """Load a VisHD TOML file.

    Relative paths are interpreted relative to the configuration file, which makes
    a project configuration portable between machines and working directories.
    """

    config_path = Path(path).expanduser().resolve()
    with config_path.open("rb") as stream:
        raw = tomllib.load(stream)

    _only_keys("root", raw, {"paths", "execution"})
    paths_raw = raw.get("paths")
    if not isinstance(paths_raw, dict):
        raise ConfigError("a [paths] section is required")
    _only_keys(
        "paths",
        paths_raw,
        {"spaceranger_root", "output_root", "positions_root"},
    )

    missing = {"spaceranger_root", "output_root"} - set(paths_raw)
    if missing:
        names = ", ".join(sorted(missing))
        raise ConfigError(f"missing required path setting(s): {names}")

    base = config_path.parent
    positions_value = paths_raw.get("positions_root")
    paths = PathConfig(
        spaceranger_root=_resolve_path(
            paths_raw["spaceranger_root"], name="paths.spaceranger_root", base=base
        ),
        output_root=_resolve_path(
            paths_raw["output_root"], name="paths.output_root", base=base
        ),
        positions_root=(
            _resolve_path(
                positions_value, name="paths.positions_root", base=base
            )
            if positions_value is not None
            else None
        ),
    )

    execution_raw = raw.get("execution", {})
    if not isinstance(execution_raw, dict):
        raise ConfigError("[execution] must be a table")
    _only_keys("execution", execution_raw, {"backend", "workers", "use_gpu", "slurm"})

    slurm_raw = execution_raw.get("slurm")
    slurm = None
    if slurm_raw is not None:
        if not isinstance(slurm_raw, dict):
            raise ConfigError("[execution.slurm] must be a table")
        _only_keys(
            "execution.slurm",
            slurm_raw,
            {"account", "memory_gb", "cpus", "time_limit", "partition"},
        )
        try:
            slurm = SlurmConfig(**slurm_raw)
        except TypeError as error:
            raise ConfigError(f"invalid [execution.slurm] settings: {error}") from error

    try:
        execution = ExecutionConfig(
            backend=execution_raw.get("backend", "local"),
            workers=execution_raw.get("workers", 1),
            use_gpu=execution_raw.get("use_gpu", False),
            slurm=slurm,
        )
    except TypeError as error:
        raise ConfigError(f"invalid [execution] settings: {error}") from error

    return VisHDConfig(paths=paths, execution=execution)
