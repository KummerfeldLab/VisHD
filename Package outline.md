# VisHD package structure (draft)

This branch starts from Lee's `uv` branch and adds a first layer for configuration and job execution. It is a structure proposal, not a choice of scientific method.

## Current files (selected)

```text
VisHD/
├── pyproject.toml
├── uv.lock
├── config/
│   └── vishd.example.toml
├── src/VisHD/
│   ├── __init__.py
│   ├── config.py
│   ├── execution/
│   │   ├── __init__.py
│   │   ├── local_backend.py
│   │   └── slurm.py
│   ├── permutation_pipeline.py
│   └── Visualize.py
└── tests/
    ├── test_config.py
    └── test_execution.py
```

These are the relevant files that exist now; Lee's other modules remain in `src/VisHD`.

## Current state

- Lee's pipeline now takes paths and compute settings from a TOML file or Python object. It runs locally by default; SLURM requires the optional `submitit` extra.
- The older experimental modules still use the `local.py` convention. Mac's scripts have not been added to this repository.
- The GPU extra is declared, but the public pipeline does not run Mac's GPU method yet.
- Tests cover configuration, local job execution, and basic imports. They do not run a full sample or check biological results.
- The README still describes the earlier setup and needs updating after the team agrees on the structure.

## Proposed next modules

Lee's pipeline uses BLADE positions and a KD-tree. Mac's separate `cpu-visium.py` and `visium.py` use Space Ranger counts and dense-grid distance transforms. Proposed `readers/` would convert platform data into one input format; proposed `methods/` would run either analysis on that input and return one result format. `results.py` could save the results. These modules do not exist yet.

Do these module boundaries fit how everyone expects to contribute? Should the first shared input representation be a table of `(sample, x, y, gene, count)`, or is another shape better for both Visium HD and CosMx?

Also, what name do we want for the package? Since we're serving both Visium HD and CosMx, we should consider something different from VisHD.
