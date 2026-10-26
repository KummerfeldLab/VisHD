# VisHD package structure (draft)

This branch starts from Lee's `uv` branch and adds a first layer for configuration and job execution. It is a structure proposal, not a choice of scientific method.

## Current files

```text
VisHD/
├── .gitignore                     # Lee's uv
├── .python-version                # Lee's uv
├── Package outline.md            # added on this branch
├── README.md                     # Lee's uv, unchanged
├── pyproject.toml                # Lee's uv, updated here
├── uv.lock                       # Lee's uv, updated here
├── config/
│   └── vishd.example.toml        # added here
├── src/VisHD/
│   ├── __init__.py               # Lee's uv, updated here
│   ├── blade.py                  # Lee's uv
│   ├── SEAC.py                   # Lee's uv
│   ├── Visualize.py              # Lee's uv, updated here
│   ├── Permutation_Test.py       # Lee's uv
│   ├── Permutation_Runner.py     # Lee's uv
│   ├── PermuteAllSamples.py      # Lee's uv
│   ├── permutation_pipeline.py   # Lee's uv, updated here
│   ├── evalPermutations.py       # Lee's uv
│   ├── make_genecat_table.py     # Lee's uv
│   ├── config.py                 # added here
│   └── execution/                # added here
│       ├── __init__.py           # added here
│       ├── local_backend.py      # added here
│       └── slurm.py              # added here
└── tests/
    ├── test_config.py            # added here
    └── test_execution.py         # added here
```

## Current state

- Lee's pipeline now takes paths and compute settings from a TOML file or Python object. It runs locally by default; SLURM requires the optional `submitit` extra.
- The older experimental modules still use the `local.py` convention. Mac's scripts have not been added to this repository.
- The GPU extra is declared but not implemented yet.
- Tests cover configuration, local job execution, and basic imports. They do not run a full sample or check biological results.

## Proposed next modules

Lee's pipeline uses BLADE positions and a KD-tree. Mac's separate `cpu-visium.py` and `visium.py` use Space Ranger counts and dense-grid distance transforms. Proposed `readers/` would convert platform data into one input format, proposed `methods/` would run either analysis on that input and return one result format. `results.py` could save the results. These modules do not exist yet.

Do these module boundaries fit how everyone expects to contribute? Should the first shared input representation be a table of `(sample, x, y, gene, count)`, or is another shape better for both Visium HD and CosMx?

Also, what name do we want for the package? Since we're serving both Visium HD and CosMx, we should consider something different from VisHD.
