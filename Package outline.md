# VisHD package structure (draft)

This branch starts from Lee's `uv` branch and adds a first layer for configuration and job execution. It is a structure proposal for the team, not a choice of scientific method.

The main README is unchanged from `uv`; it still describes the earlier `local.py` and SLURM workflow. This outline describes the proposed direction for review.

```text
VisHD/
├── pyproject.toml, uv.lock       # Lee's packaging foundation
├── config/vishd.example.toml    # example user settings
├── src/VisHD/
│   ├── config.py                # paths and execution settings
│   ├── execution/               # local or SLURM job launcher
│   ├── permutation_pipeline.py  # Lee's current pipeline, adapted to use these settings
│   ├── readers/                 # proposed: Visium HD and CosMx input adapters
│   ├── methods/                 # proposed: Lee and Mac implementations behind a common interface
│   └── results.py               # proposed: shared output writing
└── tests/                        # small package tests
```

The `readers/`, `methods/`, and `results.py` entries are proposed locations; they have not been implemented. Existing modules stay in place during the transition.

## Why these boundaries?

- **User settings live outside the installed package.** A TOML file or Python object supplies paths and compute settings. This avoids editing a machine-specific `local.py` inside package code.
- **Job execution is separate from the analysis method.** Local execution works without a cluster; SLURM is an optional extra. CPU/GPU selection is a separate concern from local/SLURM job submission.
- **Readers and methods need separate interfaces.** Lee's pipeline uses BLADE position tables and KD-tree distances. Mac's separately supplied `cpu-visium.py` and `visium.py` use Space Ranger count matrices, dense grids, and CPU/GPU distance transforms. A shared input representation and output contract would let us compare them without assuming either is the final method. Mac's scripts have not been copied into this branch.
- **GPU integration is pending.** The optional GPU dependency is declared, but the public pipeline reports that GPU execution is not connected yet.

The current tests cover configuration, local job execution, and basic package import. They do not validate biological results or run a full sample.

**Feedback requested:** Do these module boundaries fit how Lee, Mac, and Yin expect to contribute? Should the first shared input representation be a table of `(sample, x, y, gene, count)`, or is another shape better for both Visium HD and CosMx?
