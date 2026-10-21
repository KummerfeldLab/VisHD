# VIsHD

## Installation

Install directly from GitHub; dependencies are installed along with it.

```bash
uv add git+https://github.com/KummerfeldLab/VisHD     # in your own uv project
pip install git+https://github.com/KummerfeldLab/VisHD
```

Then use it like any library:

```python
from VisHD import SEAC
from VisHD.permutation_pipeline import run_pipeline
run_pipeline(samples=[...], slurm_account="...")
```

VisHD expects a machine-specific `VisHD/local.py` (not tracked in git) that
defines `bladed_d` and `spaceranger_outputs_dir` as `pathlib.Path`s.

## Development

```bash
uv sync                                          # create .venv with VisHD installed in editable mode
uv run python -m VisHD.permutation_pipeline      # run a module as a script
```
