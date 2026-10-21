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

## Using from R

With [reticulate](https://rstudio.github.io/reticulate/), point R at the
project's virtual environment before anything else loads Python:

```r
library(reticulate)
use_virtualenv("~/VisHD/.venv", required = TRUE)

SEAC     <- import("VisHD.SEAC")
pipeline <- import("VisHD.permutation_pipeline")

pipeline$run_pipeline(samples = c("H1-QYNNNRX_A1", "H1-QYNNNRX_D1"),
                      slurm_account = "...")
```

reticulate converts pandas DataFrames but not polars ones; call
`df$to_pandas()` on a polars DataFrame to bring it into R. For pipeline
outputs, it is simpler to read the parquet files directly with the `arrow`
package:

```r
thresholds <- arrow::read_parquet("<projectroot>/collated/H1-QYNNNRX_A1.parquet")
```

## Development

```bash
uv sync                                          # create .venv with VisHD installed in editable mode
uv run python -m VisHD.permutation_pipeline      # run a module as a script
```
