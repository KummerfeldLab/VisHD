import sys
from pathlib import Path

from VisHD import Permutation_Test, local

snames = [x.name.split("_bladed")[0] for x in local.bladed_d.iterdir()]
results_dir = Path.home() / "All Projects/Adipose/MaxMinSEAC/PermutationOutput"

if __name__ == "__main__":
    start = 800
    end = 900
    sample_name = snames[int(sys.argv[-1])]
    print(snames[int(sys.argv[-1])], start, end)
    sdir = local.spaceranger_outputs_dir / sample_name
    result = Permutation_Test.permute_SEAC_sample(sdir, start, end)
    result.write_parquet(results_dir / f"{sample_name}_{start}_{end}.parquet")
    print("finished")
