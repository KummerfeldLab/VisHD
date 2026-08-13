import getpass
import submitit
from datetime import date
from pathlib import Path
from VisHD import local, PermuteAllSamples
from VisHD.SEAC import SEAC
import gc
import polars as pl


samples = [
    "H1-QYNNNRX_A1",
    "H1-QYNNNRX_D1",
    "H1-TMNBZ4N_A1",
    "H1-TMNBZ4N_D1",
    "H1-TQXC9QB_A1",
    "H1-TQXC9QB_D1",
    "H1-HDMGRMD_A1",
]

gene_dict = SEAC


def collate(sample, permutations):
    for nfiles, f in enumerate(permutations.iterdir(), start=1):
        gc.collect()
        df = pl.DataFrame()
        if nfiles == 1:
            df = pl.read_parquet(f)
        else:
            df = df.join(pl.read_parquet(f), on="barcode")

    df = df.select(
        pl.col("barcode"),
        pl.concat_list(pl.exclude("barcode"))
        .list.sort()
        .list.head(50)
        .list.max()
        .alias("threshold"),
    )
    gc.collect()

    return df


if __name__ == "__main__":
    username = getpass.getuser()
    today = date.today().strftime("%Y%m%d")
    print(username, today)

    # Create Permuation Output Directories
    projectroot = Path("/sbatch.global") / username
    projectroot.mkdir(exist_ok=True)

    collated_dir = projectroot / "collated"
    collated_dir.mkdir(exist_ok=True)

    image_dir = projectroot / "images"
    image_dir.mkdir(exist_ok=True)

    log_dir = projectroot / "logs"
    log_dir.mkdir(exist_ok=True)

    for sample in samples:
        (projectroot / sample).mkdir(exist_ok=True)

    # Make Gene Category table
    count_table = PermuteAllSamples.make_gene_count_table(samples, gene_dict)

    # Do the Permutations
    executor = submitit.AutoExecutor(folder=log_dir)
    executor.update_parameters(
        timeout_min=360,
        mem_gb=64,
        slurm_account="lniedern",
    )

    # Collate

    # Generate Image Files
