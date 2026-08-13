import getpass
import submitit
import datetime
from pathlib import Path
from VisHD import local, PermuteAllSamples, SEAC
import scipy
import gc
import polars as pl
import numpy as np
from zoneinfo import ZoneInfo


samples = [
    "H1-QYNNNRX_A1",
    "H1-QYNNNRX_D1",
    "H1-TMNBZ4N_A1",
    "H1-TMNBZ4N_D1",
    "H1-TQXC9QB_A1",
    "H1-TQXC9QB_D1",
    "H1-HDMGRMD_A1",
]

GENE_DICT = SEAC.SEAC
PERMS_PER_RUN = 100
POSITIONS_DIR = local.bladed_d


def permute_sample(samplename, runnum, topdir, count_table):
    """Provides MaxMin evaluations of permutations for one sample of the genes of all samples"""

    start = runnum * PERMS_PER_RUN
    end = start + PERMS_PER_RUN

    # Set/Create output directory
    outdir = topdir / sample
    outdir.mkdir(exist_ok=True)

    # Read samples longform barcodes
    lf_path = topdir / f"longforms/{samplename}.parquet"
    longform = pl.read_parquet(lf_path)

    # Explode all the gene categories
    null_perm = []
    for cat, n in zip(count_table["category"], count_table["Total"]):
        null_perm += [cat] * n

    # Determine the indices of the permutation belonging to all samples
    indices = {}
    cumsum = 0
    for s in count_table.columns[1:-1]:
        begin = cumsum
        cumsum = begin + count_table[s].sum()
        indices[s] = (begin, cumsum)

    # Determine the indices of the permutaiton belonging to this sample
    start_idx, end_idx = indices[samplename]

    # read the positions data
    positions = pl.read_parquet(POSITIONS_DIR / f"{samplename}_bladed.parquet")
    results = positions.select("barcode")

    # Permute the base list of gene categories and save as a column
    for i in range(start, end):
        print(i, flush=True)
        permutation = np.random.default_rng(seed=i).permuted(null_perm)[
            start_idx:end_idx
        ]
        results = results.with_columns(
            pl.lit(
                PermuteAllSamples.evaluate_permutation(longform, permutation, positions)
            ).alias(str(i))
        )

    results.write_parquet(outdir / f"{samplename}_{start}_{end}.parquet")


def make_gene_count_table(samples, longformdir, gene_dict=GENE_DICT):
    """Reads each sample in `samples`, and returns a sample counting the number
    of genes in the values of each key of gene_dict as well as the number of non-category genes.
    Also outputs the longform barcodes to the longform directory """
    table = None
    for sample in samples:
        sampledir = local.spaceranger_outputs_dir / sample
        adata = SEAC.read_h5(sampledir / "outs/binned_outputs/square_002um")
        pos = pl.read_parquet(local.bladed_d / f"{sampledir.name}_bladed.parquet")

        # Create a long form table
        row, col, count = scipy.sparse.find(adata.X)
        df = pl.DataFrame(
            {
                "barcode": adata.obs_names[row],
                "gene": adata.var_names[col],
                "nreads": count,
            }
        )
        barcodes = pos.filter(pl.col("blade").eq(0))["barcode"].to_list()
        df = df.filter(pl.col("barcode").is_in(barcodes))
        long_form = df.select(pl.exclude("nreads").repeat_by("nreads").explode())

        # Save long form barcodes
        long_form.select('barcode').write_parquet(longformdir / f"{sample}.parquet")

        # Build a dynamic when expression
        when_expr = pl.coalesce(
            pl.when(pl.col("gene").is_in(v)).then(pl.lit(k))
            for k, v in gene_dict.items()
        )

        category_counts = long_form.select(
            when_expr.fill_null(pl.lit("Other")).alias("category")
        )["category"].value_counts()
        category_counts = category_counts.select(
            "category", pl.col("count").alias(sample.name)
        )

        if table is None:
            table = category_counts
        else:
            table = table.join(category_counts, on="category")

    return table


def collate(samplename, permutations, outdir):
    df = pl.DataFrame()
    for nfiles, f in enumerate(permutations.iterdir(), start=1):
        gc.collect()
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

    df.write_parquet(outdir / f"{samplename}.parquet")


if __name__ == "__main__":
    username = getpass.getuser()
    chitime = datetime.datetime.now(tz=ZoneInfo("America/Chicago"))
    datestamp = chitime.strftime("%Y%m%d")
    print(username, datestamp)

    # Create Permuation Output Directories
    userroot = Path("/sbatch.global") / username
    userroot.mkdir(exist_ok=True)

    projectroot = Path("/sbatch.global") / username
    projectroot.mkdir(exist_ok=True)

    collated_dir = projectroot / "collated"
    collated_dir.mkdir(exist_ok=True)

    image_dir = projectroot / "images"
    image_dir.mkdir(exist_ok=True)

    longformdir = projectroot / "longform"
    longformdir.mkdir(exist_ok=True)

    log_dir = projectroot / "logs"
    log_dir.mkdir(exist_ok=True)

    for sample in samples:
        (projectroot / sample).mkdir(exist_ok=True)

    # Make Gene Category table
    count_table = make_gene_count_table(samples, longformdir, GENE_DICT)

    # Do the Permutations
    executor = submitit.AutoExecutor(folder=log_dir)
    executor.update_parameters(
        timeout_min=360,
        mem_gb=64,
        slurm_account="lniedern",
        cpus_per_task=2,
        name="VisHDPermute",
    )

    # Collate

    # Generate Image Files
