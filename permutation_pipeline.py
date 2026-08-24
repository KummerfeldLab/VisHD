import getpass
import submitit
import datetime
from zoneinfo import ZoneInfo
import time
from pathlib import Path
from VisHD import local, PermuteAllSamples, SEAC
import VisHD.Visualize as viz
import scipy
import gc
import polars as pl
import numpy as np
from scipy.spatial import KDTree
from enum import Enum


samples = [
    "H1-QYNNNRX_A1",
    "H1-QYNNNRX_D1",
    "H1-TMNBZ4N_A1",
    "H1-TMNBZ4N_D1",
    "H1-TQXC9QB_A1",
    "H1-TQXC9QB_D1",
    "H1-HDMGRMD_A1",
]

class PHASE(Enum):
    ALL = 0
    SETUP = 1 
    MAKE_COUNT_TABLE = 2
    PERMUTE = 3
    COLLATE = 4
    VISUALIZE = 5

GENE_DICT = SEAC.SEAC
PERMS_PER_RUN = 100
POSITIONS_DIR = local.bladed_d



def timestamp():
    return time.strftime("(%H:%M:%S)")


def permute_sample(runinfo):
    """Provides MaxMin evaluations of multiple permutations for one sample of the genes of all samples"""

    samplename = runinfo["samplename"]
    runnum = runinfo["i"]
    topdir = runinfo["projectroot"]
    count_table = runinfo["count_table"]
    start = runnum * PERMS_PER_RUN
    end = start + PERMS_PER_RUN

    # Set/Create output directory
    outdir = topdir / samplename
    outdir.mkdir(exist_ok=True)

    # Read samples longform barcodes
    lf_path = topdir / f"longform/{samplename}.parquet"
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
    Also outputs the longform barcodes to the longform directory"""
    table = pl.DataFrame()
    for sample in samples:
        print("\t{timestamp()} {sample}: Reading Data")
        sampledir = local.spaceranger_outputs_dir / sample
        adata = SEAC.read_h5(sampledir / "outs/binned_outputs/square_002um")
        pos = pl.read_parquet(local.bladed_d / f"{sampledir.name}_bladed.parquet")

        # Create a long form table
        print("\t{timestamp()} {sample}: Reformatting")
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
        long_form.select("barcode").write_parquet(longformdir / f"{sample}.parquet")

        # Build a dynamic when expression
        when_expr = pl.coalesce(
            pl.when(pl.col("gene").is_in(v)).then(pl.lit(k))
            for k, v in gene_dict.items()
        )

        category_counts = long_form.select(
            when_expr.fill_null(pl.lit("Other")).alias("category")
        )["category"].value_counts()
        category_counts = category_counts.select(
            "category", pl.col("count").alias(sample)
        )

        if len(table) == 0:
            table = category_counts
        else:
            table = table.join(category_counts, on="category")

    table = table.with_columns(table[:, 1:].sum_horizontal().alias("Total"))
    table.write_parquet(longformdir.parent / "count_table.parquet")

    return table


def collate(samplename, topdir, outdir):
    df = pl.DataFrame()
    permutation_dir = topdir / samplename
    for nfiles, f in enumerate(permutation_dir.iterdir(), start=1):
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


def lazy_collate(jobinfo):
    samplename = jobinfo["samplename"]
    topdir = jobinfo["projectroot"]
    outdir = jobinfo["outdir"]
    permutation_dir = topdir / samplename
    print(f"{timestamp()} Collating {samplename}")
    print(f"Source dirctory: {permutation_dir}")

#   lazy_frames = [pl.scan_parquet(f) for f in permutation_dir.iterdir()]
    files = list(permutation_dir.iterdir())
    lazy_frames = [pl.scan_parquet(files[0])] + [
        pl.scan_parquet(f).drop("barcode") for f in files[1:]
    ]
    print(f"Found {len(lazy_frames)} files")

#   merged = lazy_frames[0]
#   for lf in lazy_frames[1:]:
#       merged = merged.join(lf, on="barcode", how="inner")

    merged = pl.concat(lazy_frames, how='horizontal')


    result = merged.select(
        pl.col("barcode"),
        pl.concat_list(pl.exclude("barcode"))
        .list.sort()
        .list.head(50)
        .list.max()
        .alias("threshold"),
    )

    print(f"writing {outdir}/{samplename}.parquet")
    result.sink_parquet(outdir / f"{samplename}.parquet")
    print(f"{timestamp()} Finished Collating {samplename}")


def make_annotated_image(samplename, gene_list, collated_dir, img_outdir):
    gc.collect()
    positions_filename = local.bladed_d / f"{samplename}_bladed.parquet"
    print(positions_filename, "Reading", end=" ")
    pos = pl.read_parquet(positions_filename)
    tpos = pos.filter(pl.col("MaxMin").le(100), pl.col("blade").eq(0))
    collated_permutations = pl.read_parquet(collated_dir / f"{samplename}.parquet")

    print("Calculating ROIs:", end=" ")
    # find spots whose MaxMin is below the max permutation
    # and extract their barcodes and coordinates
    bcs = {}
    collated_permutations = collated_permutations.with_columns(
        pl.col("low").list.max().alias("threshold")
    )
    joined = collated_permutations.join(tpos, on="barcode")
    hotspots_pos = joined.filter(pl.col("MaxMin").le(pl.col("threshold")))
    bcs["hotspots"] = hotspots_pos["barcode"].to_list()
    hotspot_coords = hotspots_pos.select("array_row", "array_col").to_numpy()

    # Find closest category read for each hotspot
    for category in gene_list:
        print(category, end=" ")
        category_pos = tpos.filter(pl.col(category).eq(0))
        category_coords = category_pos.select("array_row", "array_col").to_numpy()
        tree = KDTree(category_coords)
        _, idxs = tree.query(hotspot_coords, workers=-1)
        bcs[category] = category_pos[np.unique(idxs)]["barcode"].to_list()

    print("\nImaging", end=" ")
    viz.make_multi_annotated_image(
        samplename, tpos, bcs, colors=viz.COLORS, outdir=img_outdir
    )
    print("Done", end="\n\n")


if __name__ == "__main__":
    username = getpass.getuser()
    chitime = datetime.datetime.now(tz=ZoneInfo("America/Chicago"))
    datestamp = chitime.strftime("%Y%m%d")

    # Create Permuation Output Directories
    print(timestamp(), "creating directories")
    userroot = Path("/scratch.global") / username
    userroot.mkdir(exist_ok=True)

#   projectroot = userroot / f"permutation_{datestamp}"
    projectroot = userroot / "permutation_20260820"
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

    # Instantiate Job executor
    executor = submitit.AutoExecutor(folder=log_dir)
    executor.update_parameters(
        slurm_account="lniedern",
        mem_gb=64,
        slurm_time="30:00"
        )

    # Make Gene Category table
    if not (projectroot / "count_table.parquet").exists():
        print(timestamp(), "building count table")
        executor.update_parameters(
            name="Make_Count_Table",
#           slurm_output="Make_Count_Table.log",
        )

        job = executor.submit(make_gene_count_table, samples, longformdir, GENE_DICT)
        count_table = job.result()
        print(timestamp(), "count table finished")
    else:
        count_table = pl.read_parquet(projectroot / "count_table.parquet")
    print(count_table)

    # Create  Dictionary of Run Data
    runs = []
    for samplename in samples:
        for i in range(50):
            perm_d = {
                "samplename": samplename,
                "i": i,
                "projectroot": projectroot,
                "count_table": count_table,
            }
            runs.append(perm_d)

    # Make Permutations for Each Sample
    print(timestamp(), "Permuting Samples")
    executor.update_parameters(
        name="VisHDPermute",
        slurm_job_name="VisHDPermute",
#       slurm_output="permute_%a.out",
        slurm_time="06:00:00",
        cpus_per_task=4,
    )
#   jobs = executor.map_array(permute_sample, runs)
#   _ = [job.wait() for job in jobs]
    print(timestamp(), "Permuting Finished")
    print("Sanity Check:")
    for samplename in samples:
        files = [f for f in (projectroot / samplename).iterdir() if '.parquet' in f.name]
        print(samplename, len(files))

    # Collate Permutations of Each Sample
    runs = []
    for sample in samples:
        collate_d = {
            "samplename": sample,
            "projectroot": projectroot,
            "outdir": collated_dir,
        }
        runs.append(collate_d)

    print(timestamp(), "Collating Permutations")
    executor.update_parameters(
        name="CollatePermutations",
        mem_gb=64,
        cpus_per_task=2,
#       slurm_output=f"{collated_dir}/collate_%a.out",
    )
    jobs = executor.map_array(lazy_collate, runs)
    _ = [job.result() for job in jobs]
    print(timestamp(), "Collating Finished")

    # Generate Image Files
    print(timestamp(), "Generating Images")
    executor.update_parameters(name="GenerateImages")
    jobs = []
    with executor.batch():
        for samplename in samples:
            job = executor.submit(
                make_annotated_image, samplename, GENE_DICT, collated_dir, image_dir
            )
            jobs.append(job)
    _ = [job.result() for job in jobs]
    print(timestamp(), "Finished")
