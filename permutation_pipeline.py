import getpass
import submitit
import datetime
from zoneinfo import ZoneInfo
import time
from pathlib import Path
from VisHD import local, SEAC
import VisHD.Visualize as viz
import scipy
import gc
import polars as pl
import numpy as np
from scipy.spatial import KDTree
from enum import Enum
import pyarrow
import pyarrow.parquet as pq


samples = [
    "H1-QYNNNRX_A1",
    "H1-QYNNNRX_D1",
    "H1-TMNBZ4N_A1",
    "H1-TMNBZ4N_D1",
    "H1-TQXC9QB_A1",
    "H1-TQXC9QB_D1",
    "H1-HDMGRMD_A1",
    "H1-HDMGRMD_D1",
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
    perms_per_run = runinfo.get("perms_per_run", PERMS_PER_RUN)
    positions_dir = runinfo.get("positions_dir", POSITIONS_DIR)
    start = runnum * perms_per_run
    end = start + perms_per_run
    threshold = 100

    # Set/Create output directory
    outdir = topdir / samplename
    outdir.mkdir(exist_ok=True)

    # Read samples longform barcodes, this tells us where the reads were found, and how many
    lf_path = topdir / f"longform/{samplename}.parquet"
    longform = pl.read_parquet(lf_path)

    # read the positions data to get barcodes for the results file
    positions = pl.read_parquet(positions_dir / f"{samplename}_bladed.parquet")
    # restrict positions to spots not set for blade and are in tissue
    positions = positions.filter(pl.col('blade').eq(0),
                                 pl.col('in_tissue').eq(1))

    # make a reduced copy of positions to contain only candidate hotspot rows
    # hotspots must have a MaxMin distance less than threshold
    hotspot_pos = positions.filter(pl.col("MaxMin").le(threshold))
    # make a numpy array of the hotspot coordinates
    hotspot_coords = hotspot_pos.select("array_row", "array_col").to_numpy()
    # create a new table for storing distance results
    results = hotspot_pos.select("barcode")
    # barcode -> coordinates lookup for the permuted read locations
    coords_by_barcode = positions.select("barcode", "array_row", "array_col")

    # lookup how many SEAC gene hits are in the sample for each category
    nhits = {
        category: count_table.filter(pl.col("category").eq(category))[samplename].item()
        for category in SEAC.SEAC.keys()
    }

    # From start to end, make a table of permuted MaxMin values
    for i in range(start, end):
        print(i, flush=True)

        # shuffle all read locations once using i as the seed; consecutive slices
        # of the shuffle are assigned to each category (sampling without replacement)
        order = np.random.default_rng(seed=i).permutation(longform.height)
        offset = 0
        distances = []
        for category, n in nhits.items():
            # select barcodes of the read locations assigned to this category
            seacgene_idxs = order[offset:offset + n]
            offset += n
            seacgene_barcodes = longform[seacgene_idxs].select("barcode").unique()
            # select coordinates of those barcodes
            seacgene_coords = (
                coords_by_barcode.join(seacgene_barcodes, on="barcode")
                .select("array_row", "array_col")
                .to_numpy()
            )

            # build tree of seacgene_coords and get distances to them from hotspots
            tree = KDTree(seacgene_coords)
            d, _ = tree.query(hotspot_coords, workers=-1)
            distances.append(d)
        # MaxMin: the max over categories of the min distance to each category
        results = results.with_columns(pl.Series(str(i), np.max(distances, axis=0)))

    results.write_parquet(outdir / f"{samplename}_{start}_{end}.parquet")


def make_gene_count_table(
    samples,
    longformdir,
    gene_dict=GENE_DICT,
    positions_dir=POSITIONS_DIR,
    spaceranger_dir=local.spaceranger_outputs_dir,
):
    """Reads each sample in `samples`, and returns a sample counting the number
    of genes in the values of each key of gene_dict as well as the number of non-category genes.
    Also outputs the longform barcodes to the longform directory"""
    table = pl.DataFrame()
    for samplename in samples:
        print(f"\t{timestamp()} {samplename}: Reading Data")
        sampledir = spaceranger_dir / samplename
        adata = SEAC.read_h5(sampledir / "outs/binned_outputs/square_002um")
        pos = pl.read_parquet(positions_dir / f"{sampledir.name}_bladed.parquet")

        # Create a long form table
        print(f"\t{timestamp()} {samplename}: Reformatting")
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
        long_form.select("barcode").write_parquet(longformdir / f"{samplename}.parquet")

        # Build a dynamic when expression
        when_expr = pl.coalesce(
            pl.when(pl.col("gene").is_in(v)).then(pl.lit(k))
            for k, v in gene_dict.items()
        )

        category_counts = long_form.select(
            when_expr.fill_null(pl.lit("Other")).alias("category")
        )["category"].value_counts()
        category_counts = category_counts.select(
            "category", pl.col("count").alias(samplename)
        )

        if len(table) == 0:
            table = category_counts
        else:
            table = table.join(category_counts, on="category")

    table = table.with_columns(table[:, 1:].sum_horizontal().alias("Total"))
    table.write_parquet(longformdir.parent / "count_table.parquet")

    return table


def newcollate(jobinfo):
    samplename = jobinfo["samplename"]
    topdir = jobinfo["projectroot"]
    outdir = jobinfo["outdir"]
    n_lowest = jobinfo.get("n_lowest", 50)
    permutation_dir = topdir / samplename

    barcodes = []
    chunk = []
    for nfiles, f in enumerate(permutation_dir.iterdir(), start=1):
        print(f"{timestamp()} reading {f}", flush=True)
        table = pq.read_table(f)
        if nfiles == 1:
            print(f"{timestamp()} extracting barcodes", flush=True)
            barcodes = table.column("barcode")
            barcodes.to_numpy()
        print(f"{timestamp()} extracting columns", flush=True)
        array = np.column_stack(
            [table.column(name).to_numpy() for name in table.column_names[1:]]
        )
        del table
        gc.collect()
        if nfiles == 1:
            chunk = array
        else:
            chunk = np.hstack([chunk, array])
        print(f"{timestamp()} selecting lowest {n_lowest}", flush=True)
        chunk = np.partition(chunk, n_lowest, axis=1)[:, :n_lowest]
    print(f"{timestamp()} selecting row-wise max", flush=True)
    thresholds = chunk.max(axis=1)
    result = pyarrow.Table.from_arrays(
        [barcodes, thresholds], names=["barcode", "threshold"]
    )
    print(f"{timestamp()} writing results", flush=True)
    pyarrow.parquet.write_table(result, outdir / f"{samplename}.parquet")


def make_annotated_image(
    samplename, gene_list, collated_dir, img_outdir, positions_dir=POSITIONS_DIR
):
    positions_filename = positions_dir / f"{samplename}_bladed.parquet"
    print(f"{timestamp()} Reading positions data from {positions_filename}", flush=True)
    pos = pl.read_parquet(positions_filename)
    pos = pos.filter(pl.col("MaxMin").le(100), pl.col("blade").eq(0))
    print(f"{timestamp()} Positions file has columns {pos.columns}", flush=True)

    collated_f = collated_dir / f"{samplename}.parquet"
    print(f"{timestamp()} Reading collated permutations from {collated_f}", flush=True)
    collated_permutations = pl.read_parquet(collated_f)
    print(
        f"{timestamp()} Collated file has columns {collated_permutations.columns}",
        flush=True,
    )

    print("Calculating ROIs:", flush=True)
    # find spots whose MaxMin is below the max permutation
    # and extract their barcodes and coordinates
    bcs = {}
    #   collated_permutations = collated_permutations.with_columns(
    #       pl.col("low").list.max().alias("threshold")
    #   )
    print(f"{timestamp()} Joinings positions file with thresholds", flush=True)
    joined = collated_permutations.join(pos, on="barcode")
    print(f"{timestamp()} joined dataframe has columns {joined.columns}", flush=True)
    hotspots_pos = joined.filter(pl.col("MaxMin").le(pl.col("threshold")))
    bcs["hotspots"] = hotspots_pos["barcode"].to_list()
    hotspot_coords = hotspots_pos.select("array_row", "array_col").to_numpy()

    # Find closest category read for each hotspot
    for category in gene_list:
        print(category, end=" ")
        category_pos = joined.filter(pl.col(category).eq(0))
        category_coords = category_pos.select("array_row", "array_col").to_numpy()
        tree = KDTree(category_coords)
        _, idxs = tree.query(hotspot_coords, workers=-1)
        bcs[category] = category_pos[np.unique(idxs)]["barcode"].to_list()

    print("\nImaging", end=" ")
    viz.make_multi_annotated_image(
        samplename, joined, bcs, colors=viz.COLORS, outdir=img_outdir
    )
    print("Done", end="\n\n")


def run_pipeline(
    samples=samples,
    gene_dict=GENE_DICT,
    projectroot=None,
    scratch_root=Path("/scratch.global"),
    positions_dir=POSITIONS_DIR,
    spaceranger_dir=local.spaceranger_outputs_dir,
    n_runs=50,
    perms_per_run=PERMS_PER_RUN,
    n_lowest=50,
    slurm_account="lniedern",
    mem_gb=64,
    count_table_time="30:00",
    permute_time="06:00:00",
    permute_cpus=4,
    collate_time="60:00",
    collate_cpus=2,
    image_time="20:00",
    permute_all_samples=False,
):
    """Run the full permutation pipeline on slurm.

    Phases: build the gene category count table, permute each sample,
    collate the permutations into per-barcode thresholds, and image hotspots.
    Each sample gets n_runs * perms_per_run permutations; the per-barcode
    threshold is the max of the n_lowest permuted MaxMin values.

    Parameters
    ----------
    samples : list of str
        Sample names to process.
    gene_dict : dict
        Maps category name -> list of genes. Genes in no category are
        counted as "Other".
    projectroot : Path or None
        Output directory for all pipeline products. If None, defaults to
        <scratch_root>/<user>/permutation_<YYYYMMDD> (Chicago time).
    scratch_root : Path
        Parent of the per-user directory; only used when projectroot is None.
    positions_dir : Path
        Directory of <sample>_bladed.parquet position files.
    spaceranger_dir : Path
        Directory of spaceranger outputs, one subdirectory per sample.
    n_runs : int
        Number of permutation jobs per sample.
    perms_per_run : int
        Number of permutations computed by each job.
    n_lowest : int
        Number of lowest permuted MaxMin values per barcode used to set the
        threshold. Keep at ~1% of n_runs * perms_per_run.
    slurm_account : str
        Slurm account to charge jobs to.
    mem_gb : int
        Memory (GB) requested for each job.
    count_table_time : str
        Slurm time limit for the count table job.
    permute_time : str
        Slurm time limit for each permutation job.
    permute_cpus : int
        CPUs per permutation job.
    collate_time : str
        Slurm time limit for each collate job.
    collate_cpus : int
        CPUs per collate job.
    image_time : str
        Slurm time limit for each image job.
    permute_all_samples : bool
        If True, permute gene categories pooled across all samples; if False,
        permute each sample's own gene categories only.
    """
    # Create Permuation Output Directories
    print(timestamp(), "creating directories")
    if projectroot is None:
        username = getpass.getuser()
        chitime = datetime.datetime.now(tz=ZoneInfo("America/Chicago"))
        datestamp = chitime.strftime("%Y%m%d")
        userroot = Path(scratch_root) / username
        userroot.mkdir(exist_ok=True)
        projectroot = userroot / f"permutation_{datestamp}"
    projectroot = Path(projectroot)
    projectroot.mkdir(parents=True, exist_ok=True)

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
        slurm_account=slurm_account, mem_gb=mem_gb, slurm_time=count_table_time
    )

    # Make Gene Category table
    if not (projectroot / "count_table.parquet").exists():
        print(timestamp(), "building count table")
        executor.update_parameters(
            name="Make_Count_Table",
            #           slurm_output="Make_Count_Table.log",
        )

        job = executor.submit(
            make_gene_count_table,
            samples,
            longformdir,
            gene_dict,
            positions_dir,
            spaceranger_dir,
        )
        job.wait()
        count_table = pl.read_parquet(projectroot / "count_table.parquet")
        print(timestamp(), "count table finished")
    else:
        count_table = pl.read_parquet(projectroot / "count_table.parquet")
    print(count_table)

    # Create  Dictionary of Run Data
    runs = []
    for samplename in samples:
        for i in range(n_runs):
            perm_d = {
                "samplename": samplename,
                "i": i,
                "projectroot": projectroot,
                "count_table": count_table,
                "perms_per_run": perms_per_run,
                "positions_dir": positions_dir,
                "permute_all_samples": permute_all_samples,
            }
            runs.append(perm_d)

    # Make Permutations for Each Sample
    print(timestamp(), "Permuting Samples")
    executor.update_parameters(
        name="VisHDPermute",
        slurm_job_name="VisHDPermute",
        slurm_time=permute_time,
        cpus_per_task=permute_cpus,
    )
    jobs = executor.map_array(permute_sample, runs)
    _ = [job.wait() for job in jobs]
    print(timestamp(), "Permuting Finished")
    print("Sanity Check:")
    for samplename in samples:
        files = [
            f for f in (projectroot / samplename).iterdir() if ".parquet" in f.name
        ]
        print(samplename, len(files))

    # Collate Permutations of Each Sample and Determine 1% Threshold
    runs = []
    for sample in samples:
        collate_d = {
            "samplename": sample,
            "projectroot": projectroot,
            "outdir": collated_dir,
            "n_lowest": n_lowest,
        }
        runs.append(collate_d)

    print(timestamp(), "Collating Permutations")
    executor.update_parameters(
        name="CollatePermutations",
        mem_gb=mem_gb,
        cpus_per_task=collate_cpus,
        slurm_time=collate_time,
        #       slurm_output=f"{collated_dir}/collate_%a.out",
    )
    jobs = executor.map_array(newcollate, runs)
    _ = [job.wait() for job in jobs]
    print(timestamp(), "Collating Finished")

    # Generate Image Files
    print(timestamp(), "Generating Images")
    executor.update_parameters(
        name="GenerateImages",
        slurm_time=image_time,
    )
    jobs = []
    with executor.batch():
        for samplename in samples:
            job = executor.submit(
                make_annotated_image,
                samplename,
                gene_dict,
                collated_dir,
                image_dir,
                positions_dir,
            )
            jobs.append(job)
    _ = [job.wait() for job in jobs]
    print(timestamp(), "Finished")


if __name__ == "__main__":
    run_pipeline()
