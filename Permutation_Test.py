import multiprocessing as mp

import numpy as np
import polars as pl
import scipy
from VisHD import SEAC, blade, local
from scipy.spatial import cKDTree


def read_basedata(sampledir):
    sample = SEAC.read_h5(sampledir / "outs/binned_outputs/square_002um")
    positions = pl.read_parquet(local.bladed_d / f"{sampledir.name}_bladed.parquet")
    return sample, positions


def permute_SEAC_sample(sampledir, start=0, end=50):
    sample, positions = read_basedata(sampledir)
    long_form = adata2longform(sample, positions)
    result = None
    for i in range(start, end):
        if result is None:
            result = evaluate_permutation(long_form, positions, str(i), seed=i)
        else:
            result = result.join(
                evaluate_permutation(long_form, positions, str(i), seed=i), on="barcode"
            )

    return result


def permute_all_samples(sampledir, seed=None):
    long_form = None
    all_pos = None

    # Build up data structure
    for sample in sampledir.iterdir():
        if not sample.is_dir():
            continue
        print("Reading Sample {sample.name}")
        adata, pos = read_basedata(sample)
        lf = adata2longform(adata, pos)
        lf = lf.with_columns(pl.lit(sample.name).alias('sample'))
        pos = pos.with_columns(pl.lit(sample.name).alias('sample'))
        print("Concatenating longform, positions")
        if long_form is None:
            long_form = lf
            all_pos = pos
        else:
            long_form = pl.concat([long_form, lf])
            all_pos = pl.concat([all_pos, pos])
        del adata, pos
    return

    # Process permutation
    print("Permuting Long Form")
    long_form = long_form.with_columns(pl.col("gene").shuffle(seed=seed).alias("gene"))
#   allperms = None
    for sample in sampledir.iterdir():
        if not sample.is_dir():
            continue
        print("Evaluating Permutation of {sample}")
        perm = evaluate_permutation(long_form.filter(pl.col("sample").eq(sample.name)),
                                    all_pos.filter(pl.col("sample").eq(sample.name)),
                                    seed=seed,
                                    shuffle=False)
#       perm = perm.with_columns(pl.lit(sample.name).alias('sample'))
#       if allperms is None:
#           allperms = perm
#       else:
#           allperms = pl.concat([allperms, perm])

        # Save to scratch
        print("Writing File")
        perm.write_parquet(f'/scratch.global/will4379/{sample}/perms_{seed}.parquet')

 

def evaluate_permutation(long_form,
                         positions,
                         newcolname="MaxMin",
                         seed=None,
                         shuffle=True):
    """
    uses long_form to create a new permutation instance, finds the barcodes
    where the permuted genes of each category occur, and evalutes the neighbor
    distance and a new MaxMin score

    returns a new polars dataframe with columns [`barcode`, `MaxMin`]
    """
    if shuffle:
        long_form = long_form.with_columns(pl.col("gene").shuffle(seed=seed).alias("gene"))
    all_coords = positions.select("array_row", "array_col").to_numpy()
    for category, genes in SEAC.SEAC.items():
        gene_spots = long_form.filter(pl.col("gene").is_in(genes))["barcode"].unique() 

        # Build KDTree and evaluate Min distance for the seac category
        category_coords = (
            positions.filter(pl.col("barcode").is_in(gene_spots))
            .select("array_row", "array_col")
            .to_numpy()
        )
        tree = cKDTree(category_coords)
        min_distances, _ = tree.query(all_coords)
        positions = positions.with_columns(pl.lit(np.ceil(min_distances)).alias(category))

    positions = positions.with_columns(
        pl.max_horizontal(SEAC.SEAC.keys()).alias(newcolname)
    )
    return positions[["barcode", newcolname]]


def multiprocessing_neighbors(start_spots, pos, category, q):
    """ Not Tested. """
    print(category)
    cat_distance = blade.get_neighbors(
        pos, col=category, barcodes=start_spots, show_progress=False
    )
    q.put(cat_distance)


def multiprocessing_evaluate_permutation(long_form, positions):
    """
    NOT READY FOR PRIME TIME
    uses long_form to create a new permutation instance, finds the barcodes
    where the permuted genes of each category occur, and evalutes the neighbor
    distance and a new MaxMin score

    returns a new polars dataframe with columns [`barcode`, `MaxMin`]
    """
    q = mp.Queue()
    #   pool = mp.Pool(processes=4)
    long_form = long_form.with_columns(pl.col("gene").shuffle().alias("gene"))
    processes = []
    for category, genes in SEAC.SEAC.items():
        start_spots = long_form.filter(pl.col("gene").is_in(genes))["barcode"].unique()
        processes.append(
            mp.Process(
                target=multiprocessing_neighbors,
                args=[
                    list(start_spots),
                    positions.clone(),
                    category,
                    q,
                ],
            )
        )
    for p in processes:
        p.start()
    for p in processes:
        p.join()
    #   results = [q.get() for p in processes]
    return q


def adata2longform(sample, positions):
    barcodes = positions.filter(pl.col("blade").ne(1))["barcode"]
    row, col, count = scipy.sparse.find(sample.X)
    df = pl.DataFrame(
        {
            "barcode": sample.obs_names[row],
            "gene": sample.var_names[col],
            "nreads": count,
        }
    )

    long_form = df.select(pl.exclude("nreads").repeat_by("nreads").explode())
    long_form = long_form.filter(pl.col("barcode").is_in(barcodes))
    return long_form
