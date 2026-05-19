import multiprocessing as mp

import polars as pl
import scipy
from VisHD import SEAC, blade, local


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
            result = evaluate_permutation(long_form, positions, str(i))
        else:
            result = result.join(
                evaluate_permutation(long_form, positions, str(i)), on="barcode"
            )

    return result


def evaluate_permutation(long_form, positions, newcolname="MaxMin"):
    """
    uses long_form to create a new permutation instance, finds the barcodes
    where the permuted genes of each category occur, and evalutes the neighbor
    distance and a new MaxMin score

    returns a new polars dataframe with columns [`barcode`, `MaxMin`]
    """
    long_form = long_form.with_columns(pl.col("gene").shuffle().alias("gene"))
    for category, genes in SEAC.SEAC.items():
        start_spots = long_form.filter(pl.col("gene").is_in(genes))["barcode"].unique()
        positions = blade.get_neighbors(
            positions, col=category, barcodes=list(start_spots), show_progress=True
        )
    positions = positions.with_columns(
        pl.max_horizontal(SEAC.SEAC.keys()).alias(newcolname)
    )
    return positions[["barcode", newcolname]]


def multiprocessing_neighbors(start_spots, pos, category, q):
    print(category)
    cat_distance = blade.get_neighbors(
        pos, col=category, barcodes=start_spots, show_progress=False
    )
    q.put(cat_distance)


def multiprocessing_evaluate_permutation(long_form, positions):
    """
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
