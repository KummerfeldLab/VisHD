import polars as pl
import scipy
from VisHD import SEAC, blade, local


def read_basedata(sampledir):
    sample = SEAC.read_h5(sampledir / "outs/binned_outputs/square_002um")
    positions = pl.read_parquet(local.bladed_d / f"{sampledir.name}_bladed.parquet")
    return sample, positions


def permute_SEAC_sample(sampledir):
    sample, positions = read_basedata(sampledir)
    long_form = get_longform(sample, positions)
    return long_form


def evaluate_permutation(long_form, positions):
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
        pl.max_horizontal(SEAC.SEAC.keys()).alias("MaxMin")
    )
    return positions[["barcode", "MaxMin"]]


def get_longform(sample, positions):
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
