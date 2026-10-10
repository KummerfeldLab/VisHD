from pathlib import Path

import numpy as np
import polars as pl
import scipy
from scipy.spatial import KDTree
from VisHD import SEAC, local


def read_basedata(sampledir):
    """
    Given the name of a sample return the original count matrix as anndata,
    and the bladed and MaxMinned positions file as a polars dataframe
    """
    sample = SEAC.read_h5(sampledir / "outs/binned_outputs/square_002um")
    positions = pl.read_parquet(local.bladed_d / f"{sampledir.name}_bladed.parquet")
    return sample, positions


def adata2longform(adata, positions):
    barcodes = positions.filter(pl.col("blade").ne(1))["barcode"]
    row, col, count = scipy.sparse.find(adata.X)
    df = pl.DataFrame(
        {
            "barcode": adata.obs_names[row],
            "gene": adata.var_names[col],
            "nreads": count,
        }
    )

    long_form = df.select(pl.exclude("nreads").repeat_by("nreads").explode())
    long_form = long_form.filter(pl.col("barcode").is_in(barcodes))
    return long_form


def make_longforms(outdir, indir=local.spaceranger_outputs_dir):
    """
    Reads the spaceranger output of each sample in a directory,
    and writes a column of barcodes to a file of the same name as the sample.
    The barcodes are from the longform of the sample.
    """
    for sample in indir.iterdir():
        adata, pos = read_basedata(sample)
        lf = adata2longform(adata, pos)
        lf.select("barcode").write_parquet(f"{outdir}/{sample.name}.parquet")


def make_gene_count_table(samples, gene_dict=SEAC.SEAC):
    """Reads each sample in `samples`, and returns a sample counting the number
    of genes in the values of each key of gene_dict as well as the number of non-category genes"""
    table = None
    for sample in samples:
        sampledir = local.spaceranger_outputs_dir / sample
        adata, pos = read_basedata(sampledir)

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

#       gene_categories = list(gene_dict.keys())
#       genecat = gene_categories[0]

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


def evaluate_permutation(longform, permutedgenes, positions):
    """Provide the MinMax values for each spot give a permuted longform,
    in the same order as positions
    """

    tempdf = longform.with_columns(pl.lit(permutedgenes).alias("category"))
    all_coords = positions.select("array_row", "array_col").to_numpy()
    categories = [x for x in tempdf["category"].unique() if x != "Other"]
    for category in categories:
        if category == "Other":
            continue
        category_spots = tempdf.filter(pl.col("category").eq(category))[
            "barcode"
        ].to_list()
        category_coords = (
            positions.filter(pl.col("barcode").is_in(category_spots))
            .select("array_row", "array_col")
            .to_numpy()
        )
        tree = KDTree(category_coords)
        distances, _ = tree.query(all_coords)
        positions = positions.with_columns(
            pl.lit(np.ceil(distances).astype(int)).alias(category)
        )
    maxmin = positions.select(categories).max_horizontal()
    return maxmin


def permute_sample(samplename, start=0, end=100):
    """Provides MaxMin evaluations of permutations for one sample of the genes of all samples"""

    lf_path = Path.home() / f"All Projects/Adipose/longforms/{samplename}.parquet"
    longform = pl.read_parquet(lf_path)
    t = pl.read_csv(
        Path.home() / "All Projects/Adipose/260605 BigSim/SEAC_category_counts.csv"
    )
    null_perm = []
    for cat, n in zip(t["category"], t["Total"]):
        null_perm += [cat] * n

    indices = dict()
    cumsum = 0
    for s in t.columns[1:-1]:
        begin = cumsum
        cumsum = begin + t[s].sum()
        indices[s] = (begin, cumsum)

    start_idx, end_idx = indices[samplename]

    positions = pl.read_parquet(local.bladed_d / f"{samplename}_bladed.parquet")
    results = positions.select("barcode")

    for i in range(start, end):
        print(i, flush=True)
        permutation = np.random.default_rng(seed=i).permuted(null_perm)[
            start_idx:end_idx
        ]
        results = results.with_columns(
            pl.lit(evaluate_permutation(longform, permutation, positions)).alias(str(i))
        )

    return results
