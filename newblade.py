import sys
from pathlib import Path

import polars as pl
from tqdm import tqdm

neighbor_set = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]


def neighbors(pos, col):
    max_dist = pos.select(pl.col(col).max()).item()
    rights = []
    for x, y in neighbor_set:
        if x == -1 and y == -1:
            rights.append(
                pos.filter(
                    pl.col(col).eq(max_dist),
                    pl.col("array_col").gt(0),
                    pl.col("array_row").gt(0),
                )
                .select(pl.col("array_col").add(x), pl.col("array_row").add(y))
                .cast(pl.UInt32)
            )
        elif x == -1:
            rights.append(
                pos.filter(pl.col(col).eq(max_dist), pl.col("array_col").gt(0))
                .select(pl.col("array_col").add(x), pl.col("array_row").add(y))
                .cast(pl.UInt32)
            )
        elif y == -1:
            rights.append(
                pos.filter(pl.col(col).eq(max_dist), pl.col("array_row").gt(0))
                .select(pl.col("array_col").add(x), pl.col("array_row").add(y))
                .cast(pl.UInt32)
            )
        else:
            rights.append(
                pos.filter(pl.col(col).eq(max_dist))
                .select(pl.col("array_col").add(x), pl.col("array_row").add(y))
                .cast(pl.UInt32)
            )
    right = pl.concat(rights).unique()
    left = pos.filter(pl.col(col).is_null())
    merge = left.join(right, on=["array_row", "array_col"], how="inner")
    return merge


def get_neighbors(r, col="dist", barcodes=[], max_iterations=3349):
    assert type(barcodes) is list
    if len(barcodes) == 0:
        r = r.with_columns(
            pl.when(pl.col("in_tissue") == 0).then(0).cast(pl.UInt16).alias(col),
        )
    else:
        r = r.with_columns(
            pl.when(
                pl.col("barcode")
                .is_in(pl.Series(barcodes).implode())
                .then(0)
                .cast(pl.UInt16)
                .alias(col)
            )
        )

    for _ in tqdm(range(max_iterations)):
        if len(r.filter(pl.col(col).is_null())) == 0:
            continue
        max_dist = r.select(pl.col(col).max()).item()
        last_bc = neighbors(r, col).select(pl.col("barcode")).to_series()
        r = r.with_columns(
            pl.when(pl.col("barcode").is_in(last_bc))
            .then(max_dist + 1)
            .otherwise(pl.col(col))
            .alias(col)
        )
    return r


def tissue_edge(pos, maxdist=2):
    pos = pos.with_columns(pl.lit(None).cast(pl.UInt16).alias("tissue_edge_distance"))

    pos = pos.with_columns(
        pl.when(
            (pl.col("array_row").eq(pl.col("array_row").min()))
            | (pl.col("array_row").eq(pl.col("array_row").max()))
            | (pl.col("array_col").eq(pl.col("array_col").min()))
            | (pl.col("array_col").eq(pl.col("array_col").max()))
        )
        .then(1)
        .otherwise("tissue_edge_distance")
        .alias("tissue_edge_distance"),
    )

    for i in range(2, maxdist + 1):
        last_bc = (
            neighbors(pos, "tissue_edge_distance").select(pl.col("barcode")).to_series()
        )
        pos = pos.with_columns(
            pl.when(pl.col("barcode").is_in(last_bc.implode()))
            .then(maxdist)
            .otherwise("tissue_edge_distance")
            .alias("tissue_edge_distance")
        )

    return pos


def blade(pos):
    pos = tissue_edge(pos, maxdist=2)
    pos = get_neighbors(pos, 

    return pos


def setup(pos):
    """
    setup, primarily used in debugging/development
    appends a columm called cdist to pos
    """
 #  result = pos.with_columns(
 #      pl.when(pl.col("in_tissue") == 0).then(0).cast(pl.UInt16).alias("cdist"),
#                       )
    pos = pl.read_parquet(
    "tissue_positions.parquet",
    columns=["barcode", "array_row", "array_col", "in_tissue"],
    )
    return result


if __name__ == "__main__":
    if len(sys.argv) == 1:
        # basically a debug/development situation
        pos = pl.read_parquet(
            "tissue_positions.parquet",
            columns=["barcode", "array_row", "array_col", "in_tissue"],
        )

    #    bc = pos.filter(pl.col("in_tissue") == 0).select(pl.col("barcode")).to_series()
    # r = setup(pos)
    else:
        infile = Path(sys.argv[1])
        if infile.is_dir:
            infile /= "spatial/tissue_positions.parquet"
        if not infile.exists():
            ValueError("Non-existent source file or directory")
        if len(sys.argv) > 2:
            outfile = infile.parent / "tissue_positions_bladed.parquet"
        else:
            outfile = Path(sys.argv[2])
