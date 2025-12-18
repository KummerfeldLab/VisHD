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


def get_neighbors(r, col="dist", barcodes=[], max_iterations=3349, show_progress=True):
    assert type(barcodes) is list
    if len(barcodes) == 0:
        r = r.with_columns(
            pl.when(pl.col("in_tissue") == 0).then(0).cast(pl.UInt16).alias(col),
        )
    else:
        r = r.with_columns(
            pl.when(pl.col("barcode").is_in(pl.Series(barcodes).implode()))
            .then(0)
            .cast(pl.UInt16)
            .alias(col)
        )

    for _ in tqdm(range(max_iterations), disable=not show_progress):
        if len(r.filter(pl.col(col).is_null())) == 0:
            continue
        max_dist = r.select(pl.col(col).max()).item()
        last_bc = neighbors(r, col).select(pl.col("barcode")).to_series()
        r = r.with_columns(
            pl.when(pl.col("barcode").is_in(last_bc.implode()))
            .then(max_dist + 1)
            .otherwise(pl.col(col))
            .alias(col)
        )
    return r


def capture_edge(pos, maxdist=2):
    pos = pos.with_columns(pl.lit(None).cast(pl.UInt16).alias("capture_edge_distance"))

    pos = pos.with_columns(
        pl.when(
            (pl.col("array_row").eq(pl.col("array_row").min()))
            | (pl.col("array_row").eq(pl.col("array_row").max()))
            | (pl.col("array_col").eq(pl.col("array_col").min()))
            | (pl.col("array_col").eq(pl.col("array_col").max()))
        )
        .then(1)
        .otherwise("capture_edge_distance")
        .alias("capture_edge_distance"),
    )

    for i in range(2, maxdist + 1):
        last_bc = (
            neighbors(pos, "capture_edge_distance")
            .select(pl.col("barcode"))
            .to_series()
        )
        pos = pos.with_columns(
            pl.when(pl.col("barcode").is_in(last_bc.implode()))
            .then(maxdist)
            .otherwise("capture_edge_distance")
            .alias("capture_edge_distance")
        )

    return pos


def blade(pos, show_progress=False):
    pos = capture_edge(pos, maxdist=2)
    pos = get_neighbors(
        pos, col="tissue_edge_distance", max_iterations=1, show_progress=show_progress
    )
    pos = pos.with_columns(
        pl.when(
            pl.col("capture_edge_distance").gt(0) | pl.col("tissue_edge_distance").gt(0)
        )
        .then(1)
        .otherwise(0)
        .alias("blade")
    )

    return pos


def setup():
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
    return pos


if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("""
              blade.py
              Usage:
                python blade.py [infile | directory] <outfile>
                Takes either a tissue_positions.parquet file or a directory containing one.
                if no outfile is specified, will create/overwrite a tissue_positions_bladed.parquet
                file in the given directory or same directory as the input file.


              Currently only working on visium hd files
              """)
    else:
        infile = Path(sys.argv[1])
        if infile.is_dir():
            infile /= "tissue_positions.parquet"
        if not infile.exists():
            ValueError("Non-existent source file or directory")
        if len(sys.argv) > 2:
            outfile = Path(sys.argv[2])
        else:
            outfile = infile.parent / "tissue_positions_bladed.parquet"
        print(f"infile: {infile}\noutfile:{outfile}")
        pos = pl.read_parquet(infile)
        bladed = blade(pos).drop(["capture_edge_distance", "tissue_edge_distance"])
        bladed.write_parquet(outfile)
