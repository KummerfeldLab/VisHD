import gc
import sys
import time
from pathlib import Path

import polars as pl
from VisHD import local


def evalPermutations():
    for sdir in local.bladed_d.iterdir():
        gc.collect()
        blade = pl.read_parquet(sdir).filter(pl.col("blade").eq(0))
        sample_name = sdir.name.split("_bladed")[0]
        print(f"looking for {sample_name}")
#       start_time = time.time()
        results = None
        total_perms = 0
        for fname in Path("PermutationOutput").iterdir():
            if sample_name not in fname.name:
                continue
            print(fname.name)
            df = pl.read_parquet(fname) # .filter(
#               pl.col("barcode").is_in(blade["barcode"].implode())
#           )
            df = df.select(
                pl.col("barcode"),
                pl.concat_list(pl.exclude("barcode")).list.sort().alias("permutation"),
            )
            df = df.join(blade[['barcode', 'MaxMin']], on='barcode')
            df = df.with_columns(pl.lit(0).alias("pval"))
            listlen = df.select(pl.col("permutation").list.len().min()).item()
            total_perms += listlen
            for i in range(1, listlen):
                df = df.with_columns(
                    pl.when(pl.col("permutation").list.get(i).lt(pl.col("MaxMin")))
                    .then(pl.lit(i))
                    .otherwise(pl.col("pval"))
                    .alias("pval")
                )
            if results is None:
                results = df[['barcode','pval']]
            else:
                results = results.join(df[['barcode', 'pval']], on='barcode')
                results = results.with_columns(
                    pl.sum_horizontal(pl.exclude('barcode')).alias('total')
                )
                results = results[['barcode', 'total']]
        results = results.with_columns(
            pl.lit(total_perms).alias('n_permutations')
        )
        results.write_parquet(Path('Permuted_pvals') / f"{sample_name}.parquet")
#       print("Duration:", round(time.time() - start_time, 2), "seconds\n")


def permutation_summary_stats(sample_name):
    """ We want to know the .025 percentile threshold of each cell
    in each sample """
    basedir = Path.home() / "All Projects/Adipose/MaxMinSEAC"
    results_dir = basedir / "260520_PermutationSummary"
    source_dir = basedir / "PermutationOutput"

    result = None
    total_perms = 0
#   for fname in Path("PermutationOutput").iterdir():
    for fname in source_dir.iterdir():
        if sample_name not in fname.name:
            continue
        print(fname.name)
        df = pl.read_parquet(fname)
        total_perms += (df.shape[1] - 1)
        df = df.select(
            pl.col("barcode"),
            pl.concat_list(pl.exclude("barcode")).list.sort().list.head(25).alias('low'),
        )
        if result is None:
            result = df
        else:
            result = result.join(df, on='barcode')
            result = result.select(
                pl.col('barcode'),
                pl.col('low').list.concat('low_right').list.sort().list.head(25).alias('low')
            )
    nperms = (total_perms * 25) // 1000
    result = result.with_columns(pl.col('low').list.head(nperms).list.max().alias('low'))
    result.write_parquet(results_dir / f"{sample_name}.parquet")


if __name__ == "__main__":
    print("I'm main")
    snames = [x.name.split("_bladed")[0] for x in local.bladed_d.iterdir()]
    sample_name = snames[int(sys.argv[-1])]
    print(sample_name)
    permutation_summary_stats(sample_name)