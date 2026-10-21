import scipy
import polars as pl
import sys
from VisHD.SEAC import SEAC

genelists = SEAC
outdir = sys.argv[1]
samples = sys.argv[2:]
print(outdir)
print(samples)


def make_gene_count_table(samples, gene_dict=genelists, outdir=None):
    """Reads each sample in `samples`, and returns a sample counting the number
    of genes in the values of each key of gene_dict as well as the number of non-category genes"""
    table = None
    for sample in samples:
        sampledir = local.spaceranger_outputs_dir / sample
        adata, pos = read_basedata(sampledir)

        #
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

        gene_categories = list(gene_dict.keys())
        genecat = gene_categories[0]
        when_expr = pl.coalesce(
            pl.when(pl.col("gene").is_in(v)).then(pl.lit(k))
            for k, v in gene_dict.items()
        )

        category_counts = long_form.select(
            when_expr.fill_null(pl.lit("Other")).alias("category")
        )["category"].value_counts()
        category_counts = category_counts.select(
            "category", pl.col("count").alias(x.name)
        )

        if table is None:
            table = category_counts
        else:
            table = table.join(category_counts, on="category")

    if outdir is None:
        return table
    else:
        table.write_parquet(f"{outdir}/genecat_table.py")
