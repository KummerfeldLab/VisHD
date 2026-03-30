import polars as pl
from tqdm import tqdm
from VisHD import Adipose, blade, local
from pathlib import Path


SEAC260314 = {
    "SASP": ["IL1B", 
             "IL6", 
             "CXCL12",
             "CCL2", 
             "CXCL8",
             "TGFB1",
             "VEGFA"],
    "EMC": ["PLAUR", "PLAU", "PAPPA", "MMP14", "ADAMTS4", "COL6A1", "THBS1"],
    "Anti_Apoptotic": [
        "SERPINE1",
        "SERPINB2",
        "SERPINE2",
        "BCL2L1",
        "MCL1",
        "BCL2",
    ],
    "Cell_Cycle_Related": [
        "CDKN1A",
        "TP53",
        "GADD45A",
        "CDKN2A",
        "GADD45B",
    ],
}

old_SEAC = {
    "SASP": ["IL1B", "IL6", "CCL7", "CXCL8", "CCL2", "CXCL2", "CD55"],
    "EMC": ["PLAUR", "PLAU", "TNC", "VCAM1", "ADAMTS4", "TMEM217", "ICAM4"],
    "Anti_Apoptotic": [
        "SERPINE1",
        "SERPINB2",
        "BCL2L1",
        "GZMB",
        "RELT",
        "MCL1",
        "SERPINE2",
    ],
    "Cell_Cycle_Related": [
        "CDKN1A",
        "TP53",
        "GADD45A",
        "CDKN2A",
        "CCNB1",
        "CCNB2",
        "CDK1",
    ],
}


SEAC = SEAC260314

def make_maxMin(p, categories=SEAC, outdir=local.MaxMin_d, show_progress=True):
    """given path p, read in counts and positions and save a
    modified positions table to a dir"""
    print(p.name)
    subdir = "outs/binned_outputs/square_002um"
    adata, pos = Adipose.read_sample(p / subdir)
    r = blade.blade(pos)
    excluded_bc = (
        r.filter(pl.col("blade") == 1).select(pl.col("barcode")).to_series().to_list()
    )

    for k, v in categories.items():
        print('\t', k)
        goodgenes = [g for g in v if g in adata.var_names]
        print('\tGood:', goodgenes)
        badgenes = [g for g in v if g not in adata.var_names]
        print('\tBad:', badgenes)
        subtable = adata[:, goodgenes]
        bc = subtable[subtable.X.sum(axis=1) > 0].obs.index.to_list()
        bc = [x for x in bc if x not in excluded_bc]
        if len(bc) == 0:
            continue
        r = blade.get_neighbors(r, col=k, barcodes=bc, show_progress=show_progress)
    r = r.with_columns(pl.max_horizontal(categories.keys()).alias("MaxMin"))
    r.write_parquet(outdir / f"{p.name}_bladed.parquet")


def read_annotations(sdir, annotations_dir=local.annotations_dir):
    """ Find cell annotations that could be in one of several directories
        and have various names. """
    for anndir in annotations_dir.iterdir():
        if 'adults' not in anndir.name:
            continue
        for sample_ann in anndir.iterdir():
            if sdir.name in sample_ann.name:
                print('-->', sample_ann)
                return pl.read_csv(sample_ann / 'Graph-Based.csv')

def get_sample_and_annotations(sdir):
    """ Adds celltyping to the positions file """
    mapping_f = 'outs/barcode_mappings.parquet'
    mapping = pl.read_parquet(sdir / mapping_f)
    ann = read_annotations(sdir)
    
    bladed = pl.read_parquet(local.bladed_d / f'{sdir.name}_bladed.parquet')
    maxmin = bladed.with_columns(pl.lit('None').alias('Celltype'))
    ctype_col = ann.columns[-1]
    for ctype in ann[ctype_col].unique().to_list():
        cell_ids = ann.filter(pl.col(ctype_col).eq(ctype))['Barcode'].to_list()
        spots = mapping.filter(pl.col('cell_id').is_in(cell_ids))['square_002um'].to_list()
        maxmin = maxmin.with_columns(pl.when(pl.col('barcode').is_in(spots))
            .then(pl.lit(ctype))
            .otherwise(pl.col('Celltype'))
            .alias('Celltype'))
    maxmin = maxmin.with_columns(pl.lit(sdir.name).alias('Sample'))
    maxmin = maxmin.filter(pl.col('Celltype').ne('None'))
#   maxmin = maxmin[['barcode', 'MaxMin', 'Celltype', 'Sample']]
    return maxmin