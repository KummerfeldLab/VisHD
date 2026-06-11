import warnings

import polars as pl
import scanpy as sc
from VisHD import blade

SEAC260314 = {
    "SASP": ["IL1B", "IL6", "CXCL12", "CCL2", "CXCL8", "TGFB1", "VEGFA"],
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


def read_h5(sampledir):
    """reads the h5 file as anndata, handling warnings and duplicate var names"""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        adata = sc.read_10x_h5(sampledir / "filtered_feature_bc_matrix.h5")
        adata.var_names_make_unique()
    return adata

    """and the positions_file as a pandas DataFrame and returns both"""


def read_positions(sampledir):
    pos = pl.read_parquet(sampledir / "spatial/tissue_positions.parquet")
    #   pos.array_row = pos.array_row.astype(int)
    #   pos.array_col = pos.array_col.astype(int)

    #   if set_barcode_as_index:
    #       pos.set_index("barcode", inplace=True)
    return pos


def read_sample(sampledir):  # set_barcode_as_index=False):
    adata = read_h5(sampledir)
    pos = read_positions(sampledir)
    return adata, pos


def make_bladed_maxMin_positionsFiles(p, outdir, show_progress=True):
    """given path p, read in counts and positions and save a
    modified positions table to a dir"""
    print(p.name)
    subdir = "outs/binned_outputs/square_002um"
    adata, pos = read_sample(p / subdir)
    r = blade.blade(pos)
    excluded_bc = (
        r.filter(pl.col("blade") == 1).select(pl.col("barcode")).to_series().to_list()
    )

    for k, v in SEAC.items():
        print("\t", k)
        goodgenes = [g for g in v if g in adata.var_names]
        print("\tGood:", goodgenes)
        badgenes = [g for g in v if g not in adata.var_names]
        print("\tBad:", badgenes)
        subtable = adata[:, goodgenes]
        bc = subtable[subtable.X.sum(axis=1) > 0].obs.index.to_list()
#       bc = [x for x in bc if x not in excluded_bc]
        bc = list(set(bc).difference(excluded_bc))
        if len(bc) == 0:
            continue
        r = blade.get_neighbors(r, col=k, barcodes=bc, show_progress=show_progress)
    r = r.with_columns(pl.max_horizontal(SEAC.keys()).alias("MaxMin"))
    r.write_parquet(outdir / f"{p.name}_bladed.parquet")
