import warnings

import numpy as np
import polars as pl
import scanpy as sc
from scipy.spatial import KDTree

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


def make_bladed_maxMin_positionsFiles(p, outdir):
    """read original positions file  and anndata from path p,
    indicate rows flagged by BLADE,
    and calculate SEAC distances,
    and save modified positions file to  outdir"""

    # Read in Data
    print(p.name)
    subdir = "outs/binned_outputs/square_002um"
    adata, pos = read_sample(p / subdir)

    # Add a `blade` column to the positions dataframe
    # And identify blade target spots
    all_coords = pos.select("array_row", "array_col").to_numpy()
    tree = KDTree(all_coords)
    pos = pos.with_columns(pl.lit(0).alias("blade"))
    blade_start_barcodes = pos.filter(
        (pl.col("array_row").eq(pl.col("array_row").min()))
        | (pl.col("array_row").eq(pl.col("array_row").max()))
        | (pl.col("array_col").eq(pl.col("array_col").min()))
        | (pl.col("array_col").eq(pl.col("array_col").max()))
        | (pl.col("in_tissue").eq(0))
    )["barcode"]
    blade_coords = (
        pos.filter(pl.col("barcode").is_in(blade_start_barcodes.implode()))
        .select("pxl_row_in_fullres", "pxl_col_in_fullres")
        .to_numpy()
    )
    query_output = tree.query_ball_point(blade_coords, r=1, workers=-1)
    blade_idxs = set()
    for spot_neighbors in query_output:
        blade_idxs.update(spot_neighbors)
    blade_barcodes = pos["barcode"].to_numpy()[list(blade_idxs)]
    pos.with_columns(
        pl.when(pl.col("barcode").is_in(blade_barcodes))
        .then(1)
        .otherwise(0)
        .alias("blade")
    )

    # Evaluate SEAC MaxMin Distance
    for category, genes in SEAC.items():
        # Ensure SEAC genes appear in adata
        print("\t", category)
#       goodgenes = [g for g in genes if g in adata.var_names]
#       print("\tGood:", goodgenes)
#       badgenes = [g for g in genes if g not in adata.var_names]
#       print("\tBad:", badgenes)
#       subtable = adata[:, goodgenes]
        subtable = adata[:, genes]
        bc = subtable[subtable.X.sum(axis=1) > 0].obs.index.to_list()
        if len(bc) == 0:
            print(category, 'has no spots')
            continue

        # Build KDTree and evaluate Min distance for the seac category
        category_coords = pos.filter(pl.col("blade").eq(0),
                                     pl.col("barcode").is_in(bc)
                                    ).select("array_row", "array_col").to_numpy()
        tree = KDTree(category_coords)
        min_distances, _ = tree.query(all_coords)
        pos = pos.with_columns(pl.lit(np.ceil(min_distances)).alias(category))


    maxmin = pos.select(pl.col(SEAC.SEAC.keys())).max_horizontal()
    pos = pos.with_columns(pl.lit(maxmin).alias('MaxMin'))
#   pos = pos.with_columns(pl.max_horizontal(SEAC.keys()).alias("MaxMin"))
    pos.write_parquet(outdir / f"{p.name}_bladed.parquet")
