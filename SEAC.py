import warnings
import polars as pl
from tqdm import tqdm
from VisHD import Adipose, blade
import scanpy as sc

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

def read_sample(sampledir, set_barcode_as_index=False):
    """reads the h5 file as anndata and the positions_file as a pandas DataFrame and returns both"""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        adata = sc.read_10x_h5(sampledir / "filtered_feature_bc_matrix.h5")
        adata.var_names_make_unique()

    pos = pl.read_parquet(sampledir / "spatial/tissue_positions.parquet")
    #   pos.array_row = pos.array_row.astype(int)
    #   pos.array_col = pos.array_col.astype(int)

    #   if set_barcode_as_index:
    #       pos.set_index("barcode", inplace=True)
    return adata, pos

def make_bladed_maxMin_positionsFiles(p, outdir, show_progress=True):
    """given path p, read in counts and positions and save a
    modified positions table to a dir"""
    print(p.name)
    subdir = "outs/binned_outputs/square_002um"
#   adata, pos = Adipose.read_sample(p / subdir)
    adata, pos = read_sample(p / subdir)
    sc.experimental.pp.normalize_pearson_residuals(adata, inplace=True, copy=False, check_values=False)
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
        bc = [x for x in bc if x not in excluded_bc]
        if len(bc) == 0:
            continue
        r = blade.get_neighbors(r, col=k, barcodes=bc, show_progress=show_progress)
    r = r.with_columns(pl.max_horizontal(SEAC.keys()).alias("MaxMin"))
    r.write_parquet(outdir / f"{p.name}_bladed.parquet")


def seac_distance(adata, pos):
    """NOT TESTED"""
    sasp_adata = adata[:, SEAC["SASP"]]
    emc_adata = adata[:, SEAC["EMC"]]
    antiapop_adata = adata[:, SEAC["Anti_Apoptotic"]]
    ccr_adata = adata[:, SEAC["Cell_Cycle_Related"]]

    sasp_bc = sasp_adata[sasp_adata.X.sum(axis=1) > 0].obs.index.to_list()
    emc_bc = sasp_adata[emc_adata.X.sum(axis=1) > 0].obs.index.to_list()
    antiapop_bc = antiapop_adata[antiapop_adata.X.sum(axis=1) > 0].obs.index.to_list()
    ccr_bc = ccr_adata[ccr_adata.X.sum(axis=1) > 0].obs.index.to_list()

    r = blade.blade(pos)
    r = r.filter(pl.col("blade") == 0)
    r = r.with_columns(
        pl.when(pl.col("barcode").is_in(pl.Series(sasp_bc).implode()))
        .then(0)
        .cast(pl.UInt16)
        .alias("SASP"),
        pl.when(pl.col("barcode").is_in(pl.Series(emc_bc).implode()))
        .then(0)
        .cast(pl.UInt16)
        .alias("EMC"),
        pl.when(pl.col("barcode").is_in(pl.Series(antiapop_bc).implode()))
        .then(0)
        .cast(pl.UInt16)
        .alias("Anti_Apoptotic"),
        pl.when(pl.col("barcode").is_in(pl.Series(ccr_bc).implode()))
        .then(0)
        .cast(pl.UInt16)
        .alias("Cell_Cycle_Related"),
    )

    for _ in tqdm(range(8000)):
        null_sum = len(
            r.filter(
                pl.all_horizontal(
                    pl.col(
                        ["SASP", "EMC", "Anti_Apoptotic", "Cell_Cycle_Related"]
                    ).is_null()
                )
            )
        )
        if null_sum == 0:
            continue

        sasp_max = r.select(pl.col("SASP").max()).item()
        emc_max = r.select(pl.col("EMC").max()).item()
        antiapop_max = r.select(pl.col("Anti_Apoptotic").max()).item()
        ccr_max = r.select(pl.col("Cell_Cycle_Related").max()).item()

        r = r.with_columns(
            pl.when(
                pl.col("barcode").is_in(
                    blade.neighbors(r, "SASP")
                    .select(pl.col("barcode"))
                    .to_series()
                    .implode()
                )
            )
            .then(sasp_max + 1)
            .otherwise(pl.col("SASP"))
            .alias("SASP"),
            pl.when(
                pl.col("barcode").is_in(
                    blade.neighbors(r, "EMC")
                    .select(pl.col("barcode"))
                    .to_series()
                    .implode()
                )
            )
            .then(emc_max + 1)
            .otherwise(pl.col("EMC"))
            .alias("EMC"),
            pl.when(
                pl.col("barcode").is_in(
                    blade.neighbors(r, "Anti_Apoptotic")
                    .select(pl.col("barcode"))
                    .to_series()
                    .implode()
                )
            )
            .then(antiapop_max + 1)
            .otherwise(pl.col("Anti_Apoptotic"))
            .alias("Anti_Apoptotic"),
            pl.when(
                pl.col("barcode").is_in(
                    blade.neighbors(r, "Cell_Cycle_Related")
                    .select(pl.col("barcode"))
                    .to_series()
                    .implode()
                )
            )
            .then(ccr_max + 1)
            .otherwise(pl.col("Cell_Cycle_Related"))
            .alias("Cell_Cycle_Related"),
        )
    return r
