import polars as pl
from tqdm import tqdm
from VisHD import blade

SEAC = {
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


def seac_distance(adata, pos):
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
