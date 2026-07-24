import os
os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = pow(2,41).__str__()
from VisHD.SEAC import SEAC
from VisHD import local
import polars as pl
import matplotlib.pyplot as plt
from pathlib import Path
import cv2 as cv
import VisHD.Visualize as viz
#from sklearn.cluster import KMeans

# overlay = np.zeros_like(gs_img, np.uint8)
ROWDIFF = 5
COLDIFF = 5
LIGHTBLUE = (173, 216, 230)
BLUE = (65, 105, 255)
RED = (255, 20, 147)

# BGR Colors for OpenCV
COLORS = {
    "magenta": (209, 0, 181),  # Your original color
    "lime_green": (0, 255, 0),  # High visual pop
    "cyan": (255, 255, 0),  # Crisp contrast on grayscale
    "orange": (0, 128, 255),  # Distinct warm tone
    "yellow": (0, 230, 255),  # Maximum brightness
}

colors = {
    'SASP': (00, 209, 28),
    'EMC': (0, 77, 209),
    'Anti_Apoptotic': (209, 0, 181),
    'Cell_Cycle_Related': (209, 132, 0)}


def make_seac_image(ann, pos, img):
    """ returns a copy of an tissue image annotated an image with 
    dots added where gene reads are round, colored by gene category """
    overlay = img.copy()
    for k,v in SEAC.items():
        subtable = ann[:, v]
        barcodes = subtable[subtable.X.sum(axis=1) > 0].obs.index.to_list()
        seac_positions = pos.filter(pl.col('barcode').is_in(barcodes))
        viz.draw_spots(overlay, seac_positions, colors[k], sidelength=40)
    return overlay
#   return cv.addWeighted(img, 0.4, overlay, 0.55, 0.05)


def draw_spots(overlay, pos, color, sidelength=5):
    """ annotates IN PLACE an opencv2 image by drawing squares of color `color`
    over each barcode in pos.
    `pos` is a tissue_positions.parquet file limited to only the rows 
    which are to be colored. """
    assert len(color) == 3

    filtered = pos.with_columns(
        pl.col("pxl_row_in_fullres").add(sidelength / 2).alias("Y2"),
        pl.col("pxl_col_in_fullres").add(sidelength / 2).alias("X2"),
    )
    filtered = filtered.with_columns(
        pl.col("pxl_row_in_fullres").sub(sidelength / 2).alias("Y1"),
        pl.col("pxl_col_in_fullres").sub(sidelength / 2).alias("X1"),
    )
    filtered = filtered.select(
        pl.concat_list(["X1", "Y1"]).alias("p1"),
        pl.concat_list(["X2", "Y2"]).alias("p2"),
    )
    for (x1, y1), (x2, y2) in filtered.iter_rows():
        cv.rectangle(overlay, (int(x1), int(y1)), (int(x2), int(y2)), color, cv.FILLED)


def make_annotated_image(samplename, positions):
    """ returns an annotated image

    parameters:
    ----------
        samplename: the name of a sample.
        positions: the rows of a polars dataframe read from a 
            tissue_positions.parquet file that are to be highlighted

    Code example to save a file:
    ------
       ```
        import matplotlib.pyplot as plt
        plt.axis('off')
        plt.title(name + 'Region <= .025 ')
        plt.imshow(final)
        plt.savefig(outdir / f"{name}.png")
        ```
    """
 
    tif_ext = 'outs/spatial'
    tif_d = local.spaceranger_outputs_dir / samplename / tif_ext
    for tif_f in (tif_d).iterdir():
        if 'tif' != tif_f.name[-3:]:
            continue
        break
    tif = cv.cvtColor(cv.cvtColor(cv.imread(tif_f), cv.COLOR_BGR2GRAY), cv.COLOR_GRAY2RGB)
    overlay = tif.copy()
    viz.draw_spots(overlay, positions, (209, 0, 181), 40)
    xmin = int(positions.select(pl.col('pxl_row_in_fullres').min()).item())
    ymin = int(positions.select(pl.col('pxl_col_in_fullres').min()).item())
    xmax = int(positions.select(pl.col('pxl_row_in_fullres').max()).item())
    ymax = int(positions.select(pl.col('pxl_col_in_fullres').max()).item())
    final = cv.addWeighted(tif, .6, overlay, 0.4, 0)
    final = final[ymin-50 : ymax+50, xmin - 50 : xmax + 50]

    del tif, overlay
    return final
