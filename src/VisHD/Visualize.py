import os
os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = pow(2,41).__str__()
import json
from VisHD.SEAC import SEAC
import polars as pl
import matplotlib.pyplot as plt
from pathlib import Path
import cv2 as cv
#from sklearn.cluster import KMeans

ROWDIFF = 5
COLDIFF = 5
LIGHTBLUE = (173, 216, 230)
BLUE = (65, 105, 255)
RED = (255, 20, 147)

# BGR Colors for OpenCV
COLORS = {
    "hotspots": (209, 0, 181),  # Original magenta color
    "hotspots2": (0, 255, 0),
    "SASP": (0, 255, 0),  # LimeGreen for high visual pop
    "EMC": (255, 255, 0),  # Cyan gives crisp contrast on grayscale
    "Anti_Apoptotic": (0, 128, 255),  # Orange for a warm tone
    "Cell_Cycle_Related": (0, 230, 255),  # Yellow for maximum brightness
}

#colors = {
#    'SASP': (00, 209, 28),
#    'EMC': (0, 77, 209),
#    'Anti_Apoptotic': (209, 0, 181),
#    'Cell_Cycle_Related': (209, 132, 0)}


def make_seac_image(ann, pos, img):
    """ returns a copy of an tissue image annotated an image with 
    dots added where gene reads are round, colored by gene category """
    overlay = img.copy()
    for k,v in SEAC.items():
        subtable = ann[:, v]
        barcodes = subtable[subtable.X.sum(axis=1) > 0].obs.index.to_list()
        seac_positions = pos.filter(pl.col('barcode').is_in(barcodes))
        draw_spots(overlay, seac_positions, COLORS[k], sidelength=40)
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


def make_annotated_image(samplename, positions, spaceranger_dir):
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
    tif_d = Path(spaceranger_dir) / samplename / tif_ext
    for tif_f in (tif_d).iterdir():
        if 'tif' != tif_f.name[-3:]:
            continue
        break
    tif = cv.cvtColor(cv.cvtColor(cv.imread(tif_f), cv.COLOR_BGR2GRAY), cv.COLOR_GRAY2RGB)
    overlay = tif.copy()
    viz.draw_spots(overlay, positions, COLORS['magenta'], sidelength=5)

    # Determine coordinates for cropping
    xmin = int(positions.select(pl.col('pxl_row_in_fullres').min()).item())
    ymin = int(positions.select(pl.col('pxl_col_in_fullres').min()).item())
    xmax = int(positions.select(pl.col('pxl_row_in_fullres').max()).item())
    ymax = int(positions.select(pl.col('pxl_col_in_fullres').max()).item())

    # Join images
    final = cv.addWeighted(tif, .6, overlay, 0.4, 0)

    # Crop
    final = final[ymin-50 : ymax+50, xmin - 50 : xmax + 50]

    del tif, overlay
    return final


def make_multi_annotated_image(
    samplename,
    positions,
    barcodes,
    colors=COLORS,
    outdir=Path(),
    title="",
    sqr_pxls=5,
    cir_radius=10,
    cir_thickness=1,
    spaceranger_dir=None,
):
    """ returns an annotated image

    parameters:
    ----------
        samplename: the name of a sample.
        barcodes: a dictionary of where values are sets of barcodes to highlight
            on the image overlay
        colors: a dictionary of color values with the keys a subset of the keys in barcodes
    """
    
    if spaceranger_dir is None:
        raise ValueError("spaceranger_dir is required")

    # Read Positions file

    # Read in Tiff Image
    tif_ext = 'outs/spatial'
    tif_d = Path(spaceranger_dir) / samplename / tif_ext
    # The filename of the tif is not consistant, so we search for the first tif in the directory
    for tif_f in (tif_d).iterdir():
        if 'tif' != tif_f.name[-3:]:
            continue
        break

    # Greyscale the original image
    tif = cv.cvtColor(cv.cvtColor(cv.imread(tif_f), cv.COLOR_BGR2GRAY), cv.COLOR_GRAY2RGB)

    # Create the overlay with the annotations
    overlay = tif.copy()
    for cat, bc in barcodes.items():
        filtered = positions.filter(pl.col('barcode').is_in(bc))
        print(cat, colors[cat], end=' ')
        if 'hotspots'in cat:
            draw_spots(overlay, filtered, colors[cat], sidelength=sqr_pxls)
        else:
            for x,y in filtered.select("pxl_col_in_fullres", "pxl_row_in_fullres").iter_rows():
                cv.circle(overlay, center=(int(x),int(y)), radius=cir_radius, color=colors[cat], thickness=cir_thickness)


    # Determine coordinates for cropping
    xmin = int(positions.select(pl.col('pxl_row_in_fullres').min()).item())
    ymin = int(positions.select(pl.col('pxl_col_in_fullres').min()).item())
    xmax = int(positions.select(pl.col('pxl_row_in_fullres').max()).item())
    ymax = int(positions.select(pl.col('pxl_col_in_fullres').max()).item())

    # Join images
    final = cv.addWeighted(tif, .6, overlay, 0.4, 0)

    # Crop
    final = final[ymin-50 : ymax+50, xmin - 50 : xmax + 50]
    final = cv.copyMakeBorder(final, 340, 0, 0, 0, cv.BORDER_CONSTANT, value=(255, 255, 255))
    final = cv.putText(final, samplename + ' 1% Permutation Hotspots', (40, 280),
               cv.FONT_HERSHEY_TRIPLEX, 10, (0,0,0), 4)

    # Save
    cv.imwrite(outdir / f"{samplename}.png", final)
    del tif, overlay, final
