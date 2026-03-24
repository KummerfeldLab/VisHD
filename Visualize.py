import cv2 as cv
import polars as pl

# overlay = np.zeros_like(gs_img, np.uint8)
ROWDIFF = 5
COLDIFF = 5
LIGHTBLUE = (173, 216, 230)
BLUE = (65, 105, 255)
RED = (255, 20, 147)


def draw_spots(overlay, pos, color, sidelength=5):
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

