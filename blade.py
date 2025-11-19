from pathlib import Path
import pandas as pd

VisiumHD_neighbor_set = [
    (-1, -1),
    (-1, 0),
    (-1, 1),
    (0, -1),
    (0, 1),
    (1, -1),
    (1, 0),
    (1, 1),
]

def get_neighbor_coords(pos, centroid_barcodes, neighbor_set):
    """ 
    Returns a subset of a positions table containing the neighbors of the centroids. 
    Neighbors are found by the pandas merge method, matching array_row and 
        array_col.
    
    parameters:
    ----------
        allpositions: a pandas DataFrame read from a Visium tissue_positions file 
            containing columns 'barcode', 'array_row', 'array_col'
            the numerical type of array_row and array_col must not be unsigned
        centroids: the barcodes of the points to find neighbors of
    """
    points = pos.loc[pos.barcode.isin(centroid_barcodes)]
    masks = []
    for diff_x, diff_y in neighbor_set:
        mask = pd.DataFrame()
        mask['array_row'] = points.array_row.add(diff_x)
        mask['array_col'] = points.array_col.add(diff_y)
        masks.append(mask.drop_duplicates())
    mask = pd.concat(masks)
    return mask.drop_duplicates()

    
def neighbors(pos, centroid_barcodes, neighbor_set=VisiumHD_neighbor_set):
    mask = get_neighbor_coords(pos, centroid_barcodes, neighbor_set)
    restricted_set = pos[~pos.barcode.isin(centroid_barcodes)]
    return(pd.merge(left=restricted_set, right=mask, how='inner', on=['array_row', 'array_col']))


def capture_edge_barcodes(pos, excluded_bcs=[]):
    """return a list of barcodes that are on the capture edge. Assumes that any spot on a max or min row or column is on the capture edge.
    parameters:
    ----------
        pos: a pandas DataFrame from a tissue_positions.parquet file or similar containing
            columns named `barcode`, `array_row`, `array_col`
        excluded_bcs: a list of barcodes to exclude from consideration. Allows for multiple calls to capture_edge_neighbors.
    """
    # filters out excluded_bcs
    pos = pos.loc[~pos.barcode.isin(excluded_bcs)]

    # find barcodes that are in a min or max column
    ce_barcodes = pos.loc[
        (pos.array_row == pos.array_row.min())
        | (pos.array_row == pos.array_row.max())
        | (pos.array_col == pos.array_col.min())
        | (pos.array_col == pos.array_col.max())
    ].barcode.to_list()

    # return a list of unique barcodes found in the above
    return list(set(ce_barcodes))
