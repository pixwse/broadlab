
# Trying to visualize the CLL 2A example with the matching microscopy image
#
#
#

import os

import geopandas as gpd
from shapely.geometry import Polygon, Point
from scipy import sparse
import scanpy as sc
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib

import cv2


def lab_read_parquet():
    # Reads the tissue_positions parquet file
    #     

    input_path = '/media/erik/T9/run2_A/outs/binned_outputs/square_002um/spatial/tissue_positions.parquet'

    df = pd.read_parquet(input_path)
    df.keys()

    img_xs = df['pxl_col_in_fullres'].to_numpy()
    img_ys = df['pxl_row_in_fullres'].to_numpy()
    print(f'Image coords: x-range: {min(img_xs)}, {max(img_xs)}, y-range: {min(img_ys)}, {max(img_ys)}')

    print('Done')


def lab_overlay_bin_centers():

    plt.interactive(True)

    # Load image (one tile)
    tile_path = '/home/erik/encdata/pixwse/checkout/broadlab/_temp/cll2a_tiles_png'
    x0 = 50000
    y0 = 30000
    crop_width = 1000
    crop_height = 1000

    tile_file_name = f'tile_{x0}_{y0}.png'
    tile = cv2.imread(os.path.join(tile_path, tile_file_name))
    tile = cv2.cvtColor(tile, cv2.COLOR_BGR2RGB)

    tile = tile[0:crop_height, 0:crop_width]
    width = tile.shape[1]
    height = tile.shape[0]

    # Load spot centers
    tp_input_path = '/media/erik/T9/run2_A/outs/binned_outputs/square_002um/spatial/tissue_positions.parquet'
    tp_df = pd.read_parquet(tp_input_path)

    img_xs = tp_df['pxl_col_in_fullres'].to_numpy()
    img_ys = tp_df['pxl_row_in_fullres'].to_numpy()
    print(f'Image coords: x-range: {min(img_xs)}, {max(img_xs)}, y-range: {min(img_ys)}, {max(img_ys)}')

    # Transform and crop
    img_xs -= x0
    img_ys -= y0
    count = 0
    filtered_xs = []
    filtered_ys = []
    for x,y in zip(img_xs, img_ys):
        if 0 < x < width and 0 < y < height:
            filtered_xs.append(x)
            filtered_ys.append(y) 
            count += 1

    plt.imshow(tile)
    plt.scatter(filtered_xs, filtered_ys, marker='.', color='green', s=8)
    plt.show()

    while True:
        plt.pause(60)


def show_total_counts():
    # TODO: Show total read counts for each cell

    input_path = '/media/erik/T9/run2_A/outs/binned_outputs/square_002um'

    adata = sq.read.visium(
        input_path,
        counts_file='filtered_feature_bc_matrix.h5')





if __name__ == '__main__':
    # lab_read_parquet()
    # lab_read_visium()
    lab_overlay_bin_centers()
