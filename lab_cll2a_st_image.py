
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

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

import cv2
import PIL

# Helper functions

def spotname_to_xy(name: str) -> tuple[int, int]:
    # TODO: Error handling
    # TODO: Don't assume these spot names, should use the tissue_positions instead

    parts = name.split('_') # e.g. ['s', '002um', '00658', '01498-1']
   
    y = int(parts[2])
    parts2 = parts[3].split('-') # e.g. ['01498', '1']
    x = int(parts2[0])

    return (x,y)
    

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


def lab_counts_as_image():
    # Read the Visium HD using scanpy, sum all counts and visualize it as an
    # image, to understand how the datastructure works

    plt.interactive(True)
    input_path = '/media/erik/T9/run2_A/outs/binned_outputs/square_002um/'

    # Load Visium HD data and spatial coordinates of spots
    df_pos = pd.read_parquet(os.path.join(input_path, 'spatial/tissue_positions.parquet'))
    adata = sc.read_10x_h5(os.path.join(input_path, 'filtered_feature_bc_matrix.h5'))

    adata.X
    # adata.X is now a CSR matrix with shape 10732334x18085
    # Here, 18085 is the number of genes and 10732334 is the numer of spots. To
    # find the spatial location of each spot, we need to look at the
    # corresponding spot name or spot data in the parquet file.

    counts = adata.X.sum(1)  # Sum over columns

    # Put all counts in an image-like structure (slow, not the right way later
    # on, but to understand the structure)
    max_count = 0
    img = np.zeros((3350, 3350))
    for ix, c in enumerate(counts):

        if ix % 100000 == 0:
            print(f'Processed spots: {ix//100000}00k, max = {max_count}')

        max_count = max(c, max_count)

        # Find x,y coordinate
        spot_name = adata.obs_names[ix]
        x, y = spotname_to_xy(spot_name)

        # Basic normalization
        nc = c.item() / 80.0
        nc = min(nc, 1.0)

        # Set the pixel
        img[y, x] = nc

    print(f'max_count = {max_count}')

    # Convert to u8 and save
    img_u8 = (img * 255).astype(np.uint8)
    img_pil = PIL.Image.fromarray(img_u8)
    img_pil.save('_temp/temp.png')

    print('Done!')


def save_total_counts_image():
    # Read the Visium HD using scanpy, sum all counts and save it directly
    # as a png image, where the pixel value is the total count (max 132, so
    # works fine)

    input_path = '/media/erik/T9/run2_A/outs/binned_outputs/square_002um/'

    # Load Visium HD data
    adata = sc.read_10x_h5(os.path.join(input_path, 'filtered_feature_bc_matrix.h5'))

    counts = adata.X.sum(1)  # Sum over columns

    # Put all counts in an image-like structure (slow, not the right way later
    # on, but to understand the structure)
    max_count = 0
    img = np.zeros((3350, 3350), dtype=np.uint8)
    for ix, count in enumerate(counts):

        if ix % 100000 == 0:
            print(f'Processed spots: {ix//100000}00k, max = {max_count}')

        max_count = max(count, max_count)

        # Find x,y coordinate
        spot_name = adata.obs_names[ix]
        x, y = spotname_to_xy(spot_name)

        # Set the pixel
        img[y, x] = count

    print(f'max_count = {max_count}')

    # Convert to u8 and save
    img_pil = PIL.Image.fromarray(img)
    img_pil.save('_temp/temp_raw_u8.png')

    print('Done!')


def save_active_genes_image():
    # Read the Visium HD using scanpy, sum threshold at 1 (such that we only
    # measure the presence/absence of a specific gene at a specific location).
    # Then sum over all genes and save the result

    input_path = '/media/erik/T9/run2_A/outs/binned_outputs/square_002um/'

    # Load Visium HD data
    adata = sc.read_10x_h5(os.path.join(input_path, 'filtered_feature_bc_matrix.h5'))
    
    adata.X[adata.X >= 1] = 1
    counts = adata.X.sum(1)  # Sum over columns

    # Put all counts in an image-like structure (slow, not the right way later
    # on, but to understand the structure)
    max_count = 0
    img = np.zeros((3350, 3350), dtype=np.uint8)
    for ix, count in enumerate(counts):

        if ix % 100000 == 0:
            print(f'Processed spots: {ix//100000}00k, max = {max_count}')

        max_count = max(count, max_count)

        # Find x,y coordinate
        spot_name = adata.obs_names[ix]
        x, y = spotname_to_xy(spot_name)

        # Set the pixel
        img[y, x] = count

    print(f'max_count = {max_count}')

    # Convert to u8 and save
    img_pil = PIL.Image.fromarray(img)
    img_pil.save('_temp/temp_thresholded_raw.png')

    print('Done!')


def lab_xy_normalize():
    input_file = '_temp/2A_counts_raw.png'

    img = cv2.imread(input_file, cv2.IMREAD_GRAYSCALE)

    mean_rows = np.mean(img, axis=0, keepdims=True) # Mean over rows
    mean_cols = np.mean(img, axis=1, keepdims=True) # Mean over cols

    img_f = img / mean_rows
    img_f = img_f / mean_cols

    figure, axis = plt.subplots(2, 1)

    if False:
        # Save normalized image
        img_u8 = (img_f * 255).astype(np.uint8)
        img_pil = PIL.Image.fromarray(img_u8)
        img_pil.save('_temp/temp_rowcol_normalized.png')

    if True:
        # Plot row/col means

        mean_rows = mean_rows.flatten()
        mean_cols = mean_cols.flatten()

        mean_rows = np.sort(mean_rows)
        mean_cols = np.sort(mean_cols)

        axis[0].plot(mean_rows)
        axis[0].set_title('mean over rows')

        axis[1].plot(mean_cols)
        axis[1].set_title('mean over columns')

        figure.show()
        figure.waitforbuttonpress()


def lab_nof_active_genes():

    input_path = '/media/erik/T9/run2_A/outs/binned_outputs/square_002um/'

    # Load Visium HD data
    adata = sc.read_10x_h5(os.path.join(input_path, 'filtered_feature_bc_matrix.h5'))

    counts = adata.X.sum(1)  # Sum over columns
    hist_counts, hist_bins = np.histogram(adata.X.data, bins=np.arange(30))

    total_sum = np.sum(hist_counts)
    hist_rates = hist_counts / total_sum

    nof_1s = hist_rates[1]
    nof_2s = hist_rates[2]
    nof_3s = hist_rates[3]
    print(f'Rate of 1s: {nof_1s*100:.02f}%')
    print(f'Rate of 2s: {nof_2s*100:.02f}%')
    print(f'Rate of 3s: {nof_3s*100:.02f}%')
    # Conclusion: On locations where genes are found, it's typically just one
    # (99.38% of all times). In 0.59% of all times, there's two. 

    print('OK')


if __name__ == '__main__':
    # lab_read_parquet()
    # lab_read_visium()
    # lab_overlay_bin_centers()
    # lab_counts_as_image()

    # save_total_counts_image()
    # lab_xy_normalize()

    # lab_nof_active_genes()
    save_active_genes_image()
