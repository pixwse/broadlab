import os

import geopandas as gpd
from shapely.geometry import Polygon, Point
from scipy import sparse
import scanpy as sc
import pandas as pd
import json

import numpy as np
import matplotlib.pyplot as plt
import cv2
import PIL

from utils import *

# Lab stuff for understanding the spatial relationship between ST
# data and H&E images
#

# ----------------------------------------------------------------------------
# Lab stuff

def lab_read_parquet():
    # Reads the tissue_positions parquet file

    input_path = '/media/erik/T9/run2_A/outs/binned_outputs/square_002um/spatial/tissue_positions.parquet'

    df = pd.read_parquet(input_path)
    df.keys()

    img_xs = df['pxl_col_in_fullres'].to_numpy()
    img_ys = df['pxl_row_in_fullres'].to_numpy()
    print(f'Image coords: x-range: {min(img_xs)}, {max(img_xs)}, y-range: {min(img_ys)}, {max(img_ys)}')

    print('Done')


def lab_overlay_bin_centers():
    # Show the H&E image with overlayed spot centers

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
    plt.savefig('_temp/CLL_2A_HE_spotcenter_overlay.pdf')

    while True:
        plt.pause(60)


def lab_resample_image_at_spots(input_dir: str, output_dir: str):
    # Extract an image ROI sampling from the exact spot locations, to help
    # getting an idea of the alignment
    # 
    # Conclusion: The 2A image looks quite well-aligned. Could probably be
    # refined slightly, but local optimization might suffice for that purpose.
    # 
    # Note: The 'highres' image that we start with isn't that high resolution
    # really. It's around 6k x 6k, but covering the entire slide. What we would
    # like to have that would be more convenient to work with is something like
    # a 26800x26800 image (8x the visium resolution) covering the exact same
    # area. Still quite high-res, but small enough to be loadable using opencv.
    # We should also do a better sampling (bilinear interpolation)

    # Load image
    print('Load image')
    image = cv2.imread(os.path.join(input_dir, 'tissue_hires_image.png'))
    image_height, image_width = image.shape[0:2]

    # Load spot centers
    print('Load spot centers')
    tp_df = pd.read_parquet(os.path.join(input_dir, 'tissue_positions.parquet'))
    spot_ys = tp_df['array_row'].to_numpy()
    spot_xs = tp_df['array_col'].to_numpy()
    img_xs = tp_df['pxl_col_in_fullres'].to_numpy()
    img_ys = tp_df['pxl_row_in_fullres'].to_numpy()
    print(f'Image coords: x-range: {min(img_xs)}, {max(img_xs)}, y-range: {min(img_ys)}, {max(img_ys)}')

    # Load scale factors
    print('Load scale factors')
    with open(os.path.join(input_dir, 'scalefactors_json.json'), 'rt') as file:
        scale_factors = json.load(file)
    factor = scale_factors['tissue_hires_scalef']

    spots_width = max(spot_xs) + 1
    spots_height = max(spot_ys) + 1 
    spot_image = np.zeros((spots_height, spots_width, 3), np.uint8)

    print('Process image')
    for sx,sy,ix,iy in zip(spot_xs, spot_ys, img_xs, img_ys):
        ix2 = clamp(int(ix * factor + 0.5), 0, image_width-1)
        iy2 = clamp(int(iy * factor + 0.5), 0, image_height-1)
        spot_image[sy, sx] = image[iy2, ix2]
    
    print('Saving image')
    cv2.imwrite(os.path.join(output_dir, 'image_resampled_at_spots.png'), spot_image)
    print('Done')


def save_total_counts_image(input_dir: str, output_dir: str):
    # Read the Visium HD using scanpy, sum all counts and save it directly
    # as a png image, where the pixel value is the total count (max 132, so
    # works fine)

    # Load Visium HD data and sum over columns (genes)
    adata = sc.read_10x_h5(os.path.join(input_dir, 'filtered_feature_bc_matrix.h5'))
    counts = adata.X.sum(1) 

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

    # Convert to u8 and save with different gains, for easier visualization (a bit hacky for now)
    save_ndarray_as_image(img, 1, os.path.join(output_dir, 'counts.gain01.png'))
    save_ndarray_as_image(img, 2, os.path.join(output_dir, 'counts.gain02.png'))
    save_ndarray_as_image(img, 4, os.path.join(output_dir, 'counts.gain04.png'))
    save_ndarray_as_image(img, 8, os.path.join(output_dir, 'counts.gain08.png'))
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


def lab_xy_normalize(input_dir: str):
    # Lab stuff for normalizing the 2um data by row/col means, to compensate
    # for the line-like artifacts. Conclusion: This normalization seems
    # reasonable, the data looks much better this way.

    img = cv2.imread(os.path.join(input_dir, 'counts.gain01.png'), cv2.IMREAD_GRAYSCALE)

    mean_rows = np.mean(img, axis=0, keepdims=True) # Mean over rows
    mean_cols = np.mean(img, axis=1, keepdims=True) # Mean over cols

    img_f = img / mean_rows
    img_f = img_f / mean_cols

    figure, axis = plt.subplots(2, 1)

    if True:
        # Save normalized image
        save_ndarray_as_image(img_f, 256, os.path.join(input_dir, 'counts_rowcolnorm.gain256.png'))
        save_ndarray_as_image(img_f, 512, os.path.join(input_dir, 'counts_rowcolnorm.gain512.png'))

    if False:
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
    # Lab code to understand the distribution of counts in the 2um data.
    # Conclusion: On locations where genes are found, it's typically just one
    # (99.38% of all times). In 0.59% of all times, there's two. 3 or more is
    # only in 0.03% of the cases. So when thinking about this, it makes sense
    # to think of it as binary data as an approximation.

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

    print('OK')


def lab_all_genes_two_rows():
    # Compare the full gene expression of row 16 and 17 of the 1D example,
    # where there is a notable difference in total counts.
    # Conclusion: It looks like just a scaling factor between the two. Other
    # variations are not statistically significant (can be just due to the
    # Poisson-ish distributed counts). To know if there are other significant
    # differences, we would need to e.g. aggregate over all low/high rows
    # and see if these differences average out or not. 

    input_file = '_temp/1D_counts_raw.png'
    img = cv2.imread(input_file, cv2.IMREAD_GRAYSCALE)

    mean_rows = np.mean(img, axis=0, keepdims=True) # Mean over rows
    mean_rows = mean_rows.flatten()

    print(f'Mean of row 16: {mean_rows[16]}')  # ~5
    print(f'Mean of row 17: {mean_rows[17]}')  # ~10

    # Load Visium HD data
    input_path = '/media/erik/T9/run2_A/outs/binned_outputs/square_002um/'
    adata = sc.read_10x_h5(os.path.join(input_path, 'filtered_feature_bc_matrix.h5'))

    genome_sum_row16 = np.zeros((1, 18085)) # TODO: Don't hard-code
    genome_sum_row17 = np.zeros((1, 18085))

    for ix, spot_name in enumerate(adata.obs_names):
        
        # For debugging
        # if ix == 10000:
        #     break

        x, y = spotname_to_xy(spot_name)
        if y == 16:
            genome_sum_row16 += adata.X[ix, :]

        if y == 17:
            genome_sum_row17 += adata.X[ix, :]

    genome_sum_row16 = np.asarray(genome_sum_row16).flatten()
    genome_sum_row17 = np.asarray(genome_sum_row17).flatten()

    check_sum16 = np.sum(genome_sum_row16).item()
    check_sum17 = np.sum(genome_sum_row17).item()
    print(f'total counts, row 16: {check_sum16}, row 17: {check_sum17}')

    # Look for highly expressed genes, print their counts + names
    for ix in range(len(genome_sum_row16)):
        count_row16 = int(genome_sum_row16[ix])
        if count_row16 > 10:
            count_row17 = int(genome_sum_row17[ix])
            name = adata.var_names[ix]
            print(f'Index {ix}: Name: {name}, count in row 16: {count_row16}, count in row 17: {count_row17}')

    figure, axis = plt.subplots(2, 1)

    axis[0].plot(genome_sum_row16)
    # axis[0].set_ylim([0.0, 20.0])
    axis[1].plot(genome_sum_row17)
    # axis[1].set_ylim([0.0, 20.0])

    figure.show()
    figure.savefig('_temp/1D_row_16_and_17.pdf')
    figure.waitforbuttonpress()


def lab_all_genes():
    # Show the full histogram of all genes in the entire sample. To get a feeling
    # of their count distribution
    # 
    # Conclusion: Only a few genes are very active. For most of them, we find
    # something like one mRNA per 10-100 cells. This means that we need more
    # significant pooling in order to measure anything meaningful. This could be
    # pooling over tissue area, or we could hope that lots of genes are co-expressed,
    # such that we can get significance by looking at PCA components instead.
    # 
    # Assuming that we can pool over cell type and that we find some cell type
    # that occupies 1% of the sample. Then this is still 100k spots, i.e. even
    # a gene that only active at a 0.01 count/spot intensity would still get
    # 1k counts. That should be enough to start to do statisticall relevant
    # things.

    # Load Visium HD data
    input_path = '/media/erik/T9/run2_A/outs/binned_outputs/square_002um/'
    adata = sc.read_10x_h5(os.path.join(input_path, 'filtered_feature_bc_matrix.h5'))

    genome_sum = adata.X.sum(0)  # Sum over spots
    genome_sum = np.asarray(genome_sum).flatten()

    # genome_sum *= 16 * 1/(3350*3350) # Rescale to roughly expected concentration per cell (if 1 cell is 4x4 pixels)

    # Look for highly expressed genes, print their counts + names
    # for ix in range(len(genome_sum_row16)):
    #     count_row16 = int(genome_sum_row16[ix])
    #     if count_row16 > 10:
    #         count_row17 = int(genome_sum_row17[ix])
    #         name = adata.var_names[ix]
    #         print(f'Index {ix}: Name: {name}, count in row 16: {count_row16}, count in row 17: {count_row17}')

    figure, axis = plt.subplots(2, 1)

    axis[0].plot(genome_sum)
    axis[1].plot(genome_sum)
    axis[1].set_ylim([0.0, 1e5])

    figure.show()
    figure.savefig('_temp/1D_all_genes.pdf')
    figure.waitforbuttonpress()


if __name__ == '__main__':
    # lab_read_parquet()
    # lab_read_visium()
    # lab_overlay_bin_centers()
    # lab_counts_as_image()

    # lab_nof_active_genes()
    # save_active_genes_image()

    # lab_all_genes_two_rows()
    # lab_all_genes()

    # Good for visualizing alignment
    save_total_counts_image('/media/erik/T9/run2_D/outs/binned_outputs/square_002um/', '_temp/CLL2D')
    lab_xy_normalize('_temp/CLL2D')
    lab_resample_image_at_spots('/media/erik/T9/run2_D/outs/binned_outputs/square_002um/spatial', '_temp/CLL2D')
