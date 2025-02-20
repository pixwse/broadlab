import os

import geopandas as gpd
from shapely.geometry import Polygon, Point
from scipy import sparse
import scanpy as sc
import pandas as pd
import json
import anndata

import numpy as np
import matplotlib.pyplot as plt
import cv2
import PIL
from scipy.ndimage import gaussian_filter1d
from scipy.signal import medfilt

from utils import *


# Lab code for normalizing each row/column separately individually for each
# gene. Looking at some highly expressed genes first, then for others


def lab_save_gene_image(
    input_dir: str, 
    gene_name: str, 
    output_prefix: str,
    load_h5ad = False):

    if load_h5ad:
        adata = anndata.io.read_h5ad(input_dir)
    else:
        adata = sc.read_10x_h5(os.path.join(input_dir, 'filtered_feature_bc_matrix.h5'))

    col_ix = adata.var_names.get_loc(gene_name)
    this_col = adata.X[:, col_ix].toarray().flatten()
    nof_rows = this_col.shape[0] # One row is one spot

    # Put all counts in an image-like structure (slow, not the right way later
    # on, but to understand the structure)
    max_count = 0
    img = np.zeros((3350, 3350), dtype=np.uint8)
    for row_ix in range(nof_rows):

        if row_ix % 100000 == 0:
            print(f'Processed spots: {row_ix//100000}00k, max = {max_count}')

        count = this_col[row_ix]
        max_count = max(count, max_count)

        # Find x,y coordinate
        spot_name = adata.obs_names[row_ix]
        x, y = spotname_to_xy(spot_name)

        # Set the pixel
        img[y, x] = count

    save_ndarray_as_image(img, 1, output_prefix + '.gain01.png')
    save_ndarray_as_image(img, 2, output_prefix + '.gain02.png')
    save_ndarray_as_image(img, 4, output_prefix + '.gain04.png')
    save_ndarray_as_image(img, 8, output_prefix + '.gain08.png')
    save_ndarray_as_image(img, 16, output_prefix + '.gain16.png')
    save_ndarray_as_image(img, 32, output_prefix + '.gain32.png')
    save_ndarray_as_image(img, 64, output_prefix + '.gain64.png')
    print('Done')


def lab_save_count_image(
    input_path: str, 
    output_prefix: str,
    load_h5ad = False):
    # TODO: Duplicated from lab_cll2a_st_image (make util function + deduplicate)

    if load_h5ad:
        adata = anndata.io.read_h5ad(input_path)
    else:
        adata = sc.read_10x_h5(os.path.join(input_path, 'filtered_feature_bc_matrix.h5'))

    counts = adata.X.sum(1) 

    # Put all counts in an image-like structure (slow, not the right way later
    # on, but to understand the structure)
    max_count = 0
    img = np.zeros((3350, 3350), dtype=np.uint8) # TODO: Don't hard-code
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
    save_ndarray_as_image(img, 1,  output_prefix + '.gain01.png')
    save_ndarray_as_image(img, 2,  output_prefix + '.gain02.png')
    save_ndarray_as_image(img, 4,  output_prefix + '.gain04.png')
    save_ndarray_as_image(img, 8,  output_prefix + '.gain08.png')
    save_ndarray_as_image(img, 16, output_prefix + '.gain16.png')
    save_ndarray_as_image(img, 32, output_prefix + '.gain32.png')
    save_ndarray_as_image(img, 64, output_prefix + '.gain64.png')
    save_ndarray_as_image(img, 96, output_prefix + '.gain96.png')
    print('Done!')


def lab_xy_normalize(input_dir: str, gene_name: str):
    # Lab stuff for normalizing the 2um data by row/col means, to compensate
    # for the line-like artifacts. Conclusion: This normalization seems
    # reasonable, the data looks much better this way.
    #
    # Some problems with this approach:
    # - The computed mean values are sensitive to what we find in the data
    #   Both variations between regions and small peaks (like the IGKCs) will
    #   be destroyed.
    #
    # - We destroy the absolute intensity, which makes it hard to compare the
    #   overall level between two different images.

    img = cv2.imread(os.path.join(input_dir, f'onegene_{gene_name}.gain01.png'), cv2.IMREAD_GRAYSCALE)

    mean_rows = np.mean(img, axis=0, keepdims=True) # Mean over rows
    mean_cols = np.mean(img, axis=1, keepdims=True) # Mean over cols

    img_f = img.astype(float) / mean_rows
    img_f = img_f / mean_cols

    figure, axis = plt.subplots(2, 1)

    # Save normalized image
    # save_ndarray_as_image(img_f, 0.25,  os.path.join(input_dir, f'onegene_{gene_name}_rowcolnorm.gain0.25.png'))
    # save_ndarray_as_image(img_f, 0.5,  os.path.join(input_dir, f'onegene_{gene_name}_rowcolnorm.gain0.50.png'))
    # save_ndarray_as_image(img_f, 1,  os.path.join(input_dir, f'onegene_{gene_name}_rowcolnorm.gain01.png'))
    # save_ndarray_as_image(img_f, 2,  os.path.join(input_dir, f'onegene_{gene_name}_rowcolnorm.gain02.png'))
    save_ndarray_as_image(img_f, 4,  os.path.join(input_dir, f'onegene_{gene_name}_rowcolnorm.gain04.png'))
    save_ndarray_as_image(img_f, 8,  os.path.join(input_dir, f'onegene_{gene_name}_rowcolnorm.gain08.png'))
    # save_ndarray_as_image(img_f, 16, os.path.join(input_dir, f'onegene_{gene_name}_rowcolnorm.gain16.png'))


def lab_plot_profiles():

    input_dir = '_temp/CLL1D'
    gene_name = 'IGKC'
    img = cv2.imread(os.path.join(input_dir, f'onegene_{gene_name}.gain01.png'), cv2.IMREAD_GRAYSCALE)


    if True:
        # Selecting 3 columns that we know are interesting from this data
        col1 = img[:, 2104]
        col2 = img[:, 2105]
        col3 = img[:, 2106]
        sigma = 31.0
        col1f = gaussian_filter1d(col1.astype(float), sigma)
        col2f = gaussian_filter1d(col2.astype(float), sigma)
        col3f = gaussian_filter1d(col3.astype(float), sigma)

        plt.plot(col1f, color='red')
        plt.plot(col2f * 0.17, color='green')
        plt.plot(col3f, color='blue')
    
    if False:
        # Selecting 3 rows that we know are interesting from this data
        # Looks reasonably linear as well
        row1 = img[2683, :]
        row2 = img[2684, :]
        row3 = img[2685, :]
        sigma = 15.0
        row1f = gaussian_filter1d(row1.astype(float), sigma)
        row2f = gaussian_filter1d(row2.astype(float), sigma)
        row3f = gaussian_filter1d(row3.astype(float), sigma)

        plt.plot(row1f, color='red')
        plt.plot(row2f, color='green')
        plt.plot(row3f, color='blue')

    plt.show()
    plt.waitforbuttonpress()


def lab_xy_normalize2(input_dir: str, gene_name: str):
    # Lab stuff for normalizing the 2um data by row/col means, to compensate
    # for the line-like artifacts.
    # 
    # Improved version, normalizing by supressing peaks, measured by comparing
    # median-filtered row/col profile with actual profiles.
    # 
    # This works better than version 1 above with respect to maintaining overall
    # absolute intensities and keeping slowly varying trends. However, we still
    # have a problem that small high-intensity peaks affect the mean values too
    # much.

    img = cv2.imread(os.path.join(input_dir, f'onegene_{gene_name}.gain01.png'), cv2.IMREAD_GRAYSCALE)

    mean_rows = np.mean(img, axis=0, keepdims=True) # Mean over rows
    mean_cols = np.mean(img, axis=1, keepdims=True) # Mean over cols

    # Compute the normalization
    mean_rows_f = medfilt(mean_rows.flatten(), 51)
    mean_cols_f = medfilt(mean_cols.flatten(), 51)

    row_factor = mean_rows_f / (mean_rows.flatten() + 1e-6)
    col_factor = mean_cols_f / (mean_cols.flatten() + 1e-6)

    # row_factor[row_factor >= 1.0] = 1.0
    # col_factor[col_factor >= 1.0] = 1.0

    row_factor = row_factor.reshape(mean_rows.shape)
    col_factor = col_factor.reshape(mean_cols.shape)

    # Apply the normalization
    img_f = img.astype(float) * row_factor * col_factor

    # Save normalized image
    # save_ndarray_as_image(img_f, 1,  os.path.join(input_dir, f'onegene_{gene_name}_rowcolnorm2.gain01.png'))
    # save_ndarray_as_image(img_f, 2,  os.path.join(input_dir, f'onegene_{gene_name}_rowcolnorm2.gain02.png'))
    # save_ndarray_as_image(img_f, 4,  os.path.join(input_dir, f'onegene_{gene_name}_rowcolnorm2.gain04.png'))
    # save_ndarray_as_image(img_f, 8,  os.path.join(input_dir, f'onegene_{gene_name}_rowcolnorm2.gain08.png'))
    save_ndarray_as_image(img_f, 16, os.path.join(input_dir, f'onegene_{gene_name}_rowcolnorm2.gain16.png'))
    save_ndarray_as_image(img_f, 32, os.path.join(input_dir, f'onegene_{gene_name}_rowcolnorm2.gain32.png'))
    save_ndarray_as_image(img_f, 64, os.path.join(input_dir, f'onegene_{gene_name}_rowcolnorm2.gain64.png'))
    save_ndarray_as_image(img_f, 96, os.path.join(input_dir, f'onegene_{gene_name}_rowcolnorm2.gain96.png'))


def get_rowcol_normalization_factors(img: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    # Return scaling factors (rows, columns) of the input image, using 
    # the best known row/col-wise normalization that we found. Currently
    # using a large median filter etc, but this could be adjusted.
    # 
    # This is still early lab code. Parameters should be added later.
    mean_rows = np.mean(img, axis=0, keepdims=True) # Mean over rows
    mean_cols = np.mean(img, axis=1, keepdims=True) # Mean over cols

    # Compute the normalization factors
    mean_rows_f = medfilt(mean_rows.flatten(), 51)
    mean_cols_f = medfilt(mean_cols.flatten(), 51)

    row_factors = mean_rows_f / (mean_rows.flatten() + 1e-6)
    col_factors = mean_cols_f / (mean_cols.flatten() + 1e-6)

    return (row_factors, col_factors)



def lab_precompute_xy_lookup():
    
    input_dir = '/media/erik/T9/run1_D/outs/binned_outputs/square_002um/'
    adata = sc.read_10x_h5(os.path.join(input_dir, 'filtered_feature_bc_matrix.h5'))

    obs_names = adata.obs_names.to_numpy()
    nof_spots = obs_names.shape[0]

    lookup = np.zeros((nof_spots, 2), dtype=np.int32)

    for spot_ix in range(nof_spots):

        if spot_ix % 100000 == 0:
            print(f'Processed spots: {spot_ix//100000}00k, max = {nof_spots}')

        x, y = spotname_to_xy(obs_names[spot_ix])
        lookup[spot_ix, 0] = x
        lookup[spot_ix, 1] = y
    
    np.save('_temp/CLL1D/spot_xy_lookup.npy', lookup)


def lab_xy_normalize_all(input_dir: str):
    # Trying to normalize the entire data (all genes) by normalizing each one
    # individually. First, just do this and look at the final complete count
    # image, to see if it looks more correctly normalized now.
    # 
    # TODO: Need to consider how to deal with genes that are very low-expressed
    # Perhaps just skip them for now (only apply the normalization for
    # genes that are above a certain count
    # 
    # TODO: If we downweight some but not all genes, we should perhaps downweight
    # the others as well, since they will otherwise be biased in comparison.
    # Alternatively, when we adjust each gene, we should do so in a way that
    # doesn't affect the total expression level for that gene, just the balance
    # between rows/columns.
    # 
    # TODO(PRIO): This function is broken! If we extract the data for just
    # IKGC later and just compute the single-gene counts for that gene, the
    # results are not right (lines still visible)

    # Parameters (TODO: Make input struct)
    debug = False
    count_threshold = 200000 # This will include ~25 genes in the normalization for 1D, just as an ad-hoc first attempt
    input_dir = '/media/erik/T9/run1_D/outs/binned_outputs/square_002um/'
    output_dir = '_temp/CLL1D'

    # Load stuff
    adata = sc.read_10x_h5(os.path.join(input_dir, 'filtered_feature_bc_matrix.h5'))
    spot_xy_lookup = np.load(os.path.join(output_dir, 'spot_xy_lookup.npy'))

    total_gene_counts = np.asarray(adata.X.sum(0)).flatten()
    nof_spots = adata.X.shape[0]
    nof_genes = adata.X.shape[1]

    # Temp, for debugging:
    if debug:
        nof_spots = 100

    for gene_ix in range(nof_genes):
        this_count = total_gene_counts[gene_ix]

        if this_count > count_threshold:
            print(f'gene_ix = {gene_ix}, this_count = {this_count}')

            # Apply col/row normalization for this gene
            X_this_gene = adata.X[:, gene_ix].toarray().flatten() # Flat list of counts for all spots

            # TODO: Break out
            print('Computing normalization factors')
            this_count_image = np.zeros((3350, 3350), dtype=np.float32) # TODO: Don't hard-code
            for spot_ix in range(nof_spots):
                x = spot_xy_lookup[spot_ix, 0]
                y = spot_xy_lookup[spot_ix, 1]

                count = X_this_gene[spot_ix]
                this_count_image[y, x] = count

            # Compute col/row scaling factors from this_count_image
            row_factors, col_factors = get_rowcol_normalization_factors(this_count_image)
            
            # Now apply them to the counts. Requires a second pass through the data
            # TODO: Research faster ways later. Especially when we want to do
            # something like this on several genes.
            print('Applying normalization factors')
            for spot_ix in range(nof_spots):
                count = X_this_gene[spot_ix]
                if count > 0:
                    x = spot_xy_lookup[spot_ix, 0]
                    y = spot_xy_lookup[spot_ix, 1]

                    count *= row_factors[y] * col_factors[x]
                    adata.X[spot_ix, gene_ix] = count
                    
    print('Saving normalized adata matrix')
    adata.write(os.path.join(output_dir, 'normalized_adata.h5ad'))
    print('Done!')


if __name__ == '__main__':

    # lab_save_gene_image('/media/erik/T9/run2_A/outs/binned_outputs/square_002um/', 'IGKC', '_temp/CLL2A/onegene_IGKC')
    # lab_xy_normalize('_temp/CLL2A', 'IGKC')
    # lab_xy_normalize2('_temp/CLL2A', 'IGKC')

    # lab_save_gene_image('/media/erik/T9/run1_D/outs/binned_outputs/square_002um/', 'IGKC', '_temp/CLL1D/onegene_IGKC')
    # lab_xy_normalize('_temp/CLL1D', 'IGKC')
    # lab_xy_normalize2('_temp/CLL1D', 'IGKC')

    # lab_save_gene_image('/media/erik/T9/run2_D/outs/binned_outputs/square_002um/', 'IGKC', '_temp/CLL2D/onegene_IGKC')
    # lab_xy_normalize('_temp/CLL2D', 'IGKC')
    # lab_xy_normalize2('_temp/CLL2D', 'IGKC')

    # lab_plot_profiles()

    # lab_precompute_xy_lookup()
    # lab_xy_normalize_all('/media/erik/T9/run1_D/outs/binned_outputs/square_002um/')

    # Very experimental: Save using the experimental gene-by-gen row/col normalization
    # lab_save_count_image('_temp/CLL1D/normalized_adata.h5ad', '_temp/CLL1D/counts_gxynorm', load_h5ad=True)

    lab_save_gene_image('_temp/CLL1D/normalized_adata.h5ad', 'IGKC', '_temp/CLL1D/hacky_test_igkc', load_h5ad=True)
