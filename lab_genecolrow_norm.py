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


# Lab code for normalizing each row/column separately individually for each
# gene. Looking at some highly expressed genes first, then for others


def lab_save_gene_image(input_dir: str, gene_name: str, output_prefix: str):
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


# def save_total_counts_image(input_dir: str, output_dir: str):
    # Read the Visium HD using scanpy, sum all counts and save it directly
    # as a png image, where the pixel value is the total count (max 132, so
    # works fine)

    # # Load Visium HD data and sum over columns (genes)
    # adata = sc.read_10x_h5(os.path.join(input_dir, 'filtered_feature_bc_matrix.h5'))
    # counts = adata.X.sum(1) 

    # # Put all counts in an image-like structure (slow, not the right way later
    # # on, but to understand the structure)
    # max_count = 0
    # img = np.zeros((3350, 3350), dtype=np.uint8)
    # for ix, count in enumerate(counts):

    #     if ix % 100000 == 0:
    #         print(f'Processed spots: {ix//100000}00k, max = {max_count}')

    #     max_count = max(count, max_count)

    #     # Find x,y coordinate
    #     spot_name = adata.obs_names[ix]
    #     x, y = spotname_to_xy(spot_name)

    #     # Set the pixel
    #     img[y, x] = count

    # print(f'max_count = {max_count}')

    # # Convert to u8 and save with different gains, for easier visualization (a bit hacky for now)
    # img_pil = PIL.Image.fromarray(img)
    # img_pil.save(os.path.join(output_dir, 'counts.png'))

    # img *= 2
    # img_pil = PIL.Image.fromarray(img)
    # img_pil.save(os.path.join(output_dir, 'counts_gain2.png'))

    # img *= 2
    # img_pil = PIL.Image.fromarray(img)
    # img_pil.save(os.path.join(output_dir, 'counts_gain4.png'))

    # img *= 2
    # img_pil = PIL.Image.fromarray(img)
    # img_pil.save(os.path.join(output_dir, 'counts_gain8.png'))

    # print('Done!')


if __name__ == '__main__':

    # lab_save_gene_image('/media/erik/T9/run2_A/outs/binned_outputs/square_002um/', 'IGKC', '_temp/CLL2A/onegene_IGKC')
    lab_save_gene_image('/media/erik/T9/run1_D/outs/binned_outputs/square_002um/', 'IGKC', '_temp/CLL1D/onegene_IGKC')
