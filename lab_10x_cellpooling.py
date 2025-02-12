
import os

import geopandas as gpd
from shapely.geometry import Polygon, Point
from scipy import sparse
import scanpy as sc
import pandas as pd

import numpy as np
import matplotlib.pyplot as plt
import PIL

# Reproducing stuff from the 10x cell-level pooling experiment


def lab_10x_pooling():

    st_input_dir = '/media/erik/T9/run2_A/outs/binned_outputs/square_002um'
    pass


def lab_read_visium():
    # Read the Visium HD the same way that they do in the example
    input_path = '/media/erik/T9/run2_A/outs/binned_outputs/square_002um/'

    # Load Visium HD data and spatial coordinates of spots
    adata = sc.read_10x_h5(os.path.join(input_path, 'filtered_feature_bc_matrix.h5'))
    df_pos = pd.read_parquet(os.path.join(input_path, 'spatial/tissue_positions.parquet'))

    # Set the index of the dataframe to the barcodeslab_read_visium
    df_pos = df_pos.set_index('barcode')

    # Create an index in the dataframe to check joins (TODO: What happens here?)
    df_pos['index'] = df_pos.index

    # Adding the tissue positions to the meta data (TODO: What happens here?)
    adata.obs = pd.merge(adata.obs, df_pos, left_index=True, right_index=True)

    # Create a GeoDataFrame from the DataFrame of coordinates
    geometry = [Point(xy) for xy in zip(df_pos['pxl_col_in_fullres'], df_pos['pxl_row_in_fullres'])]
    gdf_coordinates = gpd.GeoDataFrame(df_pos, geometry=geometry)



if __name__ == '__main__':
    lab_read_visium2()
    # lab_10x_pooling()
