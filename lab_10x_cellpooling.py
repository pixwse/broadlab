
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


def lab_read_visium2():
    # Read the Visium HD using scanpy and play around a bit
    # TODO: Move to some scanpy_lab file!

    plt.interactive(True)

    input_path = '/media/erik/T9/run2_A/outs/binned_outputs/square_002um/'

    # Load Visium HD data and spatial coordinates of spots
    df_pos = pd.read_parquet(os.path.join(input_path, 'spatial/tissue_positions.parquet'))
    adata = sc.read_10x_h5(os.path.join(input_path, 'filtered_feature_bc_matrix.h5'))

    adata.X
    # adata.X is now a CSR matrix with shape 10732334x18085
    # Here, 18085 is the number of genes and 10732334 is some form of
    # encoded 

    counts = adata.X.sum(1)  # Sum over columns

    # Reshape to an image (should not be done like this, but just to
    # understand the structure)

    # This is not the right way!
    # There is obs_names and var_names that are connected to the
    # rows/columns (columns/rows?) in X. We need to extract the right
    # spot name and map to corresponding x/y coordinate. This
    # coordinate is also encoded in the spot name (but we should perhaps not
    # rely on that).

    img = np.zeros((3350, 3350))
    for ix, c in enumerate(counts):
        x = ix // 3350
        y = ix % 3350 
        # TODO: Or the other way around
        
        nc = c.item() / 10.0
        nc = min(nc, 1.0)
        img[y, x] = nc
    
    img_u8 = (img * 255).astype(np.uint8)
    img_pil = PIL.Image.fromarray(img_u8)
    img_pil.save('_temp/temp.png')

    plt.imshow(img)
    plt.waitforbuttonpress()

    print('Done!')


if __name__ == '__main__':
    lab_read_visium2()
    # lab_10x_pooling()
