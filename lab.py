import squidpy as sq
import h5py

import scanpy as sc
import anndata as ad

# Data retrieval
import pooch

import os

import matplotlib.pyplot as plt
import seaborn as sns

import PIL
import cv2
import large_image
import numpy as np
import io


def lab_scanpy_ts1():
    # Reproducing Spatial Transcriptomics example from
    # https://scanpy.readthedocs.io/en/stable/tutorials/spatial/basic-analysis.html

    # Prepared the data by downloading everything from 
    # https://support.10xgenomics.com/spatial-gene-expression/datasets/1.0.0/V1_Human_Lymph_Node
    # unzipped V1_Human_Lymph_Node_spatial.tar.gz and removed the intermediate
    # dir, such that everything is directly under 'spatial'

    data_dir = '/home/erik/encdata/pixwse/blobs/datasets/broad/10x_samples/human_lymph_node'

    adata = sq.read.visium(
        data_dir,
        counts_file='V1_Human_Lymph_Node_filtered_feature_bc_matrix.h5')
    

    adata.var_names_make_unique() # TODO: Vad gör denna?

    # adata.var_names.str.startswith("MT-") returns an ndarray with bools
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    # path, *, counts_file='filtered_feature_bc_matrix.h5', library_id=None, load_images=True, source_image_path=None, **kwargs)
    print('Done')


def lab_h5_clldata_original():
    # Labbing with the raw space ranger output data in HDF5 format from 
    # Johan's CLL dataset 

    data_dir = '/media/erik/T9/run2_A/outs'

    with h5py.File(os.path.join(data_dir, 'feature_slice.h5'), "r") as f:
        # Print all root level object names (aka keys) 
        # these can be group or dataset names 
        print("Keys: %s" % f.keys())
        
        # Can try things like
        # min(f['feature_slices']['1000']['row'][:])
        # qq = [int(k) for k in f['feature_slices'].keys()]

        print('Done')


def lab_h5_clldata_binned_002um():
    # Labbing with the raw space ranger output data in HDF5 format from 
    # Johan's CLL dataset 

    data_dir = '/media/erik/T9/run2_A/outs'

    with h5py.File(os.path.join(data_dir, 'binned_outputs/square_002um/raw_feature_bc_matrix.h5'), "r") as f:
        # Print all root level object names (aka keys) 
        # these can be group or dataset names 
        print("Keys: %s" % f.keys())
        
        # Can try things like
        # min(f['feature_slices']['1000']['row'][:])
        # qq = [int(k) for k in f['feature_slices'].keys()]

        print('Done')


def lab_load_large_image():
    path = '/media/erik/T9/HD_images/run2_A/2024-08-29SRParryRichterPD1460124.tif'

    # Doesn't work using directly, says:
    # Image size (9954918400 pixels) exceeds limit of 178956970 pixels
    # (meaning that it's around 100k x 100k large)
    # img = PIL.Image.open(path)

    # Testing with OpenCV
    # Only supports up to 2^30, and can't be configure to more than that
    # img = cv2.imread(path)

    # Testing large-image (with pip install large-image[tiff])
    ts = large_image.getTileSource(path)
    qq = ts.getMetadata()

    plt.interactive(True)

    for tile_info in ts.tileIterator(
            region=dict(left=40000, top=40000, width=20000, height=20000, units='base_pixels'),
            scale=dict(magnification=20),
            tile_size=dict(width=1000, height=1000),
            tile_overlap=dict(x=50, y=50),
            format=large_image.tilesource.TILE_FORMAT_PIL
            ):

        im_tile = np.array(tile_info['tile'])        
        plt.imshow(im_tile)
        plt.show()
        print('One tile ok')
        plt.waitforbuttonpress()


def split_large_image():
    # For example: tile 166 is within the ST region (top part of the
    # large crack-like structure)

    path = '/media/erik/T9/HD_images/run2_A/2024-08-29SRParryRichterPD1460124.tif'
    output_path = '_temp/cll2a_tiles_png'

    # Testing large-image (with pip install large-image[tiff])
    ts = large_image.getTileSource(path)
    qq = ts.getMetadata()

    if False:
        region, mime_type = ts.getRegion(
            output={'maxWidth': 10000, 'maxHeight': 10000, 'format': 'numpy'}
        )

        image = PIL.Image.open(io.BytesIO(region))
        im.save(os.path.join(output_path, 'lowres.jpg'))

    tile_no = 0
    for tile_info in ts.tileIterator(
            # region=dict(left=0, top=0, width=20000, height=20000, units='base_pixels'),
            # scale=dict(magnification=20),
            tile_size=dict(width=5000, height=5000),
            tile_overlap=dict(x=0, y=0),
            format=large_image.tilesource.TILE_FORMAT_PIL
            ):

        if tile_no > 150 and tile_no < 170:
            im_tile = np.array(tile_info['tile'])        
            im = PIL.Image.fromarray(im_tile)
            im.save(os.path.join(output_path, f"tile{tile_no}.png"))
            print(f'Tile {tile_no} saved')

        tile_no += 1


if __name__ == '__main__':
    # sc.logging.print_versions()
    # print("hej")

    # lab_scanpy_ts1()
    # lab_h5_clldata_binned_002um()

    # lab_load_large_image()
    split_large_image()
