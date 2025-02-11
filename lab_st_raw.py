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

    data_dir = '/home/erik/encdata/pixwse/blobs/datasets/broad/10x_samples/human_lymph_node_modded'

    adata = sq.read.visium(
        data_dir,
        counts_file='filtered_feature_bc_matrix.h5')

    adata.var_names_make_unique() # TODO: Vad gör denna?

    # adata.var_names.str.startswith("MT-") returns an ndarray with bools
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    # Example stuff that can be done
    # adata.obs['total_counts'].to_numpy()
    # sc.pp.filter_cells(adata, min_counts=5000)
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)

    # PCA and clustering
    sc.pp.pca(adata)
    # Produces a 4039x50 matrix in adata.obsm['X_pca']
    # This looks like 50 PCA components.
    # TODO: Where are the eigenvalues?
    # TODO: Where is the 50-channel image?

    # This plots 'total_counts' and 'n_genes_by_counts' (TODO: What is that exactly?)
    # overlayed on the HE image. 
    
    # According to the doc of sc.pl.spatial, assuming that coordinates are 
    # in image pixel coordinates. All images in 'spatial' (except the lowres one)
    # is 1920x2000

    # Reverse-engineering to find what is used here:
    #   adata.uns['spatial']['V1_Human_Lymph_Node']['images']['hires'] - numpy array containing the bg image
    #   The circle radius is computed like: circle_radius = size * scale_factor * spot_size * 0.5
    sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts"])

    sc.pp.neighbors(adata)

    # path, *, counts_file='filtered_feature_bc_matrix.h5', library_id=None, load_images=True, source_image_path=None, **kwargs)
    print('Done')


# ----------------------------------------------------------------------------
# Understanding the input data

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

# ----------------------------------------------------------------------------
# Understanding the binned data

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


if __name__ == '__main__':
    # sc.logging.print_versions()
    # print("hej")

    lab_scanpy_ts1()
    # lab_h5_clldata_binned_002um()
