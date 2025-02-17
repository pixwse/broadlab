# import squidpy as sq
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


# Idea:
# ---
# Even if some (most?) bins will consist of mixed cell types, perhaps the
# PCA still offers a reasonable basis to project the data into? However, 
# when we work with the 2um data, we just get individual counts, and
# it doesn't make sense to log transform. So, we would like directions in
# the data that are relevant for original data. Perhaps running PCA on the
# non-log-transformed data could give such directions? The first few PCAs
# will mostly contain individual highly expressed genes, but then we will
# get other stuff, where the PCA components might provide good sets of
# genes to aggregate over.
# 
# If so, a reasonable pipeline could be something like:
# - Start with the 2um data, normalize rows/cols instead of normalizing bins.
# - Bin down to 16um, run PCA on non-log-transformed data, save the top 50 PCs or so.
# - Revisit the 2um data, project onto PCA components, get 50 coefficients per
#   bin (2.1 Gb for the entire data if using float32)
# - Aggregate over cells. Can be done after projection as long as we're just
#   doing linear projectsions).
# - Cluster in this PCA space instead. Could include spatial location.
#
# Initial experiment:
# Goal: Find out if it is better to normalize rows/colums than to normalize
# each bin separately.
# Approach: Compare the two normalization procedures. Run some default
# clustering thingy afterwards. How to look at the data to see if the results
# are better? Just look in umap and see if things look separated?
# How to know where true cell types actually are?


def lab_scanpy_ts1():
    # Reproducing Spatial Transcriptomics example from
    # https://scanpy.readthedocs.io/en/stable/tutorials/spatial/basic-analysis.html
    # Doesn't exist anymore?

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


def lab_scanpy_pca():
    # Lab code to look at PCA vectors etc
    # 
    # Conclusion: If we don't run log transformation first, there are only
    # 3 PCs that dominate, clearly indicating a few very highly expressed genes.
    # When we do log transform, it is more spread out. This makes sense, and it
    # is reasonable to assume that the log transform is good.

    input_path = '/media/erik/T9/run1_D/outs/binned_outputs/square_016um/'

    # Load Visium HD data
    adata = sc.read_10x_h5(os.path.join(input_path, 'filtered_feature_bc_matrix.h5'))

    # Example stuff that can be done
    # adata.obs['total_counts'].to_numpy()
    # sc.pp.filter_cells(adata, min_counts=5000)
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)

    # PCA
    print('Running PCA')
    sc.pp.pca(adata) # Takes a few seconds
    # Results are stored in:
    #   - Input representation in the PCA basis:
    #     adata.obsm['X_pca'] (168778x50 matrix)
    #
    #   - Principal components: 
    #     adata.varm['PCs'] (18085x50 matrix)
    #
    #   - Variance per PC:
    #     adata.uns['pca']['variance']

    # Compute neighbor structure (needs to be done before UMAP)
    print('Running neighbors') # Take a few seconds (< 1min). Uses PCA if it exists. If it doesn't, it warns and runs PCA first (unless we explicitly say otherwise)
    sc.pp.neighbors(adata) # First need to create a neighborhood graph
    # TODO: Where is the output stored?

    # Compute UMAP structure (needs to be done before plotting)
    print('Running UMAP') 
    sc.tl.umap(adata)
    # TODO: Where is the output stored?

    # TODO: Save pre-computed stuff such that we don't need to rerun all the above

    # Plot (using the above precomputed stuff, stored in adata)
    print('Plotting UMAP')
    sc.pl.umap(adata, size=2, color='IGKC')

    # Example:
    # qq = adata.X.sum(1) # Total counts
    # adata.obs['erik_test'] = qq 
    # sc.pl.umap(adata, size=2, color='erik_test')

    plt.waitforbuttonpress()

    print('Done')


if __name__ == '__main__':
    # sc.logging.print_versions()
    # print("hej")

    # lab_scanpy_ts1()
    # lab_h5_clldata_binned_002um()
    lab_scanpy_pca()
