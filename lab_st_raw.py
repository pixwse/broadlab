import h5py
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import PIL

from utils import *

# Lab stuff for getting to know the raw HDF5 data format used in the
# spaceranger output. Discontinued, since it seems better to use the
# scanpy library

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


def lab_h5_molecule_info():
    # Labbing with the molecule_info.h5 file in the raw data

    data_dir = '/media/erik/T9/run1_A/outs'

    with h5py.File(os.path.join(data_dir, 'molecule_info.h5'), "r") as f:
        # Print all root level object names (aka keys) 
        # these can be group or dataset names 
        print("Keys: %s" % f.keys())
        
        # Actual UMI used for this specific molecule
        # qq = f['umis'][:]
        # max(qq)

        # How many times a specific molecule was found (this counts as one
        # count in the final data, since there is only one UMI).
        # f['counts'][5]
        # f['umi'][5]

        print('Done')


def lab_tissue_positions():

    input_path = '/media/erik/T9/run1_A/outs/binned_outputs/square_002um/'
    df_pos = pd.read_parquet(os.path.join(input_path, 'spatial/tissue_positions.parquet'))

    print('Done')


def lab_systematic_row_umi_stuff():
    # Simple experiment to check it there is something funky going on with
    # the UMIs, specifically, if some UMIs are encountered less frequently
    # in rows that look attenuatedin the output.
    # 
    # Preliminary conclusion: Yes, it looks as if some UMIs are always found
    # with a lower frequency than others, across all rows. Could be that some
    # UMIs are just used less frequently, or that the distribution of UMIs
    # across genes is not completely random. But since it looks like the rows
    # and colums have different gains, it might as well be that different UMIs
    # have different gains as well. However, there doesn't seem to be a
    # systematic relationshow with the rows, e.g. such that some UMIs stop
    # working for some rows. It rather looks like an overall UMI-specific
    # attenuation.

    data_dir = '/media/erik/T9/run1_A/outs'

    with h5py.File(os.path.join(data_dir, 'molecule_info.h5'), "r") as f:

        nof_rows = 3350
        max_umi = max(f['umi'][:])

        umi_row_counts = np.zeros((nof_rows, max_umi+1), dtype=np.int16)

        barcode_idx_array = f['barcode_idx'][:]
        umi_array = f['umi'][:]
        barcodes = np.array(f['barcodes'])
        nof_barcodes = barcodes.shape[0]

        # Precompute the y corresponding to each barcode ix
        bcix_y_lookup_file = '_temp/CLL_1A_bcix_y_lookup.npy'
        if not os.path.exists(bcix_y_lookup_file):
            print('Precomputing y')
            bcix_to_y_lookup = np.zeros((nof_barcodes), dtype=np.int32)
            for bc_ix in range(nof_barcodes):
                if bc_ix % 100000 == 0:
                    print(f'bc_ix={bc_ix//1000}k')

                bc = barcodes[bc_ix]
                x, y = byte_barcode_to_xy(bc)
                bcix_to_y_lookup[bc_ix] = y
            np.save(bcix_y_lookup_file, bcix_to_y_lookup)
        else:
            bcix_to_y_lookup = np.load(bcix_y_lookup_file)

        # Go through all molecules, save UMIs found in each row
        umi_row_counts_file = '_temp/CLL_1A_row_umi_counts.npy'
        if not os.path.exists(umi_row_counts_file):
            print('Processing molecules')
            nof_molecules = barcode_idx_array.shape[0]
            for mol_ix in range(nof_molecules):

                if mol_ix % 100000 == 0:
                    print(f'mol_ix={mol_ix//1000}k')

                bc_ix = barcode_idx_array[mol_ix]
                y = bcix_to_y_lookup[bc_ix]

                umi = umi_array[mol_ix]
                umi_row_counts[y][umi] += 1

            np.save(umi_row_counts_file, umi_row_counts)
        else:
            umi_row_counts = np.load(umi_row_counts_file)

        # Sum UMIs for each row in the image, corresponds to the bright/dark rows
        # that we have seen
        umisums = np.sum(umi_row_counts, axis=1)

        # Visualize row counts
        img = umi_row_counts[:, 0:50000]
        img = img.astype(np.uint8)
        img[img >= 1] = 255
        pimg = PIL.Image.fromarray(img)
        pimg.save('_temp/temp_unsorted.png')

        # Sort rows first (to rows in the output image will be rows that look
        # attenuated)
        ixs = np.argsort(umisums)
        umi_row_counts_sorted = umi_row_counts[ixs, :]
        img = umi_row_counts_sorted[:, 0:50000]
        img = img.astype(np.uint8)
        img[img >= 1] = 255
        pimg = PIL.Image.fromarray(img)
        pimg.save('_temp/temp_sorted.png')

        print('Done')


def lab_systematic_gene_umi_stuff():
    # Simple experiment to check it there is some systematic attenuation
    # going on such that specific UMIs are always attenuated compared to
    # others, regardless of which genes they are used for.

    data_dir = '/media/erik/T9/run1_A/outs'

    with h5py.File(os.path.join(data_dir, 'molecule_info.h5'), "r") as f:

        umi_array = f['umi'][:]
        gene_array = f['feature_idx'][:]

        # Looking at the first 2000 genes and 10000 UMIs, since the entire data is too large
        nof_genes = 2000
        nof_umis = 10000
        gene_umi_counts = np.zeros((nof_genes, nof_umis), dtype=np.int16)

        nof_molecules = gene_array.shape[0]
        for mol_ix in range(nof_molecules):

            if mol_ix % 100000 == 0:
                print(f'mol_ix={mol_ix//1000}k')

            gene_ix = gene_array[mol_ix]
            umi = umi_array[mol_ix]

            if gene_ix < nof_genes and umi < nof_umis:
                gene_umi_counts[gene_ix, umi] += 1

        np.save('_temp/CLL_1A_gene_umi_counts.npy', gene_umi_counts)

        img = gene_umi_counts.astype(np.uint8)
        img[img >= 1] = 255
        pimg = PIL.Image.fromarray(img)
        pimg.save('_temp/temp_gene_umi_bias.png')            

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

    # lab_h5_clldata_binned_002um()
    # lab_h5_molecule_info()
    # lab_tissue_positions()

    # lab_systematic_row_umi_stuff()
    lab_systematic_gene_umi_stuff()
