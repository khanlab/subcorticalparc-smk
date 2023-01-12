#!/usr/bin/env python
import numpy as np
import pandas as pd
import nibabel as nib
from sklearn.cluster import SpectralClustering

# Mute warning in case seed has same number of voxels as target ROIs
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

# Function to save niftis
def save_label_nii (label_img, affine, header, out_nifti):
    img = nib.Nifti1Image(label_img, affine=affine, header=header)
    nib.save(img, out_nifti)

def spectral_cluster(correlation, rois, max_k, out_nii, out_labels):
    # Load data
    data = np.load(correlation)
    correlation = data['corr_group']
    indices = data['indices']

    afile = rois
    atlas = nib.load(afile)
    atlas_data = atlas.get_fdata()

    # Reshape and concatenate subjects
    # Spectral clustering doesn't like negative input apparantly, or square
    corr = np.moveaxis(correlation, 0, 2)
    corr_concat = corr.reshape([corr.shape[0],corr.shape[1]*corr.shape[2]])
    corr_concat += 1 
    corr_concat[np.isnan(corr_concat)] = 1

    # Average correlation matrix
    # corr_avg = np.nanmean(correlation,axis=0)
    # corr_avg += 1
    # corr_avg[np.isnan(corr_avg)] = 1

    # Output
    out_nii_list = out_nii
    cluster_range = range(2, max_k+1)
    labels = np.zeros((corr_concat.shape[0], len(cluster_range)))

    # Run spectral clustering and save results to nifti
    for i, k in enumerate(cluster_range):
        clustering = SpectralClustering(
            n_clusters=k, 
            assign_labels="discretize",
            random_state=0,
            affinity='cosine'
        ).fit(corr_concat)
        labels[:,i] = clustering.labels_
        
        label_img = np.zeros(atlas_data.shape)
        for j in range(0, len(atlas_data[atlas_data==16])):
            label_img[indices[j][0], indices[j][1], indices[j][2]] = labels[j,i]+1
        print(f'i={i}, k={k},saving {out_nii_list[i]}')
        save_label_nii(label_img, atlas.affine, atlas.header, out_nii_list[i])

    # Save results to CSV file
    df = pd.DataFrame(labels, columns=cluster_range)
    df.to_csv(out_labels)

if __name__ == "__main__":
    spectral_cluster(
        correlation=snakemake.input.correlation,
        rois=snakemake.input.rois,
        max_k=snakemake.params.max_k,
        out_nii=snakemake.output.niftis,
        out_labels=snakemake.output.labels,
    )