#!/usr/bin/env python
import nibabel as nib
import numpy as np
from sklearn.cluster import SpectralClustering


# Define a function for saving niftis
def save_label_nii(labels, mask, affine, out_nifti):
    labels_vol = np.zeros(mask.shape)
    labels_vol[mask > 0] = labels + 1  # add a 1 so label 0 is diff from bgnd
    labels_nib = nib.Nifti1Image(labels_vol, affine)
    nib.save(labels_nib, out_nifti)


def spectral_cluster(connmap_group, max_k, out_nii_list):
    # Set up variables
    cluster_range = range(2, max_k + 1)
    out_nii_list = out_nii_list

    # Load data
    data = np.load(connmap_group)
    conn_group = data["conn_group"]
    mask = data["mask"]
    affine = data["affine"]

    # Concatenate subjects
    conn_group_m = np.moveaxis(conn_group, 0, 2)
    conn_concat = conn_group_m.reshape(
        [conn_group_m.shape[0], conn_group_m.shape[1] * conn_group_m.shape[2]]
    )

    # Perform spectral clustering and save output
    for i, k in enumerate(cluster_range):
        clustering = SpectralClustering(
            n_clusters=k,
            assign_labels="discretize",
            random_state=0,
            affinity="cosine",
        ).fit(conn_concat)

        print(f"i={i}, k={k}, saving {out_nii_list[i]}")
        save_label_nii(clustering.labels_, mask, affine, out_nii_list[i])


if __name__ == "__main__":
    spectral_cluster(
        connmap_group=snakemake.input.connmap_group_npz,
        max_k=snakemake.params.max_k,
        out_nii_list=snakemake.output,
    )
