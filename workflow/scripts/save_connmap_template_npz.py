#!/usr/bin/env python
import nibabel as nib
import numpy as np


def save_connmap_template(mask, connmap_3d, connmap_npz):
    # Load mask and get volume
    mask_nii = nib.load(mask)
    mask_vol = mask_nii.get_fdata()
    mask_indices = mask_vol > 0
    masked = mask_vol[mask_indices]

    # Get number of voxels in mask and connectivity map
    nvoxels = len(masked)
    ntargets = len(connmap_3d)

    # Gather and save template connectivity
    conn = np.zeros((nvoxels, ntargets))
    for i, conn_file in enumerate(connmap_3d):
        vol = nib.load(conn_file).get_fdata()
        masked = vol[mask_indices].T
        conn[:, i] = masked

    np.savez(
        connmap_npz,
        conn=conn,
        mask=mask_vol,
        affine=mask_nii.affine,
    )


if __name__ == "__main__":
    save_connmap_template(
        mask=snakemake.input.mask,
        connmap_3d=snakemake.params.connmap_3d,
        connmap_npz=snakemake.output.connmap_npz,
    )
