#!/usr/bin/env python
import nibabel as nib
import numpy as np

def merge_roi(atlas, roi):
    # Load atlas and update atlas
    atlas = nib.load(atlas)
    out_atlas = atlas.get_fdata()
    out_atlas[out_atlas == 16] = 0
    out_atlas[roi == 1] = 16

    # Load roi
    roi = nib.load(roi)
    roi = roi.get_fdata()

    # Save merged roi
    img = nib.Nifti1Image(atlas_new, affine=atlas.affine, header=atlas.header)
    nib.save(img, out_fpath))


if __name__ == "__main__":
    merge_roi(
        atlas=snakemake.input.atlas,
        roi=snakemake.input.roi,
        out_fpath=snakemake.output.out_fpath,
    )