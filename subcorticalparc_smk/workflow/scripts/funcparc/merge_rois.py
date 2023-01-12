#!/bin/env python
import nibabel as nib
import numpy as np

def merge_roi(atlas, roi, out_fpath):
    # Load atlas and roi
    atlas = nib.load(atlas)
    roi = nib.load(roi)
    roi = roi.get_fdata()

    # Update atlas
    out_atlas = atlas.get_fdata()
    out_atlas[out_atlas == 16] = 0
    out_atlas[roi == 1] = 16


    # Save merged roi
    img = nib.Nifti1Image(out_atlas, affine=atlas.affine, header=atlas.header)
    nib.save(img, out_fpath)


if __name__ == "__main__":
    merge_roi(
        atlas=snakemake.input.atlas,
        roi=snakemake.input.roi,
        out_fpath=snakemake.output.out_fpath,
    )