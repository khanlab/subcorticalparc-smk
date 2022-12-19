#!/usr/bin/env python 
import numpy as np
import pandas as pd
import nibabel as nib
from nibabel.nifti1 import Nifti1Image
from linear.input_data import NiftiLabelsMasker

def extract_confounds(rois, vol, movreg, output):
    atlas = nib.load(rois)
    labels = atlas.get_fdata()

    # Relabel and extract timeseries for CSF and WM labels
    img = np.zeros(labels.shape)
    img[(labels == 4) | (labels == 43)] = 1 # CSF
    img[labels > 5000] = 2 # WM

    tmp = Nifti1Image(img, atlas.affine, atlas.header)
    masker = NiftiLabelsMasker(labels_img=tmp, standardize=False)
    time_series1 = masker.fit_transform(vol)

    # Then for whole brain (i.e., 'global')
    brain = np.zeros(labels.shape)
    brain[labels != 0] = 1

    tmp = Nifti1Image(brain, atlas.affine, atlas.header)
    masker = NiftiLabelsMasker(labels_img=tmp, standardize=False)
    time_series2 = masker.fit_transform(vol)

    # Concatenate timeseries
    df1 = pd.DataFrame(
        {
            'CSF': time_series1[:, 0],
            'WhiteMatter': time_series1[:, 1], 
            'GlobalSignal': time_series2[:, 0],
        }
    )

    # Load movement parameters (and their derivatives)
    names = [
        'X','Y','Z','RotX','RotY','RotZ','Xd','Yd','Zd','RotXd','RotYd','RotZd'
    ]
    df2 = pd.read_csv(movreg, names=names, header=None, delim_whitespace=True)

    # Put everything together (excluding derivatives) 
    # and write to .tsv file for 'ciftify_clean_img' 
    df_concat = pd.concat([df2.iloc[:,:6], df1], axis=1)
    df_concat.to_csv(output, sep='\t')

if __name__ == '__main__':
    extract_confounds(
        rois=snakemake.input.rois,
        vol=snakemake.input.vol,
        movreg=snakemake.input.movreg,
        output=snakemake.output.confounds,
    )