#!/usr/bin/env python
import numpy as np
import pandas as pd
import nibabel as nib

def corr2_coeff(A,B):
    '''Calculate correlation between seed voxels and target ROIS'''
    ssA = (A**2).sum(1)
    ssB = (B**2).sum(1)
    return np.dot(A, B.T) / np.sqrt(np.dot(ssA[:,None], ssB[None]))


def compute_correlation(vol, seed, ptseries, output)
    # Load cleaned voxel-level timeseries but only keep seed 
    c = nib.load(vol)

    structures = {}
    for brain_model in c.header.get_index_map(1).brain_models:
        name  = brain_model.brain_structure.split("STRUCTURE_",1)[1]
        if name == "BRAIN_STEM": # CIFTI doesn't like custom names so stored as BRAIN_STEM
            name = seed
            stype = brain_model.model_type.split("TYPE_",1)[1]
            start = brain_model.index_offset
            end   = brain_model.index_offset+brain_model.index_count
            structures[name] = [stype, (start, end), c.get_fdata()[:,start:end]]
    seeddata = structures[seed][2]

    # Load cleaned parcellated timeseries
    p = nib.load(ptseries)
    pdata = p.get_fdata()

    parcels = {}
    n_parcels = p.get_fdata().shape[1]
    for ip in np.arange(0,n_parcels):
        parcel = p.header.get_index_map(1)[ip+3]
        name = parcel.name
        if name == "BRAIN_STEM":
            name = seed
            indices = parcel.voxel_indices_ijk
        if parcel.vertices:
            stype = 'CORTEX'
        else:
            stype = 'SUBCORTEX'     
        parcels[name] = [ip, stype] #indices, p.get_fdata()[:,ip]

    # Extract labels from parcellated timeseries
    labels = [d for d in parcels]

    # Calculate correlation between seed and targets
    corrmatrix = corr2_coeff(seeddata.T, pdata.T)
    fzmatrix =  np.arctanh(corrmatrix)

    # Write to pandas dataframe
    df = pd.DataFrame(fzmatrix,columns=labels)

    # Write to output files
    np.savez(output, indices=indices,corr=fzmatrix)

if __name__ == '__main__':
    compute_correlation(
        vol=snakemake.input.vol,
        seed=snakemake.params.seed[0],
        ptseries=snakemake.input.ptseries,
        output=snakemake.output.correlation,
    )