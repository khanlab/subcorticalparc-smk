<<<<<<< HEAD:workflow/scripts/funcparc/combine_correlation.py
#!/bin/env python
import numpy as np

def combine_correlation(correlation, output):
    '''Concatenate correlation matrices across subjects'''
    data = np.load(correlation[0])
    n_subjects = len(correlation)
=======
#!/usr/bin/env python
import numpy


def combine_correlation(correlation, output):
    """Concatenate correlation matrices across subjects"""
    data = np.load(correlation)
    n_subjects = len(data)
>>>>>>> 28d0383... restructure for snakebids:subcorticalparc_smk/workflow/scripts/funcparc/combine_correlation.py
    combined = np.zeros(
        [n_subjects, data["corr"].shape[0], data["corr"].shape[1]]
    )

    for i, npz in enumerate(correlation):
        data = np.load(npz)
        combined[i, :, :] = data["corr"]

    np.savez(output, corr_group=combined, indices=data["indices"])


if __name__ == "__main__":
    combine_correlation(
        correlation=snakemake.input.correlation,
        output=snakemake.output.combined_corr,
    )
