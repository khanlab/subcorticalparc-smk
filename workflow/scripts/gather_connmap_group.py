import numpy as np


def gather_connmap_group(subject_npz, group_npz):
    # Grab number of subjects
    nsubjects = len(subject_npz)

    # Load first file to extract shape
    data = np.load(subject_npz[0])
    conn_shape = data["conn"].shape

    # Get affine and mask
    affine = data["affine"]
    mask = data["mask"]

    # Aggregate date
    conn_group = np.zeros([nsubjects, conn_shape[0], conn_shape[1]])
    for i, npz in enumerate(subject_npz):
        data = np.load(npz)
        conn_group[i, :, :] = data["conn"]

    # Save conn_group
    np.savez(
        group_npz,
        conn_group=conn_group,
        mask=mask,
        affine=affine,
    )


if __name__ == "__main__":
    gather_connmap_group(
        subject_npz=snakemake.input.connmap_npz,
        group_npz=snakemake.output.connmap_group_npz,
    )
