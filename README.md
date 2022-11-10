# Snakemake workflow: zonaconn-smk
Snakemake workflow for diffusion and functional connectivity with the zona 
incerta. 

_This workflow is currently customized to run with data from the HCP1200 7T_

## Description
Using the HCP-MMP cortical parcellation (180 regions, left/right sym labels) 
as targets and performs probabilistic tracking from the zona incerta seed in 
each subject's native space. The connectivity data from seed voxels are brought
into the template space to perform spectral clustering on the concatenated 
feature vectors to parcellate into `k` regions.

<!-- To be updated -->
### Inputs
- Probabilistic segmentation(s) as 3D NIFTI for ZI on a single MNI template 
space
- participants.tsv with target subject IDs
- For each target subject:
  - Freesurfer processed data
  - Pre-processed DWI, registered to T1w space 
  (e.g. HCP-style, or from [prepdwi](https://github.com/khanlab/prepdwi))
  - BEDPOST processed data 
  - Transformations from *template* T1w space to/from each subject T1w, e.g. 
  from: 
  [ants_build_template_smk](https://github.com/akhanf/ants_build_template_smk);
  must include affine, warp and invwarp

### Singularity containers
 - Freesurfer (for `mri_convert`, `mris_convert`, `mri_info`)
 - Connectome workbench
 - Neuroglia (contains FSL, ANTS, gnu parallel etc..)
 - FSL6 with CUDA container 

_NOTE: Currently the tractography step in the workflow requires a GPU_
 

## RECOMMENDED EXECUTION: 
snakemake -np --profile cc-slurm

 - Further job grouping (beyond 1 job per subject) on graham with 
 `--group-components` is not recommended as it can lead to 
 over-subscribed GPUs when run on graham (TODO: fix this)

## Authors

* Ali Khan @akhanf 
* Sudesna Chakraborty
* Jason Kai

## Reference

If you use this workflow in a paper, don't forget to give credits to the 
authors by citing the URL of this (original) repository and, if available, 
its DOI (see above).