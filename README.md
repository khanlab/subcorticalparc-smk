# Snakemake workflow: zonaconn-smk
Snakemake workflow for diffusion connectivity with the zona 
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

 ## Tractmap processing
 Click the toggle to see recommendations for processing data on Graham (via
 Digital Alliance).
 <details>
 <summary>Recommendations</summary>
 Processing is best done using local scratch <code>$SLURM_TMPDIR</code> since 
 it is very I/O intensive. Since there isn't an easy way to achieve this 
 built-in to snakemake, another job script has been created that copies the 
 data to local scratch, runs snakemake to generate a single subject's 
 tractmaps, then copies that data back before the job finishes.
 
 <h3>Pre-requisites</h3>
 This requires <a href="https://github.com/pvandyken/kslurm"><code>kslurm
 </code></a> to be installed. You can do this by pasting the following in a
 terminal:

 <pre><code>curl -sSL https://raw.githubusercontent.com/pvandyken/kslurm/master/install_kslurm.py | python -</code></pre>

 <h3>Steps</h3>
 1. All processing up to tract maps is complete. 

 <pre><code>#this command does a dry-run to see if any jobs need to be run still
     snakemake --omit-from resample_clus_seed -npr

     #this will run all the jobs to complete processing before tractmaps
     snakemake --omit-from resample_clus_seed --profile cc-slurm</code></pre>

 2. The <code>results/tractmap</code> folder for a subject should be empty before submitting any new tractmap jobs.

 To run tract maps for one subject (e.g. sub-100307) you can use:
 <pre><code>kbatch 3:00 8 32G gpu -a ctb-akhanf_cpu . ./job_tractmaps sub-100307
</code></pre>

 To run on all subjects in a participants.tsv file, use:
 <pre><code>./submit_tractmaps config/participants.tsv</code></pre>
 </details>
 <br>

## Authors

* Ali Khan @akhanf 
* Sudesna Chakraborty
* Jason Kai

## Reference

If you use this workflow in a paper, don't forget to give credits to the 
authors by citing the URL of this (original) repository and, if available, 
its DOI (see above).