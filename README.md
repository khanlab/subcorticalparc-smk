# Snakemake workflow: diffparc-smk
General-purpose Snakemake workflow for diffusion-based subcortical parcellation

This workflow is customized for the basal forebrain, run on HCP1200 7T diffusion data.

Uses HCP-MMP cortical parcellation (180 regions, left/right sym labels) as targets and performs probabilistic tracking from the subcortical seed in each subject's space, then brings connectivity data from seed voxels into template space and performs spectral clustering on the concatenated feature vectors to parcellate into k regions.

Inputs:
- Probabilistic segmentation(s) as 3D NIFTI for subcortical structure(s) on a single *template* T1w space
- participants.tsv with target subject IDs
- For each target subject:
  - Freesurfer processed data
  - Pre-processed DWI, registered to T1w space (e.g. HCP-style, or from [prepdwi](https://github.com/khanlab/prepdwi))
  - BEDPOST processed data 
  - ANTS transformations from *template* T1w space to/from each subject T1w, e.g. from: [ants_build_template_smk](https://github.com/akhanf/ants_build_template_smk); must include affine, warp and invwarp


Singularity containers required:
 - Freesurfer (for `mri_convert`, `mris_convert`, `mri_info`)
 - Connectome workbench
 - Neuroglia (contains FSL, ANTS, gnu parallel etc..)
 - FSL6 with CUDA container 

NOTE: Currently the tractography step in the workflow requires a GPU 
 

# RECOMMENDED EXECUTION: 
snakemake -np --profile cc-slurm

 - Further job grouping (beyond 1 job per subject) on graham with `--group-components` is not recommended as it can lead to over-subscribed GPUs when run on graham (TODO: fix this)


## Authors

* Ali Khan @akhanf 
* Sudesna Chakraborty

## Usage

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and, if available, its DOI (see above).

### Step 1: Obtain a copy of this workflow

1. Create a new github repository using this workflow [as a template](https://help.github.com/en/articles/creating-a-repository-from-a-template).
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the newly created repository to your local system, into the place where you want to perform the data analysis.

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files in the `config/` folder. Adjust `config.yml` to configure the workflow execution, and `participants.tsv` to specify your subjects.

### Step 3: Install Snakemake

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 4: Execute workflow

Test your configuration by performing a dry-run via

    snakemake --use-singularity -n

Execute the workflow locally via

    snakemake --use-singularity --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --use-singularity --cluster qsub --jobs 100


If you are using Compute Canada, you can use the [cc-slurm](https://github.com/khanlab/cc-slurm) profile, which submits jobs and takes care of requesting the correct resources per job (including GPUs). Once it is set-up with cookiecutter, run:

    snakemake --profile cc-slurm

Or, with [neuroglia-helpers](https://github.com/khanlab/neuroglia-helpers) can get a 1-GPU, 8-core, 32gb node and run locally there. First, get a node with a GPU (default 8-core, 32gb, 3 hour limit):

    regularInteractive -g
    
Then, run:

    snakemake --use-singularity --cores 8 --resources gpu=1 mem=32000


See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

### Step 5: Investigate results

After successful execution, you can create a self-contained interactive HTML report with all results via:

    snakemake --report report.html

This report can, e.g., be forwarded to your collaborators.
An example (using some trivial test data) can be seen [here](https://cdn.rawgit.com/snakemake-workflows/rna-seq-kallisto-sleuth/master/.test/report.html).

### Step 6: Commit changes

Whenever you change something, don't forget to commit the changes back to your github copy of the repository:

    git commit -a
    git push

### Step 7: Obtain updates from upstream

Whenever you want to synchronize your workflow copy with new developments from upstream, do the following.

1. Once, register the upstream repository in your local copy: `git remote add -f upstream git@github.com:snakemake-workflows/{{cookiecutter.repo_name}}.git` or `git remote add -f upstream https://github.com/snakemake-workflows/{{cookiecutter.repo_name}}.git` if you do not have setup ssh keys.
2. Update the upstream version: `git fetch upstream`.
3. Create a diff with the current version: `git diff HEAD upstream/master workflow > upstream-changes.diff`.
4. Investigate the changes: `vim upstream-changes.diff`.
5. Apply the modified diff via: `git apply upstream-changes.diff`.
6. Carefully check whether you need to update the config files: `git diff HEAD upstream/master config`. If so, do it manually, and only where necessary, since you would otherwise likely overwrite your settings and samples.


### Step 8: Contribute back

In case you have also changed or added steps, please consider contributing them back to the original repository:

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the original repo to a personal or lab account.
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the fork to your local system, to a different place than where you ran your analysis.
3. Copy the modified files from your analysis to the clone of your fork, e.g., `cp -r workflow path/to/fork`. Make sure to **not** accidentally copy config file contents or sample sheets. Instead, manually update the example config files if necessary.
4. Commit and push your changes to your fork.
5. Create a [pull request](https://help.github.com/en/articles/creating-a-pull-request) against the original repository.

## Testing

TODO: create some test datasets 


