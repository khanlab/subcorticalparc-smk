# Tract map processing helper scripts

## Pre-requisites:
This requires `kslurm` (https://github.com/pvandyken/kslurm) to be installed. You can do this by pasting this at a terminal:
```
curl -sSL https://raw.githubusercontent.com/pvandyken/kslurm/master/install_kslurm.py | python -
```

## Tract map processing:
Tract map processing is best done using local scratch (`$SLURM_TMPDIR`) since it is very I/O intensive. Since there isn't an easy way to achieve this built-in to snakemake, another job script has been created that copies the data to local scratch, runs snakemake to generate a single subject's tractmaps, then copies that data back before the job finishes. 


*Before* you run this job, however, you want to make sure:
 1. all the processing up to tract maps is complete. Can do this with:
 ```
  #this command does a dry-run to see if any jobs need to be run still
  snakemake  --omit-from resample_clus_seed -npr 

  #this will run all the jobs to complete processing before tract maps 
  snakemake --omit-from resample_clus_seed --profile cc-slurm
 ```
 2. the `results/tractmap` folder for a subject should be empty submitting any new tractmap jobs. 


To run tract maps for one subject (e.g. sub-100307) you can use:
```
kbatch 3:00 8 32G gpu -a ctb-akhanf_gpu . ./job_tractmaps sub-100307
```

To run on all subjects in a participants.tsv file, use:
```
./submit_tractmaps config/hcptestparticipants.tsv
```

