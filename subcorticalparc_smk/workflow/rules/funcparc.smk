rule prepare_subcortical:
    """Prep subcortical rs-fMRI as done by HCP"""
    input:
        vol="results/hcp_func/data/{subject}/MNINonLinear/Results/{run}/{run}_hp2000_clean.nii.gz",
        rois=rules.create_atlas.output.nifti,
    params:
        sigma=config["hcp_func"]["sigma"],
        tmp="results/hcp_func/sub-{subject}/tmp",
        command="workflow/scripts/funcparc/prep_subcortical.sh",
    output:
        out="results/hcp_func/sub-{subject}/fmri/{run}_AtlasSubcortical{vox_res}_seed-{seed}.nii.gz",
    resources:
        mem_mb=8000,
        time=30,
    log:
        "logs/sub-{subject}/prepare_subcortical_{run}_{vox_res}_seed-{seed}.log",
    container:
        config["singularity"]["connectome_workbench"]
    group:
        "funcparc_participant2"
    shell:
        "bash {params.command} {input.vol} {input.rois} {params.tmp} {output} &> {log}"


rule cifti_separate:
    """Separate cortical timeseries data from CIFTI timeseries provided by HCP"""
    input:
        dtseries="results/hcp_func/data/{subject}/MNINonLinear/Results/{run}/{run}_Atlas_1.6mm_MSMAll_hp2000_clean.dtseries.nii",
    output:
        lh="results/hcp_func/sub-{subject}/fmri/{run}.L.59k_fs_LR.func.gii",
        rh="results/hcp_func/sub-{subject}/fmri/{run}.R.59k_fs_LR.func.gii",
    container:
        config["singularity"]["connectome_workbench"]
    group:
        "funcparc_participant2"
    shell:
        "wb_command -cifti-separate {input.dtseries} COLUMN -metric CORTEX_LEFT {output.lh} -metric CORTEX_RIGHT {output.rh}"


rule create_dtseries:
    """Recreate CIFTI timeseries data using updated subcortical data"""
    input:
        vol=rules.prepare_subcortical.output.out,
        rois=rules.create_atlas.output.nifti,
        lh=rules.cifti_separate.output.lh,
        rh=rules.cifti_separate.output.rh,
    output:
        dtseries="results/hcp_func/sub-{subject}/fmri/{run}_seed-{seed}_{vox_res}.59k_fs_LR.dtseries.nii",
    container:
        config["singularity"]["connectome_workbench"]
    group:
        "funcparc_participant2"
    shell:
        "wb_command -cifti-create-dense-timeseries {output} -volume {input.vol} {input.rois} -left-metric {input.lh} -right-metric {input.rh}"


rule extract_confounds:
    """Extract cofounds for cleaning rs-fMRI data"""
    input:
        vol=rules.prepare_subcortical.input.vol,
        rois="results/hcp_func/data/{subject}/MNINonLinear/ROIs/Atlas_wmparc.1.60.nii.gz",
        movreg="results/hcp_func/data/{subject}/MNINonLinear/Results/{run}/Movement_Regressors_dt.txt",
    output:
        confounds="results/hcp_func/sub-{subject}/fmri/{run}_seed-{seed}_{vox_res}_confounds.tsv",
    log:
        "logs/sub-{subject}/extract_confounds_{run}_{vox_res}_seed-{seed}.log",
    container:
        config["singularity"]["connectome_workbench"]
    group:
        "funcparc_participant2"
    script:
        "scripts/funcparc/extract_confounds.py"


rule clean_dtseries:
    """Regress CSF and WM timeseries from GM rs-fMRI"""
    input:
        dtseries=rules.create_dtseries.output.dtseries,
        confounds=rules.extract_confounds.output.confounds,
    params:
        cleaning=config["ciftify-clean"]["hcp"],  #if hcp in out else config['ciftify-clean']['general']
    output:
        cleaned_dtseries="results/hcp_func/sub-{subject}/fmri/{run}_cleaned.seed-{seed}_{vox_res}_59k_fs_LR.dtseries.nii",
    resources:
        mem_mb=12000,
    log:
        "logs/sub-{subject}/clean_dtseries_{run}_{vox_res}_seed-{seed}.log",
    container:
        config["singularity"]["ciftify"]
    group:
        "funcparc_participant2"
    shell:
        "cifitify_clean_img --output-file={output.cleaned_dtseries} --confounds-tsv={input.confounds} --clean-config={params.cleaning} --verbose {input.dtseries} &> {log}"


rule concat_clean_runs:
    """Concatenate in time, the different runs"""
    input:
        dtseries=expand(
            rules.clean_dtseries.output.cleaned_dtseries,
            run=config["hcp_func"]["runs"],
            allow_missing=True,
        ),
    output:
        concat_dtseries="results/hcp_func/sub-{subject}/fmri/rfMRI_REST_7T_cleaned.seed-{seed}_{vox_res}_59k_fs_LR.dtseries.nii",
    group:
        "funcparc_participant2"
    container:
        config["singularity"]["connectome_workbench"]
    shell:
        "wb_command -cifti-merge {output.concat_dtseries} -cifti `echo '{input.dtseries}' | sed 's/ / -cifti /g'`"


rule wishart_postfilter:
    """As used in https://rdcu.be/b7N8K (requires matlab)
    Also performs detrending and demeaning"""
    input:
        concat_dtseries=rules.concat_clean_runs.output.concat_dtseries,
    params:
        command="workflow/scripts/funcparc/wishart_filter.sh",
        script="workflow/scripts/funcparc/wishart_filter.m",
    output:
        cleaned="results/hcp_func/sub-{subject}/fmri/rfMRI_REST_7T_cleaned+.seed-{seed}_{vox_res}_59k_fs_LR.dtseries.nii",
    group:
        "funcparc_participant2"
    shell:
        "bash {params.command} {params.script} {input} {output}"


rule parcellate_tseries:
    """Parcellate timeseries using HCP-MMP atlas"""
    input:
        dtseries=rules.wishart_postfilter.output.cleaned,
        rois=rules.create_atlas.output.cifti,
    output:
        ptseries="results/hcp_func/sub-{subject}/fmri/rfMRI_REST_7T_cleaned+.seed-{seed}_{vox_res}_59k_fs_LR.ptseries.nii",
    resources:
        mem_mb=12000,
    group:
        "funcparc_participant2"
    container:
        config["singularity"]["connectome_workbench"]
    shell:
        "wb_command -cifti-parcellate {input.dtseries} {input.rois} COLUMN {output.ptseries} -only-numeric"


rule compute_correlation:
    """Compute ZIR voxels vs cortical ROIs correlation matrix 
    per subject"""
    input:
        ptseries=rules.parcellate_tseries.output.ptseries,
        vol=rules.wishart_postfilter.output.cleaned,
    params:
        seed=config["seed"]["structures"]["ZIR"],
    output:
        correlation="results/hcp_func/sub-{subject}/clustering/correlation_matrix_seed-{seed}_{vox_res}.npz",
    group:
        "funcparc_participant2"
    script:
        "scripts/funcparc/compute_correlation.py"


rule combine_correlation:
    """Concatenate correlation matrices across subject"""
    input:
        correlation=expand(
            rules.compute_correlation.output.correlation,
            subject=df["participant_id"].str.strip("sub-").to_list(),
            allow_missing=True,
        ),
    output:
        combined_corr="results/hcp_func/group/clustering/correlation_matrix_seed-{seed}_{vox_res}.npz",
    group:
        "funcparc_group2"
    script:
        "scripts/funcparc/combine_correlation.py"


rule func_clustering:
    """Perform spectral clustering based on fucntional connectivity"""
    input:
        correlation=rules.combine_correlation.output.combined_corr,
        rois=rules.create_atlas.output.nifti,
    params:
        max_k=config["max_k"],
    output:
        niftis=expand(
            "results/hcp_func/group/clustering/seed-{seed}_{vox_res}_method-spectralcosine_k-{k}_cluslabels.nii.gz",
            k=range(2, config["max_k"] + 1),
            allow_missing=True,
        ),
        labels="results/hcp_func/group/clustering/clusterlabels_seed-{seed}_{vox_res}.csv",
    group:
        "funcparc_group2"
    script:
        "scripts/funcparc/spectral_clustering.py"
