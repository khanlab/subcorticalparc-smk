# Directories
hcp_func_dir = str(Path(config["output_dir"]) / "hcp_func")

# BIDS
bids_hcpfunc = partial(
    bids,
    root=hcp_func_dir,
    **inputs['T1w'].input_wildcards,
)

bids_hcpfunc_group = partial(
    bids,
    root=hcp_func_dir,
)

rule prepare_subcortical:
    """Prep subcortical rs-fMRI as done by HCP"""
    input:
        vol=bids_hcpfunc(
            datatype="data",
            suffix="{subject}/MNINonLinear/Results/{run}/{run}_hp2000_clean.nii.gz",
            include_subject_dir=False,
            include_session_dir=False,
        ),
        rois=rules.create_atlas.output.nifti,
    params:
        sigma=config["hcp_func"]["sigma"],
        tmp=directory(
            bids_hcpfunc(datatype="tmp")
        ),
        command="../scripts/funcparc/prep_subcortical.sh",
    output:
        out=bids_hcpfunc(
            datatype="func",
            run="{run}",
            label="{seed}",
            suffix="AtlasSubcortical{vox_res}.nii.gz",
        )
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
        dtseries=bids_hcpfunc(
            datatype="data",
            suffix="{subject}/MNINonLinear/Results/{run}/{run}_Atlas_1.6mm_MSMAll_hp2000_clean.dtseries.nii",
            include_subject_dir=False,
            include_session_dir=False,
        ),
    output:
        lh=bids_hcpfunc(
            datatype="func",
            run="{run}",
            hemi="L",
            den="59k",
            desc="fs_LR",
            suffix="func.gii",
        ),
        rh=bids_hcpfunc(
            datatype="func",
            run="{run}",
            hemi="R",
            den="59k",
            desc="fs_LR",
            suffix="func.gii",
        ),
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
        dtseries=bids_hcpfunc(
            datatype="func",
            run="{run}",
            label="{seed}",
            den="59k",
            space="fs_LR",
            suffix="{vox_res}.dtseries.nii",
        )
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
        rois=bids_hcpfunc(
            datatype="data",
            suffix="{subject}/MNINonLinear/ROIs/Atlas_wmparc.1.60.nii.gz",
            include_subject_dir=False,
            include_session_dir=False,
        ),
        movreg=bids_hcpfunc(
            datatype="data",
            suffix="{subject}/MNINonLinear/Results/{run}/Movement_Regressors_dt.txt",
            include_subject_dir=False,
            include_session_dir=False,
        ),
    output:
        confounds=bids_hcpfunc(
            datatype="func",
            run="{run}",
            label="{seed}",
            desc="confounds",
            suffix="{vox_res}.tsv"
        ),
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
        cleaning=Path(workflow.basedir).parent / config["ciftify-clean"]["hcp"] #if hcp in out else config['ciftify-clean']['general']
    output:
        cleaned_dtseries=bids_hcpfunc(
            datatype="func",
            run="{run}",
            label="{seed}",
            den="59k",
            space="fsLR",
            desc="cleaned",
            suffix="{vox_res}.dtseries.nii",
        )
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
        concat_dtseries=bids_hcpfunc(
            datatype="func",
            desc="cleaned",
            task="rest",
            label="{seed}",
            den="59k",
            space="fsLR",
            suffix="{vox_res}.dtseries.nii"
        ),
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
        command="../scripts/funcparc/wishart_filter.sh",
        script="../scripts/funcparc/wishart_filter.m",
    output:
        cleaned=bids_hcpfunc(
            datatype="func",
            desc="cleaned+",
            task="rest",
            label="{seed}",
            den="59k",
            space="fsLR",
            suffix="{vox_res}.dtseries.nii",
        ),
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
        ptseries=bids_hcpfunc(
            datatype="func",
            desc="cleaned+",
            task="rest",
            label="{seed}",
            den="59k",
            space="fsLR",
            suffix="{vox_res}.ptseries.nii",
        ),
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
        correlation=bids_hcpfunc(
            datatype="clustering",
            label="{seed}",
            desc="correlationMatrix",
            suffix="{vox_res}.npz",
        ),
    group:
        "funcparc_participant2"
    script:
        "scripts/funcparc/compute_correlation.py"


rule combine_correlation:
    """Concatenate correlation matrices across subject"""
    input:
        correlation=expand(
            rules.compute_correlation.output.correlation,
            subject=inputs["T1w"].input_lists["subject"],
            allow_missing=True,
        ),
    output:
        combined_corr=bids_hcpfunc_group(
            datatype="group/clustering",
            label="{seed}",
            desc="correlationMatrix",
            suffix="{vox_res}.npz",
        )
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
            bids_hcpfunc_group(
                datatype="group/clustering",
                label="{seed}",
                desc="spectralCosine",
                k="{k}",
                suffix="{vox_res}clusLabels.nii.gz",
            ),
            k=range(2, config["max_k"] + 1),
            allow_missing=True,
        ),
        labels=bids_hcpfunc_group(
            datatype="group/clustering",
            label="{seed}",
            suffix="{vox_res}.csv",
        ),
    group:
        "funcparc_group2"
    script:
        "scripts/funcparc/spectral_clustering.py"
