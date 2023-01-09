H_to_hemi = dict({"L": "lh", "R": "rh"})

# Directories
hcp_mmp_dir = str(Path(config["output_dir"]) / "hcp_mmp")
fs_dir = str(Path(config["output_dir"]) / "freesurfer")

# BIDS
bids_hcpmmp = partial(
    bids,
    root=hcp_mmp_dir,
    datatype="surf",
    hemi="{hemi}",
    **inputs.input_wildcards['T1w'],
)

bids_fs = partial(
    bids,
    root=fs_dir
    **inputs.input_wildcards['T1w'],
)

wildcard_constraints:
    surfname="white|pial|sphere.reg",


# Diffparc
rule extract_from_tar:
    input:
        tar=config["freesurfer"]["tar"],
    params:
        out_folder=fs_dir,
        file_in_tar="{subject}/{modality}/{filename}",
    output:
        filename=bids_freesurfer(
            datatype="{modality,surf|mri}",
            suffix="{filename},
        )
    threads: 8
    shadow:
        "minimal"
    group:
        "diffparc_participant1"
    shell:
        "mkdir -p {params.out_folder} && tar --extract --file={input.tar} {params.file_in_tar} && "
        "mv {params.file_in_tar} {output.filename}"


def get_gifti_input(wildcards):
    #    if wildcards.surfname == 'pial': #add .T1 to the name (since pial is a symlink to pial.T1) so can use if extracting from tar
    #        return join(config['in_freesurfer'],'surf','{hemi}.{surfname}.T1'.format(hemi=H_to_hemi[wildcards.hemi],surfname=wildcards.surfname))
    #    else:
    #        return join(config['in_freesurfer'],'surf','{hemi}.{surfname}'.format(hemi=H_to_hemi[wildcards.hemi],surfname=wildcards.surfname))
    return bids_freesurfer(
        datatype="surf",
        suffix="{hemi}.{surfname}".format(
            hemi=H_to_hemi[wildcards.hemi], 
            surfname=wildcards.surfname,
        )
    )


rule convert_to_gifti:
    input:
        get_gifti_input,
    output:
        bids_hcpmmp(
            space="fsaverage",
            suffix="{surfname}.surf.gii",
        )
    params:
        license=config["fs_license"],
    container:
        config["singularity"]["freesurfer"]
    log:
        "logs/diffparc/convert_to_gifti/sub-{subject}_{hemi}_{surfname}.log",
    group:
        "diffparc_participant1"
    shell:
        "FS_LICENSE={params.license} mris_convert {input} {output} &> {log}"


rule convert_to_nifti:
    input:
        bids_fs(
            datatype="mri",
            suffix="T1.mgz",
        ),
    output:
        bids_hcpmmp(
            datatype="anat",
            suffix="T1.nii.gz",
        )
    params:
        license=config["fs_license"],
    container:
        config["singularity"]["freesurfer"]
    log:
        "logs/diffparc/convert_to_nifti/sub-{subject}_T1.log",
    group:
        "diffparc_participant1"
    shell:
        "FS_LICENSE={params.license} mri_convert {input} {output} &> {log}"


rule get_tkr2scanner:
    input:
        rules.convert_to_nifti.output,
    output:
        bids_hcpmmp(
            datatype="xfm",
            suffix="tkr2scanner.xfm",
        )
    params:
        license=config["fs_license"],
    container:
        config["singularity"]["freesurfer"]
    log:
        "logs/diffparc/get_tkr2scanner/sub-{subject}.log",
    group:
        "diffparc_participant1"
    shell:
        "FS_LICENSE={params.license} mri_info {input} --tkr2scanner > {output} 2> {log}"


rule apply_surf_tkr2scanner:
    input:
        surf=rules.convert_to_gifti.output,
        tkr2scanner=rules.get_tkr2scanner.output,
    output:
        surf=bids_hcpmmp(
            space="native",
            suffix="{surfname}.surf.gii",
        ),
    threads: 8
    container:
        config["singularity"]["connectome_workbench"]
    log:
        "logs/diffparc/apply_surf_tkr2scanner/sub-{subject}_{hemi}_{surfname}.log",
    group:
        "diffparc_participant1"
    shell:
        "wb_command -surface-apply-affine {input.surf} {input.tkr2scanner} {output.surf} &> {log}"


rule gen_midthickness:
    input:
        white=bids_hcpmmp(
            space="{space}",
            suffix="white.surf.gii",
        )
        pial=bids_hcpmmp(
            space="{space}",
            suffix="pial.surf.gii",
        )
    output:
        midthickness=bids_hcpmmp(
            space="{space}",
            suffix="midthickness.surf.gii",
        ),
    container:
        config["singularity"]["connectome_workbench"]
    threads: 8
    log:
        "logs/diffparc/gen_midthickness/sub-{subject}_{hemi}_{space}.log",
    group:
        "diffparc_participant1"
    shell:
        "wb_command -surface-average {output.midthickness} -surf {input.white} -surf {input.pial} &> {log}"


rule resample_subj_to_fsaverage_sphere:
    input:
        surf=bids_hcpmmp(
            space="fsaverage",
            suffix="midthickness.surf.gii",
        ),
        current_sphere=bids_hcpmmp(
            space="fsaverage",
            suffix="sphere.reg.surf.gii",
        ),
        new_sphere="resources/standard_mesh_atlases/resample_fsaverage/"
        "fs_LR-deformed_to-fsaverage.{hemi}.sphere.32k_fs_LR.surf.gii",
    params:
        method="BARYCENTRIC",
    output:
        surf=bids_hcpmmp(
            space="fsLR",
            den="32k",
            suffix="midthickness.surf.gii",
        ),
    container:
        config["singularity"]["connectome_workbench"]
    threads: 8
    log:
        "logs/diffparc/resample_subj_to_fsaverage_sphere/sub-{subject}_{hemi}.log",
    group:
        "diffparc_participant1"
    shell:
        "wb_command -surface-resample {input.surf} {input.current_sphere} {input.new_sphere} {params.method} {output.surf} &> {log}"


rule resample_labels_to_subj_sphere:
    input:
        label=lambda wildcards: "resources/standard_mesh_atlases/{hemi}.hcp-mmp.32k_fs_LR.label.gii".format(
            hemi=H_to_hemi[wildcards.hemi]
        ),
        current_sphere=lambda wildcards: "resources/standard_mesh_atlases/resample_fsaverage/"
        "fs_LR-deformed_to-fsaverage.{hemi}.sphere.32k_fs_LR.surf.gii",
        current_surf=rules.resample_subj_to_fsaverage_sphere.output.surf,
        new_sphere=rules.resample_subj_to_fsaverage_sphere.input.current_sphere,
        new_surf=rules.resample_subj_to_fsaverage_sphere.input.current_surf,
    params:
        method="ADAP_BARY_AREA",
    output:
        label=bids_hcpmmp(
            datatype="labels",
            label="hcpmmp",
            space="native",
            suffix="dseg.label.gii",
        ),
    container:
        config["singularity"]["connectome_workbench"]
    threads: 8
    log:
        "logs/diffparc/resample_labels_to_subj_sphere/sub-{subject}_{hemi}.log",
    group:
        "diffparc_participant1"
    shell:
        "wb_command -label-resample {input.label} {input.current_sphere} {input.new_sphere}"
        " {params.method} {output.label}"
        " -area-surfs {input.current_surf} {input.new_surf} &> {log}"


rule map_labels_to_volume_ribbon:
    input:
        label=rules.resample_labels_to_subj_sphere.output.label,
        vol_ref=rules.convert_to_nifti.output,
        surf=bids_hcpmmp(
            space="native",
            suffix="midthickness.surf.gii",
        ),
        white_surf=bids_hcpmmp(
            space="native",
            suffix="white.surf.gii",
        ),
        pial_surf=bids_hcpmmp(
            space="native",
            suffix="pial.surf.gii",
        ),
    output:
        label_vol=bids_hcpmmp(
            datatype="labels",
            label="hcpmmp",
            space="native",
            suffix="dseg.nii.gz",
        ),
    container:
        config["singularity"]["connectome_workbench"]
    threads: 8
    log:
        "logs/diffparc/map_labels_to_volume_ribbon/sub-{subject}_{hemi}.log",
    group:
        "diffparc_participant1"
    shell:
        "wb_command -label-to-volume-mapping {input.label} {input.surf} {input.vol_ref} {output.label_vol}"
        " -ribbon-constrained {input.white_surf} {input.pial_surf}"
        " -greedy &> {log}"


# currently optional
rule map_labels_to_volume_wmboundary:
    input:
        label=rules.resample_labels_to_subj_sphere.output.label,
        surf=rules.map_labels_to_volume_ribbon.output.white_surf,
        vol_ref=rules.convert_to_nifti.output,
        white_surf=rules.map_labels_to_volume_ribbon.output.white_surf,
        pial_surf=rules.map_labels_to_volume_ribbon.output.pial_surf,
    params:
        nearest_vertex="{wmbdy}",
    output:
        label_vol=bids_hcpmmp(
            datatype="labels",
            space="native",
            label="hcpmmp",
            desc="wmbound{wmbdy}",
            suffix="dseg.nii.gz",
        ),
    container:
        config["singularity"]["connectome_workbench"]
    threads: 8
    log:
        "logs/diffparc/map_labels_to_volume_wmboundary/sub-{subject}_{hemi}_wmbound-{wmbdy}.log",
    group:
        "diffparc_participant1"
    shell:
        "wb_command -label-to-volume-mapping {input.label} {input.surf} {input.vol_ref} {output.label_vol}"
        " -nearest-vertex {params.nearest_vertex} &> {log}"

# TODO
# Funcparc
rule extract_from_zip:
    input:
        packages=expand(
            join(config["hcp_func"]["dir"], "{subject}_{package}.zip"),
            package=config["hcp_func"]["icafix_package_dict"].values(),
            allow_missing=True,
        )
    params:
        files_in_pkg=expand(
            "{filename}",
            filename=config['hcp_func']['icafix_package_dict'].keys(),    
        ),
    output:
        files=expand(
            "results/hcp_func/data/{filename}",
            filename=config['hcp_func']['icafix_package_dict'].keys(),
        )
    group: "funcparc_participant1"
    run:
        for pkg, f in zip(input.packages, params.files_in_pkg):
            shell("unzip -n {pkg} {f} -d results/hcp_func/data")


rule merge_roi:
    '''Create custom subcortical atlas, labelling seed as 16 (brainstem) for wb'''
    input:
        atlas = "resources/funcparc/sub-MNI152NLin6Asym_desc-allFSstyle_workbench.nii.gz",
        roi = "resources/funcparc/sub-SNSX32NLin2020Asym_space-MNI152NLin6Asym_hemi-LR_desc-ZIR_res-1p6mm_mask.nii.gz",
    output:
        out_fpath = "results/hcp_func/group/atlas/seed-{seed}_BigBrain{vox_res}.nii.gz",
    container: 
        config["singularity"]["pythondeps"]
    group: 
        "funcparc_group1"
    script: 
        "../scripts/funcparc/merge_rois.py"


rule create_atlas:
    '''Combine ZIR and BigBrain subcortical labels'''
    input:
        atlas = rules.merge_roi.output.out_fpath,
        labels = "resources/funcparc/sub-MNI152NLin6Asym_desc-allFSstyle_labels.txt",
        lh_mmp = "resources/standard_mesh_atlases/lh.hcp-mmp.59k_fs_LR.label.gii",
        rh_mmp = "resources/standard_mesh_atlases/rh.hcp-mmp.59k_fs_LR.label.gii",
    output:
        nifti = 'results/hcp_func/group/atlas/seed-{seed}_HCP_MMP_BigBrain{vox_res}.nii.gz',
        cifti = 'results/hcp_func/group/atlas/seed-{seed}_HCP_MMP_BigBrain{vox_res}.dlabel.nii',
    container: config['singularity']['connectome_workbench']
    group:
        "funcparc_group1"
    shell:
        "wb_command -volume-label-import {input.atlas} {input.labels} {output.nifti} && "
        "wb_command -cifti-create-label {output.cifti} -volume {output.nifti} {output.nifti} -left-label {input.lh_mmp} -right-label {input.rh_mmp}"