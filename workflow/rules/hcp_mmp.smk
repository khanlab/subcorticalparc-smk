from functools import partial


H_to_hemi = dict({"L": "lh", "R": "rh"})


hcp_mmp_bids = partial(
    bids,
    root="results/hcp_mmp",
    subject="{subject}",
    hemi="{hemi}",
)


wildcard_constraints:
    surfname="white|pial|sphere.reg",
    volname="T1",


rule extract_from_tar:
    input:
        tar=config["freesurfer"]["tar"],
    params:
        out_folder=config["freesurfer"]["root"],
        file_in_tar="{subject}/{modality}/{filename}",
    output:
        filename=join(
            config["freesurfer"]["subject"],
            "{modality,surf|mri}",
            "{filename}",
        ),
    threads: 8
    shadow:
        "minimal"
    group:
        "participant1"
    shell:
        "mkdir -p {params.out_folder} && tar --extract --file={input.tar} {params.file_in_tar} && "
        "mv {params.file_in_tar} {output.filename}"


def get_gifti_input(wildcards):
    #    if wildcards.surfname == 'pial': #add .T1 to the name (since pial is a symlink to pial.T1) so can use if extracting from tar
    #        return join(config['in_freesurfer'],'surf','{hemi}.{surfname}.T1'.format(hemi=H_to_hemi[wildcards.hemi],surfname=wildcards.surfname))
    #    else:
    #        return join(config['in_freesurfer'],'surf','{hemi}.{surfname}'.format(hemi=H_to_hemi[wildcards.hemi],surfname=wildcards.surfname))
    return join(
        config["freesurfer"]["subject"],
        "surf",
        "{hemi}.{surfname}".format(
            hemi=H_to_hemi[wildcards.hemi], surfname=wildcards.surfname
        ),
    )


rule convert_to_gifti:
    input:
        get_gifti_input,
    output:
        hcp_mmp_bids(
            space="fsaverage",
            suffix="{surfname}.surf.gii",
        ),
    params:
        license=config["fs_license"],
    container:
        config["singularity"]["freesurfer"]
    log:
        "logs/convert_to_gifti/sub-{subject}_{hemi}_{surfname}.log",
    group:
        "participant1"
    shell:
        "FS_LICENSE={params.license} mris_convert {input} {output} &> {log}"


rule convert_to_nifti:
    input:
        join(config["freesurfer"]["subject"], "mri", "{volname}.mgz"),
    output:
        bids(
            root="results/hcp_mmp",
            subject="{subject}",
            suffix="{volname}.nii.gz",
        ),
    params:
        license=config["fs_license"],
    container:
        config["singularity"]["freesurfer"]
    log:
        "logs/convert_to_nifti/sub-{subject}_{volname}.log",
    group:
        "participant1"
    shell:
        "FS_LICENSE={params.license} mri_convert {input} {output} &> {log}"


rule get_tkr2scanner:
    input:
        t1=bids(
            root="results/hcp_mmp", subject="{subject}", suffix="T1.nii.gz"
        ),
    output:
        tkr2scanner=bids(
            root="results/hcp_mmp",
            subject="{subject}",
            suffix="tkr2scanner.xfm",
        ),
    params:
        license=config["fs_license"],
    container:
        config["singularity"]["freesurfer"]
    log:
        "logs/get_tkr2scanner/sub-{subject}.log",
    group:
        "participant1"
    shell:
        "FS_LICENSE={params.license} mri_info {input.t1} --tkr2scanner > {output.tkr2scanner} 2> {log}"


rule apply_surf_tkr2scanner:
    input:
        surf=rules.convert_to_gifti.output,
        tkr2scanner=rules.get_tkr2scanner.output.tkr2scanner,
    output:
        surf=hcp_mmp_bids(
            space="native",
            suffix="{surfname}.surf.gii",
        ),
    threads: 8
    container:
        config["singularity"]["connectome_workbench"]
    log:
        "logs/apply_surf_tkr2scanner/sub-{subject}_{hemi}_{surfname}.log",
    group:
        "participant1"
    shell:
        "wb_command -surface-apply-affine {input.surf} {input.tkr2scanner} {output.surf} &> {log}"


rule gen_midthickness:
    input:
        white=hcp_mmp_bids(
            space="{space}",
            suffix="white.surf.gii",
        ),
        pial=hcp_mmp_bids(
            space="{space}",
            suffix="pial.surf.gii",
        ),
    output:
        midthickness=hcp_mmp_bids(
            space="{space}",
            suffix="midthickness.surf.gii",
        ),
    container:
        config["singularity"]["connectome_workbench"]
    threads: 8
    log:
        "logs/gen_midthickness/sub-{subject}_{hemi}_{space}.log",
    group:
        "participant1"
    shell:
        "wb_command -surface-average {output.midthickness} -surf {input.white} -surf {input.pial} &> {log}"


rule resample_subj_to_fsaverage_sphere:
    input:
        surf=hcp_mmp_bids(
            space="fsaverage",
            suffix="midthickness.surf.gii",
        ),
        current_sphere=hcp_mmp_bids(
            space="fsaverage",
            suffix="sphere.reg.surf.gii",
        ),
        new_sphere="resources/standard_mesh_atlases/resample_fsaverage/"
        "fs_LR-deformed_to-fsaverage.{hemi}.sphere.32k_fs_LR.surf.gii",
    params:
        method="BARYCENTRIC",
    output:
        surf=hcp_mmp_bids(
            space="fsLR",
            den="32k",
            suffix="midthickness.surf.gii",
        ),
    container:
        config["singularity"]["connectome_workbench"]
    threads: 8
    log:
        "logs/resample_subj_to_fsaverage_sphere/sub-{subject}_{hemi}.log",
    group:
        "participant1"
    shell:
        "wb_command -surface-resample {input.surf} {input.current_sphere} {input.new_sphere} {params.method} {output.surf} &> {log}"


rule resample_labels_to_subj_sphere:
    input:
        label=lambda wildcards: "resources/standard_mesh_atlases/{hemi}.hcp-mmp.32k_fs_LR.label.gii".format(
            hemi=H_to_hemi[wildcards.hemi]
        ),
        current_sphere=lambda wildcards: "resources/standard_mesh_atlases/resample_fsaverage/"
        "fs_LR-deformed_to-fsaverage.{hemi}.sphere.32k_fs_LR.surf.gii",
        new_sphere=hcp_mmp_bids(
            space="fsaverage",
            suffix="sphere.reg.surf.gii",
        ),
        current_surf=rules.resample_subj_to_fsaverage_sphere.output.surf,
        new_surf=hcp_mmp_bids(
            space="fsaverage",
            suffix="midthickness.surf.gii",
        ),
    params:
        method="ADAP_BARY_AREA",
    output:
        label=hcp_mmp_bids(
            label="hcpmmp",
            space="native",
            suffix="dseg.label.gii",
        ),
    container:
        config["singularity"]["connectome_workbench"]
    threads: 8
    log:
        "logs/resample_labels_to_subj_sphere/sub-{subject}_{hemi}.log",
    group:
        "participant1"
    shell:
        "wb_command -label-resample {input.label} {input.current_sphere} {input.new_sphere}"
        " {params.method} {output.label}"
        " -area-surfs {input.current_surf} {input.new_surf} &> {log}"


rule map_labels_to_volume_ribbon:
    input:
        label=rules.resample_labels_to_subj_sphere.output.label,
        surf=hcp_mmp_bids(
            space="native",
            suffix="midthickness.surf.gii",
        ),
        vol_ref=bids(
            root="results/hcp_mmp", subject="{subject}", suffix="T1.nii.gz"
        ),
        white_surf=hcp_mmp_bids(
            space="native",
            suffix="white.surf.gii",
        ),
        pial_surf=hcp_mmp_bids(
            space="native",
            suffix="pial.surf.gii",
        ),
    output:
        label_vol=hcp_mmp_bids(
            label="hcpmmp",
            space="native",
            suffix="dseg.nii.gz",
        ),
    container:
        config["singularity"]["connectome_workbench"]
    threads: 8
    log:
        "logs/map_labels_to_volume_ribbon/sub-{subject}_{hemi}.log",
    group:
        "participant1"
    shell:
        "wb_command -label-to-volume-mapping {input.label} {input.surf} {input.vol_ref} {output.label_vol}"
        " -ribbon-constrained {input.white_surf} {input.pial_surf}"
        " -greedy &> {log}"


# currently optional
rule map_labels_to_volume_wmboundary:
    input:
        label=rules.resample_labels_to_subj_sphere.output.label,
        surf=hcp_mmp_bids(
            suffix="white.surf.gii",
            space="native",
        ),
        vol_ref=bids(
            root="results/hcp_mmp", subject="{subject}", suffix="T1.nii.gz"
        ),
        white_surf=hcp_mmp_bids(
            suffix="white.surf.gii",
            space="native",
        ),
        pial_surf=hcp_mmp_bids(
            suffix="pial.surf.gii",
            space="native",
        ),
    params:
        nearest_vertex="{wmbdy}",
    output:
        label_vol=hcp_mmp_bids(
            space="native",
            label="hcpmmp",
            desc="wmbound{wmbdy}",
            suffix="dseg.nii.gz",
        ),
    container:
        config["singularity"]["connectome_workbench"]
    threads: 8
    log:
        "logs/map_labels_to_volume_wmboundary/sub-{subject}_{hemi}_wmbound-{wmbdy}.log",
    group:
        "participant1"
    shell:
        "wb_command -label-to-volume-mapping {input.label} {input.surf} {input.vol_ref} {output.label_vol}"
        " -nearest-vertex {params.nearest_vertex} &> {log}"
