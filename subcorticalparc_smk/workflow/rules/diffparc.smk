# Directories
hcp_mmp_dir = str(Path(config["output_dir"]) / "hcp_mmp")
diffparc_dir = str(Path(config["output_dir"]) / "diffparc")

Path(hcp_mmp_dir).mkdir(parents=True, exist_ok=True)
Path(diffparc_dir).mkdir(parents=True, exist_ok=True)

# BIDS partials
bids_hcpmmp = partial(
    bids,
    root=hcp_mmp_dir,
    datatype="labels",
    label="hcpmmp",
    space="native",
    suffix="dseg.nii.gz",
    **inputs["T1w"].input_wildcards,
)

bids_subj_diffparc = partial(
    bids,
    root=diffparc_dir,
    space="individual",
    **inputs["T1w"].input_wildcards,
)

bids_gen_diffparc = partial(
    bids,
    root=diffparc_dir,
    label="{seed}",
    **inputs["T1w"].input_wildcards,
)

bids_tpl_diffparc = partial(
    bids_template,
    root=diffparc_dir,
    template="{template}",
    hemi="{hemi}",
    label="{seed}",
)

#########################################################
# added some hints for eventual bids-derivatives naming #
# (e.g. space, label, type(dseg, mask, probseg)..)      #
#########################################################


rule probseg_to_binary_template_seed:
    input:
        seed=join(config["seed"]["dir"], config["seed"]["nii"]),
    params:
        thresh=config["prob_seg_threshold"] * 100,
    output:
        mask=bids_tpl_diffparc(
            suffix="mask.nii.gz",
        ),
    log:
        "logs/diffparc/probseg_to_binary_template_seed/binary_{template}_hemi-{hemi}_{seed}.log",
    group:
        "group0"
    threads: 8
    container:
        config["singularity"]["neuroglia"]
    shell:
        "fslmaths {input.seed} -thr {params.thresh} -bin {output.mask} &> {log}"


rule dilate_seed:
    input:
        mask=rules.probseg_to_binary_template_seed.output.mask,
    output:
        mask=bids_tpl_diffparc(
            desc="dilatedsmoothed", 
            suffix="mask.nii.gz",
        ),
    container:
        config["singularity"]["neuroglia"]
    group:
        "group0"
    shell:
        "c3d {input.mask} -dilate 1 3x3x3vox -o {output.mask}"


# space-T1w (native), dseg
rule combine_lr_hcp:
    input:
        lh=bids_hcpmmp(hemi="L"),
        rh=bids_hcpmmp(hemi="R"),
    output:
        lh_rh=bids_subj_diffparc(
            label=config["target"]["cortical"]["atlas"],
            suffix="dseg.nii.gz",
        ),
    container:
        config["singularity"]["neuroglia"]
    log:
        "logs/diffparc/combine_lr_hcp/{subject}.log",
    group:
        "participant1"
    shell:
        "fslmaths {input.lh} -max {input.rh} {output.lh_rh} &> {log}"


# transform probabilistic seed to subject
# space-T1w,  probseg
rule seed_to_subject:
    input:
        seed=rules.dilate_seed.output.mask,
        ref=join(config["bedpost"]["dir"], config["bedpost"]["ref"]),
        invwarp=config["transforms"]["ants_invwarp"],
    output:
        seed=bids_subj_diffparc(
            hemi="{hemi}",
            label="{seed}",
            from_="{template}",
            suffix="mask.nii.gz",
        )
    container:
        config["singularity"]["neuroglia"]
    log:
        "logs/diffparc/seed_to_subject/{template}_sub-{subject}_hemi-{hemi}_{seed}.log",
    group:
        "participant1"
    threads: 8
    shell:
        "antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {input.seed} -o {output} -r {input.ref}  -t {input.invwarp} &> {log}"


# create brainmask from bedpost data, and resample to chosen resolution
# space-T1w res-? mask
rule resample_brainmask:
    input:
        dwi=join(config["bedpost"]["dir"], config["bedpost"]["samples"]),
    params:
        seed_resolution=config["probtrack"]["seed_resolution"],
    output:
        mask=bids_gen_diffparc(
            label="brain",
            suffix="mask.nii.gz",
        ),
        mask_res=bids_gen_diffparc(
            label="brain",
            res="dwi",
            suffix="mask.nii.gz",
        ),
    container:
        config["singularity"]["neuroglia"]
    log:
        "logs/diffparc/resample_brainmask/sub-{subject}.log",
    group:
        "participant1"
    shell:
        "fslmaths {input.dwi} -bin {output.mask} &&"
        "mri_convert {output.mask} -vs {params.seed_resolution} {params.seed_resolution} {params.seed_resolution} {output.mask_res} -rt nearest &> {log}"


# resample target segs to the match resampled dwi mask
# space-T1w res-? dseg
rule resample_targets:
    input:
        mask_res=rules.resample_brainmask.output.mask_res,
        targets=rules.combine_lr_hcp.output.lh_rh,
    output:
        targets_res=bids_subj_diffparc(
            label=config["target"]["cortical"]["atlas"],
            res="dwi",
            suffix="dseg.nii.gz",
        ),
    container:
        config["singularity"]["neuroglia"]
    log:
        "logs/diffparc/resample_targets/sub-{subject}.log",
    group:
        "participant1"
    shell:
        "reg_resample -flo {input.targets} -res {output.targets_res} -ref {input.mask_res} -inter 0  &> {log}"


# resamples seed seg to match resampled dwi mask
# space-T1w res=? probseg
rule resample_seed:
    input:
        mask_res=rules.resample_brainmask.output.mask_res,
        seed=rules.seed_to_subject.output.seed,
    output:
        seed_res=bids_subj_diffparc(
            hemi="{hemi}",
            label="{seed}",
            from_="{template}",
            res="dwi",
            suffix="mask.nii.gz",
        ),
    container:
        config["singularity"]["neuroglia"]
    log:
        "logs/diffparc/resample_seed/{template}_sub-{subject}_hemi-{hemi}_{seed}.log",
    group:
        "participant1"
    shell:
        "reg_resample -flo {input.seed} -res {output.seed_res} -ref {input.mask_res} -inter 0 &> {log}"


# space-T1w, mask
rule split_targets:
    input:
        targets=rules.resample_targets.output.targets_res,
    params:
        target_nums=[str(i + 1) for i in range(len(targets))],
        target_seg=lambda wildcards, output: expand(
            "{target_seg_dir}/sub-{subject}_label-{target}_mask.nii.gz",
            subject=wildcards.subject,
            target=targets,
            target_seg_dir=output.target_seg_dir,
        ),
    output:
        target_seg_dir=directory(
            bids(
                root=diffparc_dir,
                suffix="targets",
                **inputs["T1w"].input_wildcards,
            )
        ),
    container:
        config["singularity"]["neuroglia"]
    log:
        "logs/diffparc/split_targets/sub-{subject}.log",
    threads: 32
    group:
        "participant1"
    shell:  #TODO: could do this in c3d with less effort.. 
        "mkdir -p {output} && parallel  --jobs {threads} fslmaths {input.targets} -thr {{1}} -uthr {{1}} -bin {{2}} &> {log} ::: {params.target_nums} :::+ {params.target_seg}"


# txt
rule gen_targets_txt:
    input:
        target_seg_dir=bids(
            root=diffparc_dir,
            suffix="targets",
            **inputs["T1w"].input_wildcards,
        )
    params:
        target_seg=lambda wildcards, input: expand(
            "{target_seg_dir}/sub-{subject}_label-{target}_mask.nii.gz",
            subject=wildcards.subject,
            target_seg_dir=input.target_seg_dir,
            target=targets,
        ),
    output:
        target_txt=join(diffparc_dir, "sub-{subject}/targets.txt"),
    log:
        "logs/diffparc/get_targets_txt/sub-{subject}.log",
    group:
        "participant1"
    run:
        f = open(output.target_txt, "w")
        for s in params.target_seg:
            f.write(f"{s}\n")
        f.close()


# probtrack dir out
rule run_probtrack:
    input:
        seed_res=rules.resample_seed.output.seed_res,
        target_txt=rules.gen_targets_txt.output.target_txt,
        mask=rules.resample_brainmask.output.mask,
        target_seg_dir=rules.gen_targets_txt.input.target_seg_dir,
    params:
        bedpost_merged=join(
            config["bedpost"]["dir"], config["bedpost"]["merged_prefix"]
        ),
        probtrack_opts=config["probtrack"]["opts"],
        out_target_seg=lambda _, output: expand(
            bids(
                root=output.probtrack_dir,
                prefix="seeds_to",
                label="{target}",
                suffix="mask.nii.gz",
                include_subject_dir=False,
            ),
            target=targets,
        ),
        nsamples=config["probtrack"]["nsamples"],
        container=config["singularity"]["fsl_cuda"],
    output:
        probtrack_dir=directory(
            bids_gen_diffparc(
                hemi="{hemi}",
                from_="{template}",
                suffix="probtrack",
            )
        ),
    threads: 32
    resources:
        mem_mb=128000,
        time=30,  #30 mins
        gpus=1,  #1 gpu
    log:
        "logs/diffparc/run_probtrack/sub-{subject}_hemi-{hemi}_label-{seed}_{template}.log",
    group:
        "participant1"
    shell:
        "mkdir -p {output.probtrack_dir} && singularity exec -e --nv {params.container} "
        "probtrackx2_gpu --samples={params.bedpost_merged}  --mask={input.mask} --seed={input.seed_res} "
        "--targetmasks={input.target_txt} --seedref={input.seed_res} --nsamples={params.nsamples} "
        "--dir={output.probtrack_dir} {params.probtrack_opts} -V 2  &> {log}"


# check bids-deriv dwi draft (not in main bids yet)
# space-{template}
rule transform_conn_to_template:
    input:
        probtrack_dir=bids_gen_diffparc(
            hemi="{hemi}",
            from_="{template}",
            suffix="probtrack",
        ),
        warp=config["transforms"]["ants_warp"],
        ref=rules.probseg_to_binary_template_seed.output.mask,
    params:
        in_connmap_3d=lambda _, input: expand(
            bids(
                root=input.probtrack_dir,
                subject="{subject}",
                prefix="seeds_to",
                label="{target}",
                suffix="mask.nii.gz",
                include_subject_dir=False,
            ),
            target=targets,
            allow_missing=True,
        ),
        out_connmap_3d=lambda _, output: expand(
            bids(
                root=output.probtrack_dir,
                subject="{subject}",
                prefix="seeds_to",
                label="{target}",
                suffix="mask.nii.gz",
                include_subject_dir=False,
            ),
            target=targets,
            allow_missing=True,
        ),
    output:
        probtrack_dir=directory(
            bids_gen_diffparc(
                hemi="{hemi}",
                space="{template}",
                suffix="probtrack",
            )
        ),
    container:
        config["singularity"]["neuroglia"]
    threads: 32
    resources:
        mem_mb=128000,
    log:
        "logs/diffparc/transform_conn_to_template/sub-{subject}_hemi-{hemi}_label-{seed}_{template}.log",
    group:
        "participant1"
    shell:
        "mkdir -p {output} && ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1 parallel  --jobs {threads} antsApplyTransforms -d 3 --interpolation Linear -i {{1}} -o {{2}}  -r {input.ref} -t {input.warp}  &> {log} :::  {params.in_connmap_3d} :::+ {params.out_connmap_3d}"


# check bids-deriv -- connectivity?
# space-{template}
rule save_connmap_template_npz:
    input:
        mask=rules.probseg_to_binary_template_seed.output.mask,
        probtrack_dir=bids_gen_diffparc(
            hemi="{hemi}",
            space="{template}",
            suffix="probtrack",
        ),
    params:
        connmap_3d=lambda wildcards, input: expand(
            bids(
                root=input.probtrack_dir,
                subject="{subject}",
                prefix="seeds_to",
                label="{target}",
                suffix="mask.nii.gz",
                include_subject_dir=False,
            ),
            target=targets,
            allow_missing=True,
        ),
    output:
        connmap_npz=bids_gen_diffparc(
            hemi="{hemi}",
            space="{template}",
            suffix="connMap.npz",
        ),
    log:
        "logs/diffparc/save_connmap_to_template_npz/sub-{subject}_hemi-{hemi}_label-{seed}_{template}.log",
    container:
        config["singularity"]["pythondeps"]
    group:
        "participant1"
    script:
        "../scripts/diffparc/save_connmap_template_npz.py"


# space-{template}
rule gather_connmap_group:
    input:
        connmap_npz=expand(
            rules.save_connmap_template_npz.output.connmap_npz,
            subject=inputs["T1w"].input_lists["subject"],
            allow_missing=True,
        ),
    output:
        connmap_group_npz=bids_tpl_diffparc(
            desc="concat",
            from_="group",
            suffix="connMap.npz",
        ),
    log:
        "logs/diffparc/gather_connmap_group/{hemi}_{seed}_{template}.log",
    container:
        config["singularity"]["pythondeps"]
    group:
        "diffparc_group1"
    script:
        "../scripts/diffparc/gather_connmap_group.py"


# space-{template},  dseg
rule spectral_clustering:
    input:
        connmap_group_npz=rules.gather_connmap_group.output.connmap_group_npz,
    params:
        max_k=config["max_k"],
    output:
        cluster_k=expand(
            bids_tpl_diffparc(
                from_="group",
                method="spectralcosine",
                k="{k}",
                suffix="dseg.nii.gz",
            ),
            k=range(2, config["max_k"] + 1),
            allow_missing=True,
        ),
    resources:
        mem_mb=64000,
    log:
        "logs/diffparc/spectral_clustering/{hemi}_{seed}_{template}.log",
    container:
        config["singularity"]["pythondeps"]
    group:
        "diffparc_group1"
    script:
        "../scripts/diffparc/spectral_clustering.py"


# sorting the cluster label by AP, sorted = desc
rule sort_cluster_label:
    input:
        seg=bids_tpl_diffparc(
            from_="group",
            method="spectralcosine",
            k="{k}",
            suffix="dseg.nii.gz",
        ),
        shell_script= Path(workflow.basedir).parent / "workflow/scripts/diffparc/sort_labels_by_ap.sh",
    output:
        seg=bids_tpl_diffparc(
            from_="group",
            method="spectralcosine",
            k="{k}",
            desc="sorted",
            suffix="dseg.nii.gz",
        ),
    container:
        config["singularity"]["neuroglia"]
    group:
        "diffparc_group1"
    log:
        "logs/diffparc/sort_cluster_label/{hemi}_{seed}_{template}_{k}.log",
    shadow:
        "minimal"  #run in shadow dir to avoid conflicts between centroids txt files
    shell:
        "{input.shell_script} {input.seg} {output.seg} {wildcards.k} &> {log}"
