from functools import partial

# BIDS partials
hcp_mmp_bids = partial(
    bids,
    root="results/hcp_mmp",
    subject="{subject}",
    label="hcpmmp",
    space="native",
    suffix="dseg.nii.gz",
)

diffparc_subject = partial(
    bids,
    root="results/diffparc",
    subject="{subject}",
    space="individual",
)

diffparc_template = partial(
    bids,
    root="results/diffparc",
    template="{template}",
    label="{seed}",
)

diffparc_probtrack = partial(
    bids,
    root="results/diffparc",
    subject="{subject}",
    label="{seed}",
)

#########################################################
# added some hints for eventual bids-derivatives naming #
# (e.g. space, label, type(dseg, mask, probseg)..)      #
#########################################################


# space-T1w (native), dseg
rule combine_lr_hcp:
    input:
        lh=hcp_mmp_bids(hemi="L"),
        rh=hcp_mmp_bids(hemi="R"),
    output:
        lh_rh=diffparc_subject(
            label=config["target"]["cortical"]["atlas"],
            suffix="dseg.nii.gz",
        ),
    container:
        config["singularity"]["neuroglia"]
    log:
        "logs/combine_lr_hcp/{subject}.log",
    group:
        "participant1"
    shell:
        "fslmaths {input.lh} -max {input.rh} {output.lh_rh} &> {log}"


# mask
rule get_binary_template_seed:
    input:
        seed=join(config["seed"]["dir"], config["seed"]["nii"]),
    params:
        seg_thresh=config["prob_seg_threshold"],
    output:
        mask=diffparc_template(
            desc=f"thresh{str(config['prob_seg_threshold'])}",
            hemi="{hemi}",
            suffix="mask.nii.gz",
        ),
    container:
        config["singularity"]["neuroglia"]
    group:
        "group0"
    shell:
        "fslmaths {input} -thr {params.seg_thresh} {output}"


rule dilate_seed:
    input:
        mask=rules.get_binary_template_seed.output.mask,
    output:
        seed=diffparc_template(
            hemi="{hemi}",
            desc="dilatedsmoothed",
            suffix="mask.nii.gz",
        ),
    container:
        config["singularity"]["neuroglia"]
    group:
        "group0"
    shell:
        "c3d {input} -dilate 1 3x3x3vox -o {output}"


# transform probabilistic seed to subject
# space-T1w,  probseg
rule transform_to_subject:
    input:
        seed=rules.dilate_seed.output.seed,
        ref=rules.combine_lr_hcp.output.lh_rh,
        invwarp=config["transforms"]["ants_invwarp"],
    output:
        seed=diffparc_subject(
            hemi="{hemi}",
            label="{seed}",
            from_="{template}",
            suffix="mask.nii.gz",
        ),
    container:
        config["singularity"]["neuroglia"]
    log:
        "logs/transform_to_subject/{template}_sub-{subject}_hemi-{hemi}_{seed}.log",
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
        mask=bids(
            root="results/diffparc",
            subject="{subject}",
            label="brain",
            suffix="mask.nii.gz",
        ),
        mask_res=bids(
            root="results/diffparc",
            subject="{subject}",
            label="brain",
            res="dwi",
            suffix="mask.nii.gz",
        ),
    container:
        config["singularity"]["neuroglia"]
    log:
        "logs/resample_brainmask/sub-{subject}.log",
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
        targets_res=diffparc_subject(
            label=config["target"]["cortical"]["atlas"],
            res="dwi",
            suffix="dseg.nii.gz",
        ),
    container:
        config["singularity"]["neuroglia"]
    log:
        "logs/resample_targets/sub-{subject}.log",
    group:
        "participant1"
    shell:
        "reg_resample -flo {input.targets} -res {output.targets_res} -ref {input.mask_res} -NN 0  &> {log}"


# resamples seed seg to match resampled dwi mask
# space-T1w res=? probseg
rule resample_seed:
    input:
        mask_res=rules.resample_brainmask.output.mask_res,
        seed=rules.transform_to_subject.output.seed,
    output:
        seed_res=diffparc_subject(
            hemi="{hemi}",
            label="{seed}",
            from_="{template}",
            res="dwi",
            suffix="mask.nii.gz",
        ),
    container:
        config["singularity"]["neuroglia"]
    log:
        "logs/resample_seed/{template}_sub-{subject}_hemi-{hemi}_{seed}.log",
    group:
        "participant1"
    shell:
        #linear interp here now, since probabilistic seg
        "reg_resample -flo {input.seed} -res {output.seed_res} -ref {input.mask_res} -NN 0 &> {log}"


# space-T1w, mask
rule split_targets:
    input:
        targets=rules.resample_targets.output.targets_res,
    params:
        target_nums=lambda wildcards: [str(i + 1) for i in range(len(targets))],
        target_seg=lambda wildcards, output: expand(
            "{target_seg_dir}/sub-{subject}_label-{target}_mask.nii.gz",
            target_seg_dir=output.target_seg_dir,
            subject=wildcards.subject,
            target=targets,
        ),
    output:
        target_seg_dir=directory(
            bids(
                root="results/diffparc",
                subject="{subject}",
                suffix="targets",
            )
        ),
    container:
        config["singularity"]["neuroglia"]
    log:
        "logs/split_targets/sub-{subject}.log",
    threads: 32
    group:
        "participant1"
    shell:  #TODO: could do this in c3d with less effort.. 
        "mkdir -p {output} && parallel  --jobs {threads} fslmaths {input.targets} -thr {{1}} -uthr {{1}} -bin {{2}} &> {log} ::: {params.target_nums} :::+ {params.target_seg}"


# txt
rule gen_targets_txt:
    input:
        target_seg_dir=rules.split_targets.output.target_seg_dir,
    params:
        target_seg=lambda wildcards, input: expand(
            "{target_seg_dir}/sub-{subject}_label-{target}_mask.nii.gz",
            target_seg_dir=input.target_seg_dir,
            subject=wildcards.subject,
            target=targets,
        ),
    output:
        target_txt="results/diffparc/sub-{subject}/targets.txt",
    log:
        "logs/get_targets_txt/sub-{subject}.log",
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
        target_txt=rules.gen_targets_txt.output,
        mask=rules.resample_brainmask.output.mask,
        target_seg_dir=rules.split_targets.output.target_seg_dir,
    params:
        bedpost_merged=join(
            config["bedpost"]["dir"],
            config["bedpost"]["merged_prefix"],
        ),
        probtrack_opts=config["probtrack"]["opts"],
        out_target_seg=lambda wildcards, output: expand(
            bids(
                root=output.probtrack_dir,
                include_subject_dir=False,
                prefix="seeds_to",
                hemi=wildcards.hemi,
                label="{target}",
                suffix="mask.nii.gz",
            ),
            target=targets,
        ),
        nsamples=config["probtrack"]["nsamples"],
        container=config["singularity"]["fsl_cuda"],
    output:
        probtrack_dir=directory(
            diffparc_probtrack(
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
        "logs/run_probtrack/sub-{subject}_hemi-{hemi}_{seed}_{template}.log",
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
        probtrack_dir=rules.run_probtrack.output.probtrack_dir,
        ref=rules.get_binary_template_seed.output.mask,
        warp=config["transforms"]["ants_warp"],
    params:
        in_connmap_3d=lambda wildcards, input: expand(
            bids(
                root=input.probtrack_dir,
                include_subject_dir=False,
                subject=wildcards.subject,
                prefix="seeds_to",
                label="{target}",
                suffix="mask.nii.gz",
            ),
            target=targets,
        ),
        out_connmap_3d=lambda wildcards, output: expand(
            bids(
                root=output.probtrack_dir,
                include_subject_dir=False,
                subject=wildcards.subject,
                prefix="seeds_to",
                label="{target}",
                suffix="mask.nii.gz",
            ),
            target=targets,
        ),
    output:
        probtrack_dir=directory(
            diffparc_probtrack(
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
        "logs/transform_conn_to_template/sub-{subject}_hemi-{hemi}_{seed}_{template}.log",
    group:
        "participant1"
    shell:
        "mkdir -p {output.probtrack_dir} && ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1 parallel  --jobs {threads} antsApplyTransforms -d 3 --interpolation Linear -i {{1}} -o {{2}}  -r {input.ref} -t {input.warp}  &> {log} :::  {params.in_connmap_3d} :::+ {params.out_connmap_3d}"


# check bids-deriv -- connectivity?
# space-{template}
rule save_connmap_template_npz:
    input:
        mask=rules.get_binary_template_seed.output.mask,
        probtrack_dir=rules.transform_conn_to_template.output.probtrack_dir,
    params:
        connmap_3d=lambda wildcards, input: expand(
            bids(
                root=input.probtrack_dir,
                include_subject_dir=False,
                subject=wildcards.subject,
                prefix="seeds_to",
                label="{target}",
                suffix="mask.nii.gz",
            ),
            target=targets,
        ),
    output:
        connmap_npz=diffparc_subject(
            label="{seed}",
            space="{template}",
            hemi="{hemi}",
            suffix="connMap.npz",
        ),
    log:
        "logs/save_connmap_to_template_npz/sub-{subject}_hemi-{hemi}_{seed}_{template}.log",
    group:
        "participant1"
    script:
        "../scripts/save_connmap_template_npz.py"


# space-{template}
rule gather_connmap_group:
    input:
        connmap_npz=expand(
            rules.save_connmap_template_npz.output.connmap_npz,
            subject=subjects,
            allow_missing=True,
        ),
    output:
        connmap_group_npz=diffparc_template(
            desc="concat",
            hemi="{hemi}",
            from_="group",
            suffix="connMap.npz",
        ),
    log:
        "logs/gather_connmap_group/hemi-{hemi}_{seed}_{template}.log",
    group:
        "group1"
    script:
        "../scripts/gather_connmap_group.py"


# space-{template},  dseg
rule spectral_clustering:
    input:
        connmap_group_npz=rules.gather_connmap_group.output.connmap_group_npz,
    params:
        max_k=config["max_k"],
    output:
        cluster_k=expand(
            diffparc_template(
                hemi="{hemi}",
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
        "logs/spectral_clustering/hemi-{hemi}_{seed}_{template}.log",
    group:
        "group1"
    script:
        "../scripts/spectral_clustering.py"


# sorting the cluster label by AP, sorted = desc
rule sort_cluster_label:
    input:
        seg=diffparc_template(
            hemi="{hemi}",
            from_="group",
            method="spectralcosine",
            k="{k}",
            suffix="dseg.nii.gz",
        ),
        shell_script="workflow/scripts/sort_labels_by_ap.sh",
    output:
        seg=diffparc_template(
            hemi="{hemi}",
            from_="group",
            method="spectralcosine",
            k="{k}",
            desc="sorted",
            suffix="dseg.nii.gz",
        ),
    container:
        config["singularity"]["neuroglia"]
    group:
        "group1"
    log:
        "logs/sort_cluster_label/hemi-{hemi}_{seed}_{template}_{k}.log",
    shadow:
        "minimal"  #run in shadow dir to avoid conflicts between centroids txt files
    shell:
        "{input.shell_script} {input.seg} {output.seg} {wildcards.k} &> {log}"
