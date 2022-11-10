rule extract_warp_from_zip:
    input:
        hcp_zip=config["transforms"]["zip"],
    params:
        out_folder=config["transforms"]["root"],
        warp_in_zip=join(
            "{subject}/MNINonLinear/xfms", config["transforms"]["warp"]
        ),
    output:
        warp_name=join(
            config["transforms"]["subject"], config["transforms"]["warp"]
        ),
    threads: 8
    shadow:
        "minimal"
    group:
        "participant1"
    shell:
        "mkdir -p {params.out_folder} && unzip {input.hcp_zip} {params.warp_in_zip} && "
        "mv {params.file_in_zip} {output.warpname}"


rule extract_invwarp_from_zip:
    input:
        hcp_zip=config["transforms"]["zip"],
    params:
        out_folder=config["transforms"]["root"],
        invwarp_in_zip=join(
            "{subject}/MNINonLinear/xfms", config["transforms"]["invwarp"]
        ),
    output:
        invwarp_name=join(
            config["transforms"]["subject"], config["transforms"]["invwarp"]
        ),
    threads: 8
    shadow:
        "minimal"
    group:
        "participant1"
    shell:
        "mkdir -p {params.out_folder} && unzip {input.hcp_zip} {params.warp_in_zip} && "
        "mv {params.file_in_zip} {output.warpname}"


# warp group-based clusters back to each subject
# space-T1w, desc-groupclus?, dseg
rule transform_clus_to_subj:
    input:
        cluster_k=expand(
            bids(
                root="results/diffparc",
                template="{template}",
                label="{seed}",
                from_="group",
                method="spectralcosine",
                k="{k}",
                desc="sorted",
                suffix="dseg.nii.gz",
            ),
            k=range(2, config["max_k"] + 1),
            allow_missing=True,
        ),
        invwarp=rules.extract_invwarp_from_zip.output.invwarp_name,
        ref=bids(
            root="results/diffparc",
            subject="{subject}",
            space="individual",
            label=config["target"]["cortical"]["atlas"],
            suffix="dseg.nii.gz",
        ),
    output:
        cluster_k=expand(
            bids(
                root="results/tractmap",
                subject="{subject}",
                space="individual",
                label="{seed}",
                method="spectralcosine",
                k="{k}",
                from_="{template}",
                desc="sorted",
                suffix="dseg.nii.gz",
            ),
            k=range(2, config["max_k"] + 1),
            allow_missing=True,
        ),
    container:
        config["singularity"]["neuroglia"]
    log:
        "logs/transform_clus_to_subject/sub-{subject}_template-{template}_{seed}.log",
    group:
        "participant2"
    threads: 8
    shell:
        "applywarp --rel -i {input.cluster_k} -r {input.ref} -w {input.invwarp} --interp nn -o {output.cluster_k}"


# create brainmask from bedpost data, and resample to chosen resolution
# space-T1w res-? mask
rule resample_brainmask_tractmaps:
    input:
        dwi=join(config["bedpost"]["dir"], config["bedpost"]["samples"]),
    params:
        seed_resolution=config["probtrack_tractmap"]["seed_resolution"],
    output:
        mask=bids(
            root="results/tractmap",
            subject="{subject}",
            label="brain",
            suffix="mask.nii.gz",
        ),
        mask_res=bids(
            root="results/tractmap",
            subject="{subject}",
            label="brain",
            res="super",
            suffix="mask.nii.gz",
        ),
    container:
        config["singularity"]["neuroglia"]
    log:
        "logs/resample_brainmask_tractmaps/sub-{subject}.log",
    group:
        "participant2"
    shell:
        "fslmaths {input.dwi} -bin {output.mask} &&"
        "mri_convert {output.mask} -vs {params.seed_resolution} {params.seed_resolution} {params.seed_resolution} {output.mask_res} -rt nearest &> {log}"


# resamples clus seg to dwi brainmask resolution
# space-T1w, res=?   dseg
rule resample_clus_seed:
    input:
        seed=bids(
            root="results/tractmap",
            subject="{subject}",
            space="individual",
            label="{seed}",
            method="spectralcosine",
            k="{k}",
            from_="{template}",
            desc="sorted",
            suffix="dseg.nii.gz",
        ),
        mask_res=rules.resample_brainmask_tractmaps.output.mask_res,
    output:
        seed_res=bids(
            root="results/tractmap",
            subject="{subject}",
            space="individual",
            label="{seed}",
            method="spectralcosine",
            k="{k}",
            from_="{template}",
            res="super",
            desc="sorted",
            suffix="dseg.nii.gz",
        ),
    container:
        config["singularity"]["neuroglia"]
    log:
        "logs/resample_clus_seed/sub-{subject}_template-{template}_{seed}_k-{k}.log",
    group:
        "participant2"
    shell:
        #linear interp here now, since probabilistic seg
        "reg_resample -flo {input.seed} -res {output.seed_res} -ref {input.mask_res} -NN 0 &> {log}"


"""

#split segmentation into binary masks
# space-T1w   mask
rule subj_split_clus_to_binary_masks:
    input: 
        cluster_k = bids(root='results/tractmap',subject='{subject}',space='individual',label='{seed}',method='spectralcosine',k='{k}',from_='{template}',res='super',desc='sorted',suffix='dseg.nii.gz'),
    params:
        mask_file = lambda wildcards, output: bids(root=output.cluster_k_splitdir,subject=wildcards.subject,label='%02d',suffix='mask.nii.gz',include_subject_dir=False),
        mask_bg = lambda wildcards, output: bids(root=output.cluster_k_splitdir,subject=wildcards.subject,label='00',suffix='mask.nii.gz',include_subject_dir=False) #we remove this file 
    output:
        seed_label = bids(root='results/tractmap',subject='{subject}',space='individual',label='{seed}',method='spectralcosine',k='{k}',from_='{template}',res='super',desc='sorted',kindex='{kindex}',suffix='seed.nii.gz')
    container: config['singularity']['neuroglia']
    log: 'logs/subj_split_clus_to_binary_masks/sub-{subject}_template-{template}_{seed}_k-{k}.log'
    group: 'participant2'
    threads: 8
    shell:
        #use c3d's split command to go from discrete seg image to multiple binary images.. we remove image 00 since that is the background label
        'mkdir {output.cluster_k_splitdir} && c3d {input.cluster_k} -split -oo {params.mask_file} &>> {log}  && rm -f {params.mask_bg}'


"""


# perform tracking from each cluster in subj space to get tract maps
rule track_from_clusters:
    input:
        cluster_k=bids(
            root="results/tractmap",
            subject="{subject}",
            space="individual",
            label="{seed}",
            method="spectralcosine",
            k="{k}",
            from_="{template}",
            res="super",
            desc="sorted",
            suffix="dseg.nii.gz",
        ),
        mask=rules.resample_brainmask_tractmaps.output.mask,
    params:
        bedpost_merged=join(
            config["bedpost"]["dir"], config["bedpost"]["merged_prefix"]
        ),
        probtrack_opts=config["probtrack_tractmap"]["opts"],
        nsamples=config["probtrack_tractmap"]["nsamples"],
        extract_seed_cmd=lambda wildcards, input, output: f"fslmaths {input.cluster_k} -thr {wildcards.kindex} -uthr {wildcards.kindex} {output.probtrack_dir}/in_seed.nii.gz",
        out_track_dirs=lambda wildcards, output: expand(
            "{}/label-{{k_index:02d}}".format(output.probtrack_dir),
            k_index=range(1, int(wildcards.k) + 1),
        ),
        container=config["singularity"]["fsl_cuda"],
    output:
        probtrack_dir=directory(
            bids(
                root="results/tractmap",
                subject="{subject}",
                space="individual",
                label="{seed}",
                method="spectralcosine",
                k="{k}",
                from_="{template}",
                res="super",
                suffix="probtrack",
                kindex="{kindex}",
            )
        ),
        tractmap=bids(
            root="results/tractmap",
            subject="{subject}",
            space="individual",
            label="{seed}",
            method="spectralcosine",
            k="{k}",
            from_="{template}",
            res="super",
            suffix="probtrack/fdt_paths.nii.gz",
            kindex="{kindex}",
        ),
    threads: 1
    resources:
        mem_mb=64000,  #set to 64000 to ensure only 2 instances can run per node 
        time=360,
        gpus=1,  #1 gpu
    log:
        "logs/track_from_clusters/sub-{subject}_template-{template}_{seed}_k-{k}_kindex-{kindex}.log",
    group:
        "participant2"
    shell:
        "mkdir -p {output.probtrack_dir} && "
        "{params.extract_seed_cmd} && singularity exec -e --nv {params.container} "
        "probtrackx2_gpu --samples={params.bedpost_merged}  --mask={input.mask} --seed={output.probtrack_dir}/in_seed.nii.gz "
        "--seedref={output.probtrack_dir}/in_seed.nii.gz --nsamples={params.nsamples} "
        "--dir={output.probtrack_dir} {params.probtrack_opts} -V 2  &> {log}"


# check bids-deriv dwi - tractography ?
# space-T1w, res-?
rule combine_tractmaps:
    input:
        tractmaps=lambda wildcards: expand(
            bids(
                root="results/tractmap",
                subject=f"{wildcards.subject}",
                space="individual",
                label=f"{wildcards.seed}",
                method="spectralcosine",
                k=f"{wildcards.k}",
                from_=f"{wildcards.template}",
                res="super",
                suffix="probtrack/fdt_paths.nii.gz",
                kindex="{kindex}",
            ),
            kindex=range(1, int(wildcards.k) + 1),
        ),
    output:
        tractmaps_4d=bids(
            root="results/tractmap",
            subject="{subject}",
            space="individual",
            label="{seed}",
            method="spectralcosine",
            k="{k}",
            from_="{template}",
            res="super",
            suffix="tractmap4d.nii.gz",
        ),
    log:
        "logs/combine_tractmaps/sub-{subject}_template-{template}_{seed}_k-{k}.log",
    container:
        config["singularity"]["neuroglia"]
    group:
        "participant2"
    resources:
        mem_mb=32000,
    shell:
        "fslmerge -t {output.tractmaps_4d} {input.tractmaps} &> {log}"


# transform tract maps back to template space
# space-{template}, tractogrpahy?
rule transform_tractmaps_to_template:
    input:
        tractmap=bids(
            root="results/tractmap",
            subject="{subject}",
            space="individual",
            label="{seed}",
            method="spectralcosine",
            k="{k}",
            from_="{template}",
            res="super",
            suffix="probtrack/fdt_paths.nii.gz",
            kindex="{kindex}",
        ),
        warp=config["transforms"]["warp"],
        ref=config["transforms"]["ref_nii"],
    output:
        tractmap=bids(
            root="results/tractmap",
            subject="{subject}",
            label="{seed}",
            method="spectralcosine",
            k="{k}",
            from_="{template}",
            res="super",
            space="{template}",
            suffix="tractmap.nii.gz",
            kindex="{kindex}",
        ),
    container:
        config["singularity"]["neuroglia"]
    log:
        "logs/transform_tractmaps_to_template/sub-{subject}_{seed}_{template}_k-{k}_kindex-{kindex}.log",
    group:
        "participant2"
    threads: 8
    resources:
        mem_mb=32000,
        time=10,
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1  antsApplyTransforms -d 3 --interpolation Linear -i {input.tractmap} -o {output.tractmap}  -r {input.ref} -t {input.warp}  &> {log}"


# space-{template}, tractography 4d?
rule combine_tractmaps_warped:
    input:
        tractmaps=lambda wildcards: expand(
            bids(
                root="results/tractmap",
                subject=f"{wildcards.subject}",
                label=f"{wildcards.seed}",
                method="spectralcosine",
                k=f"{wildcards.k}",
                from_=f"{wildcards.template}",
                res="super",
                space=f"{wildcards.template}",
                suffix="tractmap.nii.gz",
                kindex="{kindex}",
            ),
            kindex=range(1, int(wildcards.k) + 1),
        ),
    output:
        tractmaps_4d=bids(
            root="results/tractmap",
            subject="{subject}",
            label="{seed}",
            method="spectralcosine",
            k="{k}",
            space="{template}",
            suffix="tractmap4d.nii.gz",
        ),
    log:
        "logs/combine_tractmaps_warped/sub-{subject}_template-{template}_{seed}_k-{k}.log",
    container:
        config["singularity"]["neuroglia"]
    resources:
        mem_mb=32000,
    group:
        "participant2"
    shell:
        "fslmerge -t {output.tractmaps_4d} {input.tractmaps} &> {log}"


# now, average the tract maps (over subjects) in template space
# space-template, desc-avg , tractography 4d?
rule avg_tractmaps_template:
    input:
        tractmaps_4d=expand(
            rules.combine_tractmaps_warped.output.tractmaps_4d,
            subject=subjects,
            allow_missing=True,
        ),
    output:
        average=bids(
            root="results/tractmap",
            template="{template}",
            label="{seed}",
            method="spectralcosine",
            k="{k}",
            desc="average",
            suffix="tractmap4d.nii.gz",
        ),
    container:
        config["singularity"]["neuroglia"]
    group:
        "group2"
    threads: 8
    resources:
        mem_mb=32000,
    shell:
        "AverageImages 4 {output} 0 {input}"


# use voting to get a discrete segmentation of tract map
# space-template, desc-avgtractmap, dseg
rule vote_tractmap_template:
    input:
        tractmaps=rules.avg_tractmaps_template.output.average,
    params:
        bg_th=100,  # set only if avg streamline count > bg_th 
    output:
        vote_seg=bids(
            root="results/tractmap",
            template="{template}",
            label="{seed}",
            method="spectralcosine",
            k="{k}",
            desc="avgtractmap",
            suffix="dseg.nii.gz",
        ),
    container:
        config["singularity"]["neuroglia"]
    log:
        "logs/vote_tractmap_template/{template}_{seed}_{k}.log",
    group:
        "group2"
    threads: 8
    resources:
        mem_mb=32000,
    shell:
        #get first vol, set it to bg_th value, this becomes image 0
        # then load all the tractmaps as images 1-k
        # then vote amongst the 0-k images, those where bg (bg_th) is highest get set to 0, otherwise set to the index (from 1-k)
        "c4d -verbose {input.tractmaps} -slice w 0:0 -threshold -inf +inf {params.bg_th} 0 {input.tractmaps} -slice w 0:-1  -vote -type uchar {output.vote_seg} &> {log}"
