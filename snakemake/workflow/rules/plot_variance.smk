"""
Plot mean variance images.
"""

rule strip_skull_of_variance:
    """
    Mask out the skull from our variance images.
    """
    input:
        mask="../data/misc/kastner_cortex_masks/MNI152_T1_2.5mm_full_mask.nii.gz",
        image="../processed/variance_whole_brain/startvolume-{startvolume}_variable-{variable}_baselined-{baselined}_{analysis}.nii.gz"
    output:
        image="results/plot_variance/skull_stripped/startvolume-{startvolume}_variable-{variable}_baselined-{baselined}_{analysis}.nii.gz",
    log:
        "logs/plot_variance/skull_stripped/startvolume-{startvolume}_variable-{variable}_baselined-{baselined}_{analysis}.log"
    conda:
        "../envs/afni.yaml"
    shell:
        "3dcalc -float -a {input.image} -b {input.mask} -expr 'a*step(b)' -prefix {output.image} 2> {log}"


rule plot_occipital_thresholded_variance:
    """
    Plot the occipital pole of the thresholded variances.
    """
    input:
        underlay="../data/misc/kastner_cortex_masks/MNI152_T1_1mm_masked.nii.gz",
        overlay="results/plot_variance/skull_stripped/startvolume-{startvolume}_variable-{variable}_baselined-{baselined}_{analysis}.nii.gz",
    output:
        plot="results/plot_variance/permutation_thresholded_plots/startvolume-{startvolume}_variable-{variable}_baselined-{baselined}_percentile-{percentile}_{analysis}_occipital.png",
    params:
        threshold=lambda wildcards: [number ** 2 for number in config["thresholds"][wildcards.analysis][wildcards.percentile]],
        coordinates=config["occipital coordinates"],
        title=lambda wildcards: f"{wildcards.analysis} {wildcards.variable} mean correlation\ntime lag: {wildcards.startvolume} volumes\nbaseline subtracted: {wildcards.baselined}\ncritical values={config['thresholds'][wildcards.analysis][wildcards.percentile]}\npercentile: {wildcards.percentile}",
    conda:
        "../envs/neuroplotting.yaml"
    script:
        "../scripts/plot_fmri.py"