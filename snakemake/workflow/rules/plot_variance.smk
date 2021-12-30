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
        threshold=lambda wildcards: (-round(max([threshold ** 2 for threshold in config['thresholds'][wildcards.analysis][wildcards.percentile]]), 4), round(max([threshold ** 2 for threshold in config['thresholds'][wildcards.analysis][wildcards.percentile]]), 4)),
        coordinates=config["occipital coordinates"],
        title=lambda wildcards: f"{wildcards.analysis} {wildcards.variable} mean variance\ntime lag: {wildcards.startvolume} volumes\nbaseline subtracted: {wildcards.baselined}\ncritical values={(-round(max([threshold ** 2 for threshold in config['thresholds'][wildcards.analysis][wildcards.percentile]]), 4), round(max([threshold ** 2 for threshold in config['thresholds'][wildcards.analysis][wildcards.percentile]]), 4))}\npercentile: {wildcards.percentile}",
    conda:
        "../envs/neuroplotting.yaml"
    script:
        "../scripts/plot_fmri.py"


rule stitch_thresholded_variance_plots:
    """
    Stitch thresholded plots together in nice, easy to read order.
    """
    input:
        plots=lambda wildcards: [f"results/plot_variance/permutation_thresholded_plots/startvolume-{startvolume}_variable-{wildcards.variable}_baselined-{wildcards.baselined}_percentile-{wildcards.percentile}_{wildcards.analysis}_occipital.png" for startvolume in config["analyses"][wildcards.analysis]["continuous"]["volumes"]]
    output:
        stitched_plot="results/plot_variance/stitched_thresholded_plots/variable-{variable}_baselined-{baselined}_percentile-{percentile}_{analysis}_occipital.png"
    conda:
        "../envs/imagemagick.yaml"
    shell:
        "convert {input.plots} -append {output.stitched_plot}"
