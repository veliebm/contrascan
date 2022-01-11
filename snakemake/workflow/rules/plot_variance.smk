"""
Plot mean variance images.
"""

rule strip_skull_from_variance:
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
    Plot the occipital pole of the thresholded variances using EXACT variance calculations.
    """
    input:
        underlay="../data/misc/kastner_cortex_masks/MNI152_T1_1mm_masked.nii.gz",
        overlay="results/plot_variance/skull_stripped/startvolume-{startvolume}_variable-{variable}_baselined-{baselined}_{analysis}.nii.gz",
        mins_distribution="../processed/maxes_and_mins_distributions/startvolume-{startvolume}_variable-{variable}_baselined-{baselined}_{analysis}_mins_table.csv",
        maxes_distribution="../processed/maxes_and_mins_distributions/startvolume-{startvolume}_variable-{variable}_baselined-{baselined}_{analysis}_maxes_table.csv",
    output:
        plot="results/plot_variance/exact_permutation_thresholded_plots/startvolume-{startvolume}_variable-{variable}_baselined-{baselined}_percentile-{percentile}_{analysis}_occipital.png",
    params:
        coordinates=config["occipital coordinates"],
    conda:
        "../envs/neuroplotting.yaml"
    script:
        "../scripts/plot_permutation_thresholded_variance.py"


rule plot_rougher_occipital_thresholded_variance:
    """
    Plot variance using the estimates I found.
    """
    input:
        underlay="../data/misc/kastner_cortex_masks/MNI152_T1_1mm_masked.nii.gz",
        overlay="results/plot_variance/skull_stripped/startvolume-{startvolume}_variable-{variable}_baselined-{baselined}_{analysis}.nii.gz",
    output:
        plot="results/plot_variance/permutation_thresholded_plots/startvolume-{startvolume}_variable-{variable}_baselined-{baselined}_percentile-{percentile}_{analysis}_occipital.png",
    params:
        threshold=lambda wildcards: config['thresholds']['variance'][f"baselined {wildcards.baselined}"],
        coordinates=config['occipital coordinates'],
        title=lambda wildcards: f"{wildcards.analysis} {wildcards.variable} mean correlation\ntime lag: {wildcards.startvolume} volumes\nbaseline subtracted: {wildcards.baselined}\ncritical values={config['thresholds']['variance'][f'baselined {wildcards.baselined}']}\npercentile: {wildcards.percentile}",
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
