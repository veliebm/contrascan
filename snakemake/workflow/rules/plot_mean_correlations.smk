from pathlib import Path
from os import PathLike
from typing import List


def get_poorly_named_input(wildcards) -> Path:
    """
    Given some wildcards, return the corresponding poorly named input file to standardize.
    """
    startvolume = wildcards.startvolume
    variable = wildcards.variable
    analysis = wildcards.analysis

    if analysis == "alpha" and startvolume.isdigit() and variable == "amplitude":
        filename = f"startvolume-{startvolume}_data-values_correlations_alpha_ttest+tlrc.HEAD"
    elif analysis == "alpha" and startvolume.isdigit() and variable == "SNR":
        filename = f"startvolume-{startvolume}_data-SNRs_correlations_alpha_ttest+tlrc.HEAD"
    elif analysis == "alpha" and startvolume == "na" and variable == "amplitude":
        filename = "variable-alpha_trials_correlations_ttest+tlrc.HEAD"
    elif analysis == "ssvep" and startvolume.isdigit() and variable == "amplitude":
        filename = f"startvolume-{startvolume}_variable-amplitudes_improved_correlations_ttest+tlrc.HEAD"
    elif analysis == "ssvep" and startvolume.isdigit() and variable == "SNR":
        filename = f"startvolume-{startvolume}_variable-SNRs_improved_correlations_ttest+tlrc.HEAD"
    elif analysis == "ssvep" and startvolume == "na" and variable == "amplitude":
        filename = f"variable-slidewinamp_trials_correlations_ttest+tlrc.HEAD"
    elif analysis == "ssvep" and startvolume == "na" and variable == "SNR":
        filename = f"variable-slidewinSNR_trials_correlations_ttest+tlrc.HEAD"

    try:
        return Path("../processed/correlation_whole_brain_ttest") / filename
    except NameError:
        raise NameError(f"No matching filename found: {wildcards}")


rule standardize_filenames:
    """
    Standardize our many poorly named filenames into a common template.
    """
    input:
        get_poorly_named_input
    output:
        "results/plot_mean_correlations/standardize_filenames/startvolume-{startvolume}_variable-{variable}_{analysis}_correlations_ttest.nii.gz"
    log:
        "logs/plot_mean_correlations/standardize_filenames/startvolume-{startvolume}_variable-{variable}_{analysis}_correlations_ttest.log"
    conda:
        "../envs/afni.yaml"
    shell:
        "3dcopy '{input}' '{output}' 2> {log}"


rule extract_mean:
    """
    Extract just the means from our 3dttest++ results.
    """
    input:
        image="results/plot_mean_correlations/standardize_filenames/startvolume-{startvolume}_variable-{variable}_{analysis}_correlations_ttest.nii.gz",
    output:
        image="results/plot_mean_correlations/mean_only/startvolume-{startvolume}_variable-{variable}_baselined-false_{analysis}.nii.gz",
    log:
        "logs/plot_mean_correlations/mean_only/startvolume-{startvolume}_variable-{variable}_baselined-false_{analysis}.log"
    conda:
        "../envs/afni.yaml"
    shell:
        "3dTcat '{input}[0]' -prefix '{output}' 2> {log}"


rule copy_baseline_corrected_means:
    """
    Copy our baseline corrected means to be with the rest of our files.
    """
    input:
        "../processed/compare_to_canonical/startvolume-{startvolume}_variable-{variable}_{analysis}_baselinecorrected.nii"
    output:
        "results/plot_mean_correlations/mean_only/startvolume-{startvolume}_variable-{variable}_baselined-true_{analysis}.nii.gz"
    log:
        "logs/plot_mean_correlations/mean_only/startvolume-{startvolume}_variable-{variable}_baselined-true_{analysis}.log"
    conda:
        "../envs/afni.yaml"
    shell:
        "3dcopy '{input}' '{output}' 2> {log}"


rule strip_skull:
    """
    Mask out the skull from our means.
    """
    input:
        mask="../data/misc/kastner_cortex_masks/MNI152_T1_2.5mm_full_mask.nii.gz",
        image="results/plot_mean_correlations/mean_only/startvolume-{startvolume}_variable-{variable}_baselined-{baselined}_{analysis}.nii.gz"
    output:
        image="results/plot_mean_correlations/skull_stripped/startvolume-{startvolume}_variable-{variable}_baselined-{baselined}_{analysis}.nii.gz",
    log:
        "logs/plot_mean_correlations/skull_stripped/startvolume-{startvolume}_variable-{variable}_baselined-{baselined}_{analysis}.log"
    conda:
        "../envs/afni.yaml"
    shell:
        "3dcalc -float -a {input.image} -b {input.mask} -expr 'a*step(b)' -prefix {output.image} 2> {log}"


rule plot_occipital:
    """
    Plot the occipital pole for each of our means.
    """
    input:
        underlay="../data/misc/kastner_cortex_masks/MNI152_T1_1mm_masked.nii.gz",
        overlay="results/plot_mean_correlations/skull_stripped/startvolume-{startvolume}_variable-{variable}_baselined-{baselined}_{analysis}.nii.gz",
    output:
        plot=report(
            "results/plot_mean_correlations/plots/startvolume-{startvolume}_variable-{variable}_baselined-{baselined}_{analysis}_occipital.png",
            caption="../report/mean_correlations.rst",
            category="startvolume {startvolume}",
            subcategory="{analysis} {variable}",
        ),
    params:
        threshold=(-.05, .05),
        coordinates=config["occipital coordinates"],
        title="{analysis} {variable}, {startvolume} volumes removed, baselined={baselined}, threshold={params.threshold}"
    conda:
        "../envs/neuroplotting.yaml"
    script:
        "../scripts/plot_fmri.py"


rule stitch_together_plots:
    """
    Stitch plots together in nice, easy to read order.
    """
    input:
        plots=lambda wildcards: [f"results/plot_mean_correlations/plots/startvolume-{startvolume}_variable-{wildcards.variable}_baselined-{wildcards.baselined}_{wildcards.analysis}_occipital.png" for startvolume in config["analyses"][wildcards.analysis]["continuous"]["volumes"]]
    output:
        stitched_plot="results/plot_mean_correlations/stitched_plots/variable-{variable}_baselined-{baselined}_{analysis}_occipital.png"
    conda:
        "../envs/imagemagick.yaml"
    shell:
        "convert {input.plots} -append {output.stitched_plot}"
