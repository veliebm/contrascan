from pathlib import Path

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
        "results/plot_mean_correlations/standardize_filenames/startvolume-{startvolume}_variable-{variable}_{analysis}_correlations_ttest.log"
    shell:
        "3dcopy '{input}' '{output}' 2> {log}"


rule extract_mean:
    """
    Extract just the means from our 3dttest++ results.
    """
    input:
        image="results/plot_mean_correlations/standardize_filenames/startvolume-{startvolume}_variable-{variable}_{analysis}_correlations_ttest.nii.gz",
    output:
        image="results/plot_mean_correlations/mean_only/startvolume-{startvolume}_variable-{variable}_{analysis}.nii.gz",
    conda:
        "../envs/neuroimaging.yaml"
    log:
        "logs/plot_mean_correlations/mean_only/startvolume-{startvolume}_variable-{variable}_{analysis}.log"
    shell:
        "3dTcat '{input}[0]' -prefix '{output}' 2> {log}"


rule strip_skull:
    """
    Mask out the skull from our means.
    """
    input:
        mask="../data/misc/kastner_cortex_masks/MNI152_T1_2.5mm_full_mask.nii.gz",
        image="results/plot_mean_correlations/mean_only/startvolume-{startvolume}_variable-{variable}_{analysis}.nii.gz"
    output:
        image="results/plot_mean_correlations/skull_stripped/startvolume-{startvolume}_variable-{variable}_{analysis}.nii.gz",
    conda:
        "../envs/neuroimaging.yaml"
    log:
        "logs/plot_mean_correlations/skull_stripped/startvolume-{startvolume}_variable-{variable}_{analysis}.log"
    shell:
        "3dcalc -float -a {input.image} -b {input.mask} -expr 'a*step(b)' -prefix {output.image} 2> {log}"


rule plot_occipital:
    """
    Plot the occipital pole for each of our means.
    """
    input:
        underlay="../data/misc/kastner_cortex_masks/MNI152_T1_1mm_masked.nii.gz",
        overlay="results/plot_mean_correlations/skull_stripped/startvolume-{startvolume}_variable-{variable}_{analysis}.nii.gz",
    output:
        plot="results/plot_mean_correlations/plots/startvolume-{startvolume}_variable-{variable}_{analysis}_occipital.png"
    params:
        threshold=0.1,
        coordinates=config["occipital coordinates"],
    conda:
        "../envs/neuroimaging.yaml"
    script:
        "../scripts/plot_fmri.py"
 