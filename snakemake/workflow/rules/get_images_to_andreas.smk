"""
This simple snakefile moves and renames the images Andreas wants to quality check.
"""


rule move_ssvep_amplitudes:
    """
    Standardize our many poorly named filenames into a common template.
    """
    input:
        head="../processed/correlation_whole_brain_ttest/startvolume-{startvolume}_variable-amplitudes_improved_correlations_ttest+tlrc.HEAD",
        brik="../processed/correlation_whole_brain_ttest/startvolume-{startvolume}_variable-amplitudes_improved_correlations_ttest+tlrc.BRIK"
    output: "results/images_for_andreas/ssvep_amplitude_correlation_lag-{startvolume}.nii.gz"
    conda: "../envs/afni.yaml"
    shell: "3dTcat '{input.head}[0]' -prefix '{output}'"


rule move_ssvep_SNRs:
    """
    Standardize our many poorly named filenames into a common template.
    """
    input:
        head="../processed/correlation_whole_brain_ttest/startvolume-{startvolume}_variable-SNRs_improved_correlations_ttest+tlrc.HEAD",
        brik="../processed/correlation_whole_brain_ttest/startvolume-{startvolume}_variable-SNRs_improved_correlations_ttest+tlrc.BRIK"
    output: "results/images_for_andreas/ssvep_SNR_correlation_lag-{startvolume}.nii.gz"
    conda: "../envs/afni.yaml"
    shell: "3dTcat '{input.head}[0]' -prefix '{output}'"


rule move_alpha_amplitudes:
    """
    Standardize our many poorly named filenames into a common template.
    """
    input:
        head="../processed/correlation_whole_brain_ttest/startvolume-{startvolume}_data-values_correlations_alpha_ttest+tlrc.HEAD",
        brik="../processed/correlation_whole_brain_ttest/startvolume-{startvolume}_data-values_correlations_alpha_ttest+tlrc.BRIK"
    output: "results/images_for_andreas/alpha_amplitude_correlation_lag-{startvolume}.nii.gz"
    conda: "../envs/afni.yaml"
    shell: "3dTcat '{input.head}[0]' -prefix '{output}'"


rule move_alpha_SNRs:
    """
    Standardize our many poorly named filenames into a common template.
    """
    input:
        head="../processed/correlation_whole_brain_ttest/startvolume-{startvolume}_data-SNRs_correlations_alpha_ttest+tlrc.HEAD",
        brik="../processed/correlation_whole_brain_ttest/startvolume-{startvolume}_data-SNRs_correlations_alpha_ttest+tlrc.BRIK"
    output: "results/images_for_andreas/alpha_SNR_correlation_lag-{startvolume}.nii.gz"
    conda: "../envs/afni.yaml"
    shell: "3dTcat '{input.head}[0]' -prefix '{output}'"
