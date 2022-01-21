"""
This simple snakefile moves and renames the images Andreas wants to quality check.
"""


rule move_ssvep_amplitudes:
    """
    Standardize our many poorly named filenames into a common template.
    """
    input:
        head="../processed/correlation_whole_brain_ttest/startvolume-{startvolume}_variable-amplitudes_improved_correlations_ttest+tlrc.HEAD",
        brik="../processed/correlation_whole_brain_ttest/startvolume-{startvolume}_variable-amplitudes_improved_correlations_ttest+tlrc.HEAD"
    output: "results/images_for_andreas/ssvep_amplitude_lag-{startvolume}_correlation.nii.gz"
    conda: "../envs/afni.yaml"
    shell: "3dTcat '{input.head}[0]' -prefix '{output}'"

