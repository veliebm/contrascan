"""
Workflow to remove motion regressors from BOLD signal without modeling the HRF at the same time.
"""

rule align_residuals_made_by_afniproc:
    """
    Align the residuals made by afni_proc.py to the Kastner space.
    """
    input:
        image="../processed/afniproc/sub-{id}/{id}.results/errts.{id}+tlrc.HEAD",
        alignment_matrix="../processed/IRF_alignment/MNI152_T1_2009c_aligned_mat.aff12.1D",
    output:
        image="results/prepare_residuals/align_residuals/sub-{id}.nii.gz",
    shell: "3dAllineate -cubic -1Dmatrix_apply {input.alignment_matrix} -prefix {output.image} {input.image}"


rule resample_aligned_residuals:
    """
    Resample our aligned residuals.
    """
    input:
        image="results/prepare_residuals/align_residuals/sub-{id}.nii.gz",
        template="../processed/kastner_resample/MNI152_T1_1mm_resampled+tlrc.HEAD",
    output:
        image="results/prepare_residuals/resample_residuals/sub-{id}.nii.gz",
    shell: "3dresample -master {input.template} -input {input.image} -prefix {output.image}"


rule trim_resampled_residuals:
    """
    Truncate our resampled residuals to begin such that the first stimulus is presented somewhere within the first volume.

    Hardcode 4 volumes removed because it's the same for all images.
    """
    input:
        image="results/prepare_residuals/resample_residuals/sub-{id}.nii.gz",
    output:
        image="results/prepare_residuals/trim_residuals/sub-{id}.nii.gz",
    shell: "3dTcat {input.image}[4..$] -prefix {output.image}"


rule add_lags_to_trimmed_residuals:
    """
    Lag our trimmed residuals. Andreas says these should be lagged by 3 to 5 TRs.
    """
    input:
        image="results/prepare_residuals/trim_residuals/sub-{id}.nii.gz",
    output:
        image="results/prepare_residuals/lag_residuals/sub-{id}_lag-{lag}.nii.gz",
    shell: "3dTcat {input.image}[{wildcards.lag}..$] -prefix {output.image}"
