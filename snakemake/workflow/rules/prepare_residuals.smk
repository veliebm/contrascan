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
    conda: "../envs/afni.yaml"
    shell: "3dAllineate -cubic -1Dmatrix_apply {input.alignment_matrix} -prefix {output.image} {input.image}"
