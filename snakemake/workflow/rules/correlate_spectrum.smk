"""
Workflow to gather quadratic parameters of the spectrums in each TR and correlate those with the BOLD signal.
"""


rule get_quadratic_estimates_of_spectrums:
    """
    Get time series of quadratic estimate coefficients for each TR.
    """
    input:
        pows_timeseries="../processed/eeg_alphas/sub-{id}_pows_timeseries.mat",
        freqs="../processed/eeg_alphas/sub-{id}_freqs.mat",
    output:
        zero_order="results/spectrum_quadratics/estimates/sub-{id}_order-0_estimate.csv",
        first_order="results/spectrum_quadratics/estimates/sub-{id}_order-1_estimate.csv",
        second_order="results/spectrum_quadratics/estimates/sub-{id}_order-2_estimate.csv",
    conda: "../envs/neuroimaging.yaml"
    script: "../scripts/get_quadratic_estimates.py"


rule correlate_quadratic_estimates_with_BOLD:
    """
    Trim BOLD to be same length as quadratic estimates then correlate estimates with BOLD.
    Input BOLD series and quadratic estimates should begin at about the same time.
    """
    input:
        vector="results/spectrum_quadratics/estimates/sub-{id}_order-{order}_estimate.csv",
        bold="../processed/funcs_trimmed_to_first_stim/sub-{id}_func_trimmed+tlrc.HEAD",
    output:
        correlations="results/spectrum_quadratics/correlations/sub-{id}_order-{order}_correlations.nii.gz"
    conda: "../envs/neuroimaging.yaml"
    script: "../scripts/correlate_with_bold.py"


rule correlate_quadratic_estimates_with_residuals:
    """
    Trim residuals to be same length as quadratic estimates then correlate estimates with BOLD.
    Input BOLD series and quadratic estimates should begin at about the same time.
    """
    input:
        vector="results/spectrum_quadratics/estimates/sub-{id}_order-{order}_estimate.csv",
        bold="results/prepare_residuals/lag_residuals/sub-{id}_lag-{lag}.nii.gz",
    output:
        correlations="results/spectrum_quadratics/correlation_with_residuals/sub-{id}_lag-{lag}_order-{order}.nii.gz"
    conda: "../envs/neuroimaging.yaml"
    script: "../scripts/correlate_with_bold.py"
