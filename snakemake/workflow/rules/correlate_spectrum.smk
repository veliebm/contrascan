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
        zero_order="results/correlate_spectrum/estimates/sub-{id}_order-0_estimate.csv",
        first_order="results/correlate_spectrum/estimates/sub-{id}_order-1_estimate.csv",
        second_order="results/correlate_spectrum/estimates/sub-{id}_order-2_estimate.csv",
    conda: "../envs/neuroimaging.yaml"
    script: "../scripts/get_quadratic_estimates.py"


rule plot_estimates_and_spectrums:
    """
    For a TR, plot a spectrum and its fitted quadratic curve.
    """
    input:
        zero_order="results/correlate_spectrum/estimates/sub-{id}_order-0_estimate.csv",
        first_order="results/correlate_spectrum/estimates/sub-{id}_order-1_estimate.csv",
        second_order="results/correlate_spectrum/estimates/sub-{id}_order-2_estimate.csv",
        pows_timeseries="../processed/eeg_alphas/sub-{id}_pows_timeseries.mat",
        freqs="../processed/eeg_alphas/sub-{id}_freqs.mat",
    output:
        plot="results/correlate_spectrum/estimates/sub-{id}_TR-{TR}_plot.png",
    conda: "../envs/neuroimaging.yaml"
    script: "../scripts/plot_estimates_and_spectrums_qc.py"


rule correlate_quadratic_estimates_with_residuals:
    """
    Trim residuals to be same length as quadratic estimates then correlate estimates with BOLD.
    Input BOLD series and quadratic estimates should begin at about the same time.
    """
    input:
        vector="results/correlate_spectrum/estimates/sub-{id}_order-{order}_estimate.csv",
        bold="results/prepare_residuals/lag_residuals/sub-{id}_lag-{lag}.nii.gz",
    output:
        correlations="results/correlate_spectrum/correlation_with_residuals/sub-{id}_lag-{lag}_order-{order}.nii.gz"
    wildcard_constraints:
        lag="[3-5]"
    conda: "../envs/neuroimaging.yaml"
    script: "../scripts/correlate_with_bold.py"


rule compute_mean_correlation:
    """
    Compute the mean correlation across subjects.
    """
    input:
        images=expand("results/correlate_spectrum/correlation_with_residuals/sub-{id}_lag-{{lag}}_order-{{order}}.nii.gz", id=config["subjects"])
    output:
        mean="results/correlate_spectrum/correlation_with_residuals_mean/lag-{lag}_order-{order}.nii.gz"
    shell: "3dMean -prefix {output.mean} {input.images}"
