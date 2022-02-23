"""
Workflow to gather quadratic parameters of the spectrums in each TR and correlate those with the BOLD signal.

Created 2/23/22 by Ben Velie.
"""


rule get_quadratic_parameters_of_spectrums:
    """
    Get time series of quadratic parameters for each TR.
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

