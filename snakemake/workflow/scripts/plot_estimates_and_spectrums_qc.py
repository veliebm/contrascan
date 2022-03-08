"""
Plot quadratic curve and spectrum for a TR.

Created 3/7/22 by Benjamin Velie.
"""
import scipy.io
import numpy
import matplotlib.pyplot as plt


def main():
    """
    Entry point of the script.
    """
    # Read the spectrum and its freqs for the given TR.
    TR = int(snakemake.wildcards.TR)
    OZ_ELECTRODE = 19
    pows_timeseries = scipy.io.loadmat(snakemake.input.pows_timeseries)["pows"][OZ_ELECTRODE]
    TR_spectrum = pows_timeseries[:, TR]
    TR_freqs = scipy.io.loadmat(snakemake.input.freqs)["freqs"].reshape(-1)

    # Trime freqs and spectrum so we only see first 30 Hz.
    TR_trimmed_spectrum = TR_spectrum[:61]
    TR_trmmed_freqs = TR_freqs[:61]

    # Read quadratic coefficients for the given TR.
    TR_quadratic_estimates = [
        numpy.genfromtxt(snakemake.input.zero_order)[TR],
        numpy.genfromtxt(snakemake.input.first_order)[TR],
        numpy.genfromtxt(snakemake.input.second_order)[TR],
    ]

    # Reconstruct quadratic model.
    model = numpy.polynomial.Polynomial(TR_quadratic_estimates)
    TR_quadratic_fit = model(TR_trmmed_freqs)

    # Using the coefficients, plot a regression over the original data.
    fig, ax = plt.subplots()
    ax.plot(TR_trmmed_freqs, TR_trimmed_spectrum, 's', color='#377eb8', marker='o')
    ax.plot(TR_trmmed_freqs, TR_quadratic_fit, color='#ff7f00')
    ax.set_ylabel('Power')
    ax.set_xlabel('Frequency (Hz)')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.set_size_inches(10, 9)
    plt.savefig(snakemake.output.plot)


if __name__ == "__main__":
    main()
