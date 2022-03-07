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

    # Read quadratic coefficients for the given TR.
    TR_quadratic_estimates = [
        numpy.genfromtxt(snakemake.input.zero_order)[TR],
        numpy.genfromtxt(snakemake.input.first_order)[TR],
        numpy.genfromtxt(snakemake.input.second_order)[TR],
    ]

    # Reconstruct quadratic model.

    # Using the coefficients, plot a regression over the original data.
    # Should work if I use the same input data that was used to create the original data.
    # i.e., f(x) and f_est(x)
    # What's the original data?
    # When I made the polynomial, I used the freqs as the X input, and the spectrum as the Y input.
    model = numpy.polynomial.Polynomial(TR_quadratic_estimates)
    TR_quadratic_fit = model(TR_freqs)

    # Make the plot.
    fig, ax = plt.subplots()
    ax.plot(TR_freqs, TR_spectrum, 's', color='#377eb8', marker='o')
    ax.plot(TR_freqs, TR_quadratic_fit, color='#ff7f00')
    ax.set_ylabel('Power')
    ax.set_xlabel('Frequency (Hz)')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.set_size_inches(10, 9)
    plt.savefig(snakemake.output.plot)


if __name__ == "__main__":
    main()
