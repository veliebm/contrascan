"""
Get time series of quadratic parameters for each TR.

Created 2/23/22 by Benjamin Velie.
"""
import scipy.io
import numpy


def main():
    """
    Entry point of the script.
    """
    # Read time series and freqs from .mat files.
    freqs = scipy.io.loadmat(snakemake.input.freqs)["freqs"].reshape(-1)
    oz_electrode = 19
    pows_timeseries = scipy.io.loadmat(snakemake.input.pows_timeseries)["pows"][oz_electrode]

    # Trim freqs and time series to only get first 30 Hz.
    trimmed_freqs = freqs[:61]
    trimmed_pows_timeseries = pows_timeseries[:61, :]

    # Fit quadratic model for each TR. Example pows_timeseries shape: (500, 345)
    def fit_quadratic(spectrum: numpy.array) -> numpy.array:
        return numpy.polynomial.polynomial.polyfit(x=trimmed_freqs, y=spectrum, deg=2)
    coefficients = numpy.apply_along_axis(fit_quadratic, 0, trimmed_pows_timeseries)

    # Output quadratic parameters.
    numpy.savetxt(snakemake.output.zero_order, coefficients[0])
    numpy.savetxt(snakemake.output.first_order, coefficients[1])
    numpy.savetxt(snakemake.output.second_order, coefficients[2])


if __name__ == "__main__":
    main()
