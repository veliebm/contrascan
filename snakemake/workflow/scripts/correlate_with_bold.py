"""
Correlate a list of numbers with each voxel in a functional image.

Created 2/23/22 by Ben Velie.
"""
from typing import Iterable
from scipy.stats import spearmanr
import nibabel
import numpy


def main() -> None:
    """
    Entrypoint of script.
    """
    # Load data.
    vector = numpy.genfromtxt(snakemake.input.vector, delimiter=',')
    func_image = nibabel.load(snakemake.input.bold)

    # Correlate data.
    correlation_results = get_correlation(func_image, vector)

    # Save results.
    image = nibabel.Nifti1Image(correlation_results, func_image.affine)
    nibabel.save(image, snakemake.output.correlations)


def get_correlation(func_image: nibabel.Nifti1Image, vector: Iterable) -> numpy.array:
    """
    Speaman correlate a vector of numbers with each voxel in a func image.

    Removes volumes from the end of functional image to make it the same length as the vector.
    """
    trimmed_array = func_image.dataobj[:, :, :, :len(vector)]
    print(f"Trimmed image from {func_image.dataobj.shape} to {trimmed_array.shape}")
    print("Running Spearman correlation")
    correlation_results = numpy.apply_along_axis(spearmanr, 3, trimmed_array, b=vector)
    correlation_results_with_p_vals_removed = correlation_results[:, :, :, 0]

    return correlation_results_with_p_vals_removed


if __name__ == "__main__":
    main()
