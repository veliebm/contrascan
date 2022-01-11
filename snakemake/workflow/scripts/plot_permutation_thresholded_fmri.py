"""
The point of this snakemake script is to let you programmatically view and save screenshots of fMRI variance using permutation-based thresholds.

Created 12/7/2021 by Benjamin Velie.
veliebm@gmail.com
"""
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from os import PathLike
from typing import List
import nilearn.plotting
import nilearn.image
import nibabel
import matplotlib.pyplot as plt
import pandas


def main():
    """
    Entrypoint of the script.
    """
    conservative_threshold = get_conservative_threshold(
        percentile=float(snakemake.wildcards.percentile),
        mins_distribution=snakemake.input.mins_distribution,
        maxes_distribution=snakemake.input.maxes_distribution
    )

    make_plot(
        underlay_path=snakemake.input.underlay,
        overlay_path=snakemake.input.overlay,
        save_to_path=snakemake.output.plot,
        threshold=conservative_threshold,
        coordinates=snakemake.params.coordinates,
        title=f"{snakemake.wildcards.analysis} {snakemake.wildcards.variable} mean variance\ntime lag: {snakemake.wildcards.startvolume} volumes\nbaseline subtracted: {snakemake.wildcards.baselined}\ncritical values={conservative_threshold}\npercentile: {snakemake.wildcards.percentile}",
    )


def get_conservative_threshold(percentile: float, mins_distribution: PathLike, maxes_distribution: PathLike) -> float:
    """
    Get the variance thresholds from both the mins and maxes distribution. Return the more conservative one.

    Args:
        percentile (float): Percentile to read.
        mins_distribution (PathLike): Path to the mins distribution.
        maxes_distribution (PathLike): Path to the maxes distribution.

    Returns:
        float: The more conservative threshold of the two.
    """    
    upper_percentile = percentile
    lower_percentile = 100 - upper_percentile
    upper_threshold = get_threshold(maxes_distribution, upper_percentile)
    lower_threshold = get_threshold(mins_distribution, lower_percentile)
    conservative_threshold = min([upper_threshold, lower_threshold])

    print(f"2 available thresholds: {(upper_threshold, lower_threshold)}")
    print(f"Choosing {conservative_threshold}")

    return conservative_threshold


def get_threshold(path_to_table: PathLike, percentile: float) -> float:
    """
    Read a table of percentiles and return wanted threshold value after squaring it to turn it into variance.

    Args:
        path_to_table (PathLike): Path to distribution table to read.
        percentile (float): What percentile you want to get from the table.

    Returns:
        float: Threshold value.
    """
    table = pandas.read_csv(path_to_table, index_col="percentile")
    correlation_threshold = table["value"][percentile]
    variance_threshold = correlation_threshold ** 2

    return variance_threshold


def make_plot(underlay_path: PathLike, overlay_path: PathLike, save_to_path: PathLike, threshold: float, coordinates: List[int], title: str) -> None:
    """
    Plot an underlay and overlay.

    Args:
        underlay_path (PathLike): What image to use as an underlay.
        overlay_path (PathLike): What image to use as an overlay.
        save_to_path (PathLike): Where to save our plot.
        threshold (float): What to use as our threshold.
        coordinates (List[int]): Underlay coordinates where we want to make our plot.
        title (str): What title to put on the plot.
    """
    view = plt.figure(figsize=(10, 5))
    view = nilearn.plotting.plot_anat(underlay_path, cut_coords=coordinates, figure=view, draw_cross=False, annotate=False)
    overlay_image = load_as_nifti(overlay_path)
    view.add_overlay(overlay_image, colorbar=True, cmap="cold_hot", threshold=threshold)
    view.title(title)

    view.savefig(save_to_path)
    view.close()


def load_as_nifti(path_to_image: PathLike) -> nibabel.Nifti1Image:
    """
    Load an image and convert it into NIfTI format.
    
    For some reason, our AFNI headers don't play nice with nilearn's plotting functions.
    By transforming them into NIFTI's, we circumvent that problem.
    The affine, databj, header, and filemap are preserved.

    Args:
        path_to_image (PathLike): Path to an AFNI image.

    Returns:
        nibabel.Nifti1Image: A nibabel NIfTI image object containing the image data.
    """
    raw_image = nibabel.load(path_to_image)
    nifti_image = nibabel.Nifti1Image(raw_image.dataobj, raw_image.affine, header=raw_image.header, file_map=raw_image.file_map)

    return nifti_image


if __name__ == "__main__":
    main()
