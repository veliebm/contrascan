"""
The point of this snakemake script is to let you programmatically view and save screenshots of fMRI data.

Created 12/7/2021 by Benjamin Velie.
veliebm@gmail.com
"""
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from os import PathLike
from typing import List, Tuple
import nilearn.plotting
import nilearn.image
import nibabel
import numpy


def main():
    make_plot(
        underlay_path=snakemake.input.underlay,
        overlay_path=snakemake.input.overlay,
        save_to_path=snakemake.output.plot,
        threshold=snakemake.params.threshold,
        coordinates=snakemake.params.coordinates,
        title=snakemake.params.title,
    )


def make_plot(underlay_path: PathLike, overlay_path: PathLike, save_to_path: PathLike, threshold: Tuple[float, float], coordinates: List[int], title: str) -> None:
    """
    Plot an underlay and overlay.
    """
    view = nilearn.plotting.plot_anat(underlay_path, cut_coords=coordinates, draw_cross=False, annotate=False)
    overlay_image = load_as_nifti(overlay_path)
    thresholded_overlay = apply_threshold(overlay_image, threshold[0], threshold[1])
    view.add_overlay(thresholded_overlay, colorbar=True, cmap="cold_hot")
    view.title(title)

    view.savefig(save_to_path)
    view.close()


def apply_threshold(image: nibabel.Nifti1Image, bottom_threshold: float, top_threshold: float) -> nibabel.Nifti1Image:
    """
    Apply a bottom threshold and top threshold to an image.
    """
    image_array = numpy.array(image.dataobj)
    image_array[(image_array>bottom_threshold) & (image_array<top_threshold)] = numpy.nan
    thresholded_image = nibabel.Nifti1Image(image_array, image.affine, image.header)
    return thresholded_image


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
