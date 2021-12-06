"""
Use our permutations to calculate the p-value of our correlation.
"""
from os import PathLike
from pathlib import Path
from typing import Tuple
import nibabel
import scipy.stats
import numpy


def main(catenated_permutations: PathLike, actual_correlation: PathLike, out_percentiles: PathLike) -> None:
    """
    Use our permutation correlations to calculate the percentile of our actual correlation.

    Args:
        catenated_permutations (PathLike): Path to an image of catenated permutation correlations.
        actual_correlation (PathLike): Path to an image of actual correlations.
        out_percentiles (PathLike): Where to write our percentiles image to.
    """
    make_parent_dir(out_percentiles)

    permutations_image = nibabel.load(catenated_permutations)
    actual_image = nibabel.load(actual_correlation)

    combined_array = get_combined_array(permutations_image, actual_image)
    percentile_array = get_percentiles(combined_array, axis=3)
    percentiles_of_actual_correlation = percentile_array[:, :, :, -1]

    save_array_to_nifti(actual_image.affine, percentiles_of_actual_correlation, out_percentiles)


def get_combined_array(permutations_image: nibabel.Nifti1Image, actual_image: nibabel.Nifti1Image) -> numpy.array:
    """
    Returns an array containing the combined data from the permutations and the actual image.

    Final volume in 4d array contains the acual correlation data.

    Args:
        permutations_image (nibabel.Nifti1Image): [description]
        actual_image (nibabel.Nifti1Image): [description]

    Returns:
        numpy.array: [description]
    """
    permutations_array = numpy.asarray(permutations_image.dataobj)[:, :, :, :]
    actual_3d_array = numpy.asarray(actual_image.dataobj)[:, :, :, 0]
    actual_4d_array = numpy.expand_dims(actual_3d_array, 3)
    combined_array = numpy.concatenate((permutations_array, actual_4d_array), 3)

    return combined_array


def get_percentiles(some_array: numpy.array, axis: int) -> numpy.array:
    """
    Get the percentiles of an array along a dimension.

    Percentiles are in terms of the following: percent of permutations (including the actual data) below the value of a voxel.
    Percentile formula used: Percentile = (number of values below score) ÷ (total number of scores) x 100

    Args:
        some_array (numpy.array): Array to rank percentile-wise.
        axis (int): What dimension to rank along.

    Returns:
        numpy.array: Array containing percentile results.
    """
    length = some_array.shape[axis]
    ranked_high_to_low_array = scipy.stats.rankdata(some_array, axis=axis)
    ranked_low_to_high_array = length - ranked_high_to_low_array
    percentile_array = (ranked_low_to_high_array - 1) / length * 100

    return percentile_array


def make_parent_dir(path: PathLike) -> None:
    """
    Create the parent directory for a path.
    """
    Path(path).parent.mkdir(exist_ok=True, parents=True)


def save_array_to_nifti(affine: Tuple, array: numpy.array, out_path: PathLike) -> None:
    """
    Save a numpy array as a nifti image.

    Args:
        affine (Tuple): The affine of the image.
        array (numpy.array): Array of data to save.
        out_path (PathLike): Where to write the image. Can write gunziped niftis if you name the file appropriately.
    """
    image = nibabel.Nifti1Image(array, affine)
    image.to_filename(out_path)


def _test_module() -> None:
    """
    Test this module.
    """
    kwargs = {'catenated_permutations': './processed/catenated_permutations/variable-amplitude_startvolume-5_analysis-ssvep_catenated.nii.gz', 'actual_correlation': './processed/correlation_whole_brain_ttest/startvolume-5_variable-amplitudes_improved_correlations_ttest+tlrc.HEAD', 'out_percentiles': './processed/pvalues_thresholding/variable-amplitude_startvolume-5_analysis-ssvep_pvalues.nii.gz'}
    main(**kwargs)
