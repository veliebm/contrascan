"""
Calculate Cohen's d for each voxel in a whole brain ttest.

Created 10/26/2021 by Benjmain Velie.
"""
# Import external libraries and modules.
from os import PathLike
from typing import Tuple
from pathlib import Path
import nibabel
import numpy


def main(in_ttest: PathLike, out_cohen: PathLike, n_samples: int) -> None:
    """
    Execute the module.

    Load a ttest image, calculate cohen's d for each voxel, then save them.
    """
    func_image = nibabel.load(in_ttest)

    transformed_array = apply_cohens(
        a=func_image.dataobj[:, :, :, :, 1],
        n_samples=n_samples)

    save_array_to_nifti(
        affine=func_image.affine,
        array=transformed_array,
        out_path=out_cohen
    )


def apply_cohens(a: numpy.array, n_samples: int) -> numpy.array:
    """
    Calculate Cohen's d for each voxel in an array of t-test values.

    To convert from a t value to d: d = t/sqrt(n-1)
    (from https://www.researchgate.net/post/Is_cohen_d_equal_to_z_statistics_How_can_I_calculate_cohen_d_using_Z_scores)
    """  # nopep8
    cohens_results = a / numpy.sqrt(n_samples - 1)

    return cohens_results


def save_array_to_nifti(
    affine: Tuple,
    array: numpy.array,
    out_path: PathLike
) -> None:
    """
    Save a numpy array as a nifti image.
    """
    make_parent_directory(out_path)
    image = nibabel.Nifti1Image(array, affine)
    image.to_filename(out_path)


def make_parent_directory(path: PathLike) -> None:
    """
    Create the parent directory for a path.
    """
    Path(path).parent.mkdir(exist_ok=True, parents=True)


def _test_module() -> None:
    """
    Test this module.
    """
    kwargs = {'in_ttest': './processed/correlation_whole_brain_fisher_ttest/startvolume-1_variable-values_alpha_fisher_ttest.nii', 'out_cohen': './processed/correlation_whole_brain_fisher_cohens/startvolume-1_variable-values_alpha_fisher_cohens.nii', 'n_samples': 18}
    return main(**kwargs)
