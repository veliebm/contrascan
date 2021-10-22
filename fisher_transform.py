#!/usr/bin/env python3
"""
Apply the Fisher Z transform to our whole brain correlation so we can do more valid t-tests.

Created 10/22/2021 by Ben Velie, veliebm@ufl.edu
"""
# Import external libraries and modules.
from os import PathLike
from typing import Tuple
from pathlib import Path
import nibabel
import numpy


def main(in_correlation_path: PathLike, out_fisher_path: PathLike) -> None:
    """
    Execute the module.
    """
    # Load data.
    func_image = nibabel.load(in_correlation_path)

    # Transform data.
    fisher_array = apply_fisher(a=func_image.dataobj[:,:,:,0])

    # Save data.
    save_array_to_nifti(affine=func_image.affine, array=fisher_array, out_path=out_fisher_path)


def apply_fisher(a: numpy.array) -> numpy.array:
    """
    Run a spearman correlation for a list of numbers for every voxel in a func image in 2mm MNI space.

    The vector of numbers must be the same length as the time domain of the func image.
    """
    fisher_results = numpy.arctanh(a)

    return fisher_results


def save_array_to_nifti(affine: Tuple, array: numpy.array, out_path: PathLike) -> None:
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
    return main(in_correlation_path="processed/correlation_whole_brain/sub-110_startvolume-8_variable-amplitudes_improved.nii", out_fisher_path="TESTFISHER.nii")
