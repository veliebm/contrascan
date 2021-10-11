#!/usr/bin/env python3
"""
Correlate a list of numbers with each voxel in a functional image.

The moment we've all been waiting for!

Created 7/28/2021 by Ben Velie, veliebm@ufl.edu
"""
# Import external libraries and modules.
from os import PathLike
from typing import Iterable, Tuple
from pathlib import Path
from scipy.stats import spearmanr
import nibabel
import numpy

# Import CSEA stuff
from read_matlab import get_amplitudes


def main(in_image_path: PathLike, in_eeg_path: PathLike, out_image_path: PathLike) -> None:
    """
    Correlate our EEG data with our fMRI data.

    The Oz electrode is number 20 in MatLab.
    """
    # Load data.
    amplitudes = get_amplitudes(in_eeg_path)
    func_image = nibabel.load(in_image_path)

    # Make out directory if it doesn't exist.
    Path(out_image_path).parent.mkdir(exist_ok=True, parents=True)

    # Correlate data.
    correlate_subject(func_image, amplitudes, out_image_path)


def correlate_subject(func_image: nibabel.brikhead.AFNIImage, data_to_correlate: Iterable, out_path: PathLike) -> None:
    """
    Run a spearman correlation for a list of numbers for every voxel in a func image in 2mm MNI space.

    The vector of numbers must be the same length as the time domain of the func image.
    """
    # Get func image as dataframe. Trim so it's the same length as our EEG time series. Keys = coordinates, values = series.
    print(f"Trimming array to length {len(data_to_correlate)}")
    trimmed_array = func_image.dataobj[:,:,:,:len(data_to_correlate)]

    # For each voxel, record correlation value.
    print("Running Spearman correlation")
    correlation_results = numpy.apply_along_axis(spearmanr, 3, trimmed_array, b=data_to_correlate)

    # Convert numpy array into a nibabel image.
    save_array_to_nifti(func_image.affine, correlation_results, out_path)


def save_array_to_nifti(affine: Tuple, array: numpy.array, out_path: PathLike) -> None:
    """
    Save a numpy array as a nifti image.
    """
    image = nibabel.Nifti1Image(array, affine)
    image.to_filename(out_path)


def _test_module() -> None:
    """
    Test this module.
    """
    return main(**{'in_image_path': './processed/func_final_trim/sub-104_startvolume-1_func+tlrc.HEAD', 'in_eeg_path': './processed/eeg_sliding_sliding_window/sub-104_frequency-12_oz_amplitudes.mat', 'out_image_path': './processed/correlation_whole_brain/sub-104_frequency-12_startvolume-1_correlation.nii'})
