#!/usr/bin/env python3
"""
Correlate a list of numbers with each voxel in a functional image.

The moment we've all been waiting for!

Created 7/28/2021 by Ben Velie, veliebm@ufl.edu
"""
# Import external libraries and modules.
from os import PathLike
from typing import Iterable
from pathlib import Path
import pandas
import nibabel
import numpy

# Import CSEA stuff
from read_matlab import get_amplitudes


def main(in_image_path: PathLike, in_eeg_path: PathLike, out_image_path: PathLike) -> None:
    """
    Correlate our EEG data with our fMRI data.

    The Oz electrode is number 20 in MatLab.
    """
    amplitudes = get_amplitudes(in_eeg_path)
    func_image = nibabel.load(in_image_path)
    Path(out_image_path).parent.mkdir(exist_ok=True, parents=True)

    correlate_subject(func_image, amplitudes, out_image_path)


def correlate_subject(func_image: nibabel.brikhead.AFNIImage, data_to_correlate: Iterable, out_path: PathLike) -> None:
    """
    Run a spearman correlation for a list of numbers for every voxel in a func image in 2mm MNI space.

    The vector of numbers must be the same length as the time domain of the func image.
    """

    # Get func image as dataframe. Trim so it's the same length as our EEG time series. Keys = coordinates, values = series.
    print(f"Converting func image into a DataFrame")
    func_dataframe = convert_to_dataframe(func_image, trim_volumes=len(data_to_correlate))

    # For each voxel, record correlation value.
    print("Running Spearman correlation")
    correlation_results = func_dataframe.corrwith(data_to_correlate, method='spearman')

    # Convert correlation results into a numpy array.
    print(f"Correlations complete! Converting into a numpy array")
    x_length, y_length, z_length, __ = func_image.dataobj[:,:,:,:len(data_to_correlate)].shape
    correlation_array = numpy.empty((x_length, y_length, z_length))
    for point, correlation in correlation_results.iteritems():
        x, y, z = point
        correlation_array[x,y,z] = correlation

    # Convert numpy array into a nibabel image.
    print(f"Saving image to {out_path}")
    affine = func_image.affine
    correlation_image = nibabel.Nifti1Image(correlation_array, affine)
    correlation_image.to_filename(out_path)


def convert_to_dataframe(func_image: nibabel.brikhead.AFNIImage, trim_volumes: int) -> pandas.DataFrame:
    """
    Convert a func image to a dataframe. Trim time axis to specified point.
    """
    print(f"Func shape: {func_image.dataobj.shape}")
    print(f"Length of data to correlate with func: {trim_volumes}")
    trimmed_func = func_image.dataobj[:,:,:,:trim_volumes]
    print(f"New func shape: {trimmed_func.shape}")

    future_dataframe = {}
    x_length, y_length, z_length, __ = trimmed_func.shape
    for x in range(x_length):
        for y in range(y_length):
            for z in range(z_length):
                future_dataframe[(x, y, z)] = trimmed_func[x, y, z, :]
    func_dataframe = pandas.DataFrame(future_dataframe)

    return func_dataframe


def _test_module() -> None:
    """
    Test this module.
    """
    return main(**{'in_image_path': './processed/func_final_trim/sub-104_startvolume-1_func+tlrc.HEAD', 'in_eeg_path': './processed/eeg_sliding_sliding_window/sub-104_frequency-12_oz_amplitudes.mat', 'out_image_path': './processed/correlation_whole_brain/sub-104_frequency-12_startvolume-1_correlation.nii'})
