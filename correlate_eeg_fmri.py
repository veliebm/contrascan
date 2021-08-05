#!/usr/bin/env python3
"""
Correlate our EEG data with our fMRI data.

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


def main(in_image_path: PathLike, in_eeg_path: PathLike, out_image_path: PathLike) -> None:
    """
    Correlate our EEG data with our fMRI data.

    The Oz electrode is number 20 in MatLab.
    """
    all_amplitudes = get_all_amplitudes(in_eeg_path)
    amplitudes = all_amplitudes[20]
    func_image = nibabel.load(in_image_path)
    Path(out_image_path).parent.mkdir(exist_ok=True, parents=True)
    correlate_subject(func_image, amplitudes, out_image_path)


def correlate_subject(func_image: nibabel.brikhead.AFNIImage, data_to_correlate: Iterable, out_path: PathLike) -> None:
    """
    Run a spearman correlation for a list of numbers for every voxel in a func image in 2mm MNI space.

    The vector of numbers must be the same length as the time domain of the func image.
    """

    # Get func image as dataframe. Trim so it's the same length as our EEG time series. Keys = coordinates, values = series.
    print(f"Building dataframe")
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
    trimmed_func = func_image.dataobj[:,:,:,:trim_volumes]

    future_dataframe = {}
    x_length, y_length, z_length, t_length = trimmed_func.shape
    for x in range(x_length):
        for y in range(y_length):
            for z in range(z_length):
                future_dataframe[(x, y, z)] = trimmed_func[x, y, z, :]
    func_dataframe = pandas.DataFrame(future_dataframe)

    return func_dataframe


def get_all_amplitudes(eeg_amplitudes_path: PathLike) -> pandas.DataFrame:
    """
    Returns the column of amplitudes for a channel from an eeg amplitudes tsv file.
    """
    all_amplitudes = pandas.read_csv(eeg_amplitudes_path, sep="\t", header=None).T
    all_amplitudes.columns += 1     # Make the channel numbers correspond to the MatLab channel numbers
    all_amplitudes.index *= 2

    return all_amplitudes
