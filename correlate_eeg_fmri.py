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
import scipy.io


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


def get_amplitudes(mat_path: PathLike) -> pandas.DataFrame:
    """
    Gets amplitudes from a mat file. Mat file must contain only an Nx1 vector.

    Parameters
    ----------
    mat_path : PathLike
        Path to a .mat file.
    
    Returns
    -------
    pandas.Series
        Series of numbers extracted from the mat file.
    
    Raises
    ------
    TypeError
        When the mat file contains data that isn't an Nx1 vector.
    """
    mat = load_mat_file(mat_path)
    
    _check_array_is_vector(mat)

    return pandas.Series(mat.reshape(-1))


def _check_array_is_vector(array: numpy.array) -> None:
    """
    Raises TypeError if the array is not 1xN or Nx1.

    Parameters
    ----------
    array : numpy.array
        Should be 2 dimensional and be shaped either (1,N) or (N,1)

    Raises
    ------
    TypeError
        Array is not 2 dimensional or is not shaped correctly.
    """
    shape = array.shape
    
    try:
        assert len(shape) == 2
    except AssertionError:
        raise TypeError(f"Array is constructed inappropriately: {shape}")
    
    try:
        assert shape[0] == 1 or shape[1] == 1
    except AssertionError:
        raise TypeError(f"Array is not vector shaped: {shape}")


def load_mat_file(path: PathLike) -> numpy.array:
    """
    Load a mat file into a DataFrame. Can only handle mat files with 1 variable inside.

    Parameters
    ----------
    path : PathLike
        Path to a mat file.
    
    Returns
    -------
    numpy.array
        Array of data inside the mat file.

    Raises
    ------
    ValueError
        When the mat file contains more than one variable.
    """
    mat = scipy.io.loadmat(path)

    values = []
    for key, value in mat.items():
        if "__" not in key:
            values.append(value)
    if len(values) > 1:
        raise ValueError(f"There is more than one variable in the mat file: {path}")
    value = values[0]

    return value
