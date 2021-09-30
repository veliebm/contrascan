"""
Functions to read matlab files from Python.
"""

from os import PathLike
import pandas
import scipy.io
import numpy


def get_amplitudes(mat_path: PathLike) -> pandas.Series:
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
