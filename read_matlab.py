"""
Functions to read matlab files from Python.
"""

from os import PathLike
from pathlib import Path
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
    Raises TypeError if the array is not 1xN or Nx1 or 1-dimensional.

    Parameters
    ----------
    array : numpy.array
        Should be 2 dimensional and be shaped either (1,N) or (N,1)

    Raises
    ------
    TypeError
        Array is not 2 dimensional or is not shaped correctly.
    """
    if array.ndim > 1:
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


def save_series_to_mat(data: pandas.Series, variable_name: str, out_path: PathLike, make_parent_dir: bool = True) -> None:
    """
    Save a pandas series to a mat file.

    Args:
        data (pandas.Series): Series to save.
        variable_name (str): What to name the series within the mat file.
        out_path (PathLike): Where to save the mat file.
        make_parent_dir (bool, optional): Whether to make the parent dir of the mat file or not. Defaults to True.
    """
    if make_parent_dir is True:
        make_parent_dir(out_path)

    numpy_data = data.to_numpy()

    _check_array_is_vector(numpy_data)

    scipy.io.savemat(out_path, {variable_name: numpy_data}, appendmat=False)


def make_parent_dir(path: PathLike) -> None:
    """
    Create the parent directory for a path.
    """
    Path(path).parent.mkdir(exist_ok=True, parents=True)
