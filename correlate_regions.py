#!/usr/bin/env python3
"""
Correlate the time series for a microregion with its Oz data.

Created 8/23/2021 by Ben Velie, veliebm@ufl.edu
"""

import numpy
from os import PathLike
from typing import Dict, List
import pandas
import argparse
import matplotlib.pyplot as plt
from pathlib import Path
import scipy.io


def main(load_microROI_from: PathLike, load_amplitudes_from: PathLike, save_scatter_to: PathLike, save_table_to: PathLike, save_spearman_to: PathLike) -> None:
    """
    Correlate an average microROI time series with its corresponding Oz amplitudes.
    """
    Path(save_spearman_to).parent.mkdir(parents=True, exist_ok=True)

    amplitudes = get_amplitudes(load_amplitudes_from).tolist()
    microROI = read_average(load_microROI_from)[:len(amplitudes)]

    table = make_table(
        dict(
            microROI_signal=microROI,
            amplitude=amplitudes,
        ),
        out_path=save_table_to,
    )
    get_spearman(
        table["microROI_signal"].tolist(),
        table["amplitude"].tolist(),
        out_path=save_spearman_to,
    )
    make_scatter_plot(
        table=table,
        out_path=save_scatter_to,
    )
 

def make_table(dictionary: Dict, out_path: PathLike) -> pandas.DataFrame:
    """
    Make a table and also save it to disk.
    """
    table = pandas.DataFrame(dictionary)
    table.to_csv(out_path)

    return table


def make_scatter_plot(table: pandas.DataFrame, out_path: PathLike) -> None:
    """
    Make a scatter plot of our data. Save plot and table used to produce the plot.
    """
    table.plot(x="microROI_signal", y="amplitude", kind="scatter")

    x = table["microROI_signal"]
    y = table["amplitude"]
    m, b = numpy.polyfit(x=x, y=y, deg=1)

    label = f"Line of best fit: y = {m:.3f}x + {b:.3f}"
    print(f"\n{label}\n")

    plt.plot(x, m*x + b, label=label)
    plt.legend()    
    plt.savefig(out_path)


def get_spearman(list1: List[float], list2: List[float], out_path: PathLike) -> float:
    """
    Runs a spearman correlation on 2 lists of data.
    """
    table = pandas.DataFrame([list1, list2]).T
    correlation_results = table.corr(method="spearman")
    correlation = correlation_results[0][1]
    print(f" \nCorrelation: \n{correlation} \n")

    with open(out_path, "w") as io:
        io.write(str(correlation))

    return correlation


def read_average(path: PathLike) -> None:
    """
    Get averages from the target file.
    """
    with open(path, "r") as io:
        lines = io.readlines()

    averages = [float(line.split()[0]) for line in lines]

    return averages


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


if __name__ == "__main__":
    """
    This code only runs if this file is called as a script.
    """
    parser = argparse.ArgumentParser(description="Correlate the time series for a microregion with its Oz data.")

    parser.add_argument("--load_microROI_from", required=True, help="Where to get the average time series of a microROI")
    parser.add_argument("--load_amplitudes_from", required=True, help="Where to get out moving moving window results for said microROI")
    parser.add_argument("--save_scatter_to", required=True, help="Where to save our scatter plot")
    parser.add_argument("--save_table_to", required=True, help="Where to save the table of data containing our microROIs and amplitudes")
    parser.add_argument("--save_spearman_to", required=True, help="Where to save our Spearman results")

    parsed_args = parser.parse_args()
    parsed_args_dict = vars(parsed_args)
    main(**parsed_args_dict)
