"""
Get the mean of a list of mean correlation files.

Created 11/30/2021 by Benjamin Velie.
veliebm@gmail.com
"""
import pandas
from os import PathLike
from typing import List
from pathlib import Path


def main(from_text_files: List[PathLike], to_table: PathLike, to_result: PathLike) -> None:
    """
    Read text files containing correlation averages, then average those averages.

    Args:
        from_text_files (List[PathLike]): List of text files to read.
        to_table (PathLike): Where to write table of correlations we find. (Use to check our work if you like.)
        to_result (PathLike): Where to write the average average to.
    """
    # Make out dirs.
    for path in (to_table, to_result):
        Path(path).parent.mkdir(exist_ok=True, parents=True)

    # Extract correlation average from each text file.
    correlations = {source_file: extract_correlation(source_file) for source_file in from_text_files}

    # Concatenate the correlations into a summary pandas table.
    labelled_correlations = dict()
    labelled_correlations["filename"] = correlations.keys()
    labelled_correlations["correlation"] = correlations.values()
    correlations_table = pandas.DataFrame(labelled_correlations)

    # Calculate mean.
    results = correlations_table["correlation"].mean()

    # Save results.
    results_dict = {"mean": results}
    results_table = pandas.DataFrame(results_dict, index=[0])
    results_table.to_csv(to_result)

    # Save summary table.
    correlations_table.to_csv(to_table)


def extract_correlation(path_to_file: PathLike) -> float:
    """
    Extract a correlation from a text file.

    Reads the first line of the text file and uses whatever number is finds there.
    """
    with open(path_to_file, "r") as io:
        lines = io.readlines()
        average = float(lines[0].split()[0])
        return average
