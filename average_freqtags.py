#!/usr/bin/env python3
"""
Average our freqtag stuff across subjects.

Created 8/6/21 by Benjamin Velie.
veliebm@ufl.edu
"""

from os import PathLike
from typing import List
import pandas
from pathlib import Path


def main(in_csvs: List[PathLike], out_csv: PathLike) -> None:
    """
    Average our freqtag stuff across subjects.
    """
    tables = [read_csv(csv) for csv in in_csvs]

    summed_table = tables[0]
    for table in tables[1:]:
        summed_table += table

    averaged_table = summed_table / len(tables)

    save_tsv(averaged_table, out_csv)


def save_tsv(table: pandas.DataFrame, out_path: PathLike) -> None:
    """
    Save a TSV file to disk.
    """
    Path(out_path).parent.mkdir(exist_ok=True, parents=True)
    table.to_csv(out_path, sep='\t', index=None)


def read_csv(csv_path: PathLike) -> pandas.DataFrame:
    """
    Read a csv file. Returns a DataFrame.
    """
    csv_info = pandas.read_table(
        csv_path,
        sep=",",
        index_col=False,
        header=None,
    ).T

    csv_info.columns += 1     # Make the channel numbers correspond to the MatLab channel numbers

    return csv_info
