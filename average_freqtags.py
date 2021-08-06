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


def main(in_tsvs: List[PathLike], out_tsv: PathLike) -> None:
    """
    Average our freqtag stuff across subjects.
    """
    tables = [read_tsv(tsv) for tsv in in_tsvs]

    summed_table = tables[0]
    for table in tables[1:]:
        summed_table += table

    averaged_table = summed_table / len(tables)

    save_tsv(averaged_table, out_tsv)


def save_tsv(table: pandas.DataFrame, out_path: PathLike) -> None:
    """
    Save a TSV file to disk.
    """
    Path(out_path).parent.mkdir(exist_ok=True, parents=True)
    table.to_csv(out_path, sep='\t', index=None)


def read_tsv(tsv_path: PathLike) -> pandas.DataFrame:
    """
    Read a TSV file. Returns a DataFrame.
    """
    tsv_info = pandas.read_table(
        tsv_path,
        sep="\t",
        index_col=False,
        header=None,
    ).T

    tsv_info.columns += 1     # Make the channel numbers correspond to the MatLab channel numbers

    return tsv_info
