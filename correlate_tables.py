#!/usr/bin/env python3
"""
Given tables of data containing 2 variables of interest, combine all tables and then compute one big spearman correlation.

Created 10/4/2021 by Benjamin Velie.
veliebm@ufl.edu
"""

import numpy
from os import PathLike
from typing import Dict, List
import pandas
import matplotlib.pyplot as plt
from pathlib import Path


def main(load_tables_from: List[PathLike], save_scatter_to: PathLike, save_table_to: PathLike, save_spearman_to: PathLike) -> None:
    """
    Catenate our microROI+amplitude tables together, then run a correlation on the uber-table.
    """
    Path(save_spearman_to).parent.mkdir(parents=True, exist_ok=True)

    tables = [load_table(path) for path in load_tables_from]
    
    table = catenate_tables(tables, save_table_to)
    get_spearman(table.iloc[:, 0].tolist(), table.iloc[:, 1].tolist(), save_spearman_to)
    make_scatter_plot(table, save_scatter_to)


def load_table(in_path: PathLike) -> pandas.DataFrame:
    """
    Load a table we've stored as a CSV file.
    """
    table = pandas.read_csv(in_path, index_col=0)

    return table

    
def catenate_tables(tables: List[pandas.DataFrame], out_path: PathLike) -> pandas.DataFrame:
    """
    Make a table and also save it to disk.
    """
    table = pandas.concat(tables)
    table.to_csv(out_path)

    return table


def make_scatter_plot(table: pandas.DataFrame, out_path: PathLike) -> None:
    """
    Make a scatter plot of our data. Save plot and table used to produce the plot.
    """
    x_name = table.columns[0]
    y_name = table.columns[1]
    table.plot(x=x_name, y=y_name, kind="scatter")

    x = table[x_name]
    y = table[y_name]
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
    print(table)
    print(f" \nCorrelation: \n{correlation} \n")

    with open(out_path, "w") as io:
        io.write(str(correlation))

    return correlation
