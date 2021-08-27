#!/usr/bin/env python3
"""
We've already correlated each microROI with its Oz data. Now let's correlate ALL microROIs
with ALL Oz data!

Created 8/27/2021 by Benjamin Velie.
veliebm@ufl.edu
"""

import numpy
from os import PathLike
from typing import Dict, List
import pandas
import argparse
import matplotlib.pyplot as plt
from pathlib import Path


def main(load_tables_from: List[PathLike], save_scatter_to: PathLike, save_table_to: PathLike, save_spearman_to: PathLike) -> None:
    """
    Catenate our microROI+amplitude tables together, then run a correlation on the uber-table.
    """
    Path(save_spearman_to).parent.mkdir(parents=True, exist_ok=True)

    tables = [load_table(path) for path in load_tables_from]
    
    table = catenate_tables(tables, save_table_to)
    get_spearman(table["microROI_signal"].tolist(), table["amplitude"].tolist(), save_spearman_to)
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


if __name__ == "__main__":
    """
    This code only runs if this file is called as a script.
    """
    parser = argparse.ArgumentParser(description="Correlate the time series for a microregion with its Oz data.")

    parser.add_argument("--load_tables_from", required=True, nargs="+", help="Where to get each of our microROI+amplitude tables")
    parser.add_argument("--save_scatter_to", required=True, help="Where to save our scatter plot")
    parser.add_argument("--save_table_to", required=True, help="Where to save the table of data we make from our tiny tables")
    parser.add_argument("--save_spearman_to", required=True, help="Where to save our Spearman results")

    parsed_args = parser.parse_args()
    parsed_args_dict = vars(parsed_args)
    main(**parsed_args_dict)
