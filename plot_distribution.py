"""
Plot a distribution of maxes and mins.

Created 12/16/2021 by Benjamin Velie.
"""
from pathlib import Path
from os import PathLike
from typing import Tuple
import matplotlib.pyplot as plt
import pandas


def main(in_distribution: PathLike, out_plot: PathLike, title: str, xlim: Tuple[float, float], ylim: Tuple[float, float]) -> None:
    """
    Plot a distribution of maxes and mins.
    """
    # Make figure.
    table = pandas.read_csv(in_distribution, index_col=0)
    fig = plt.plot(table["percentile"], table["value"])

    # Adjust the scale of the axes.
    plt.xlim(*xlim)
    plt.ylim(*ylim)

    # Label the plot.
    plt.xlabel("Percentile")
    plt.ylabel("Value")
    plt.title(title)

    # Save figure.
    make_parent_dir(out_plot)
    plt.savefig(out_plot)
    plt.close()


def make_parent_dir(path: PathLike) -> None:
    """
    Create the parent directory for a path.
    """
    Path(path).parent.mkdir(exist_ok=True, parents=True)


def _test_module():
    """
    Test this module.
    """
    kwargs = {'in_distribution': './processed/maxes_and_mins_distributions/startvolume-5_variable-amplitude_ssvep_mins_table.csv', 'out_plot': './processed/maxes_and_mins_distributions/startvolume-5_variable-amplitude_ssvep_mins_plot.png', 'title': 'Distribution of ssvep mins', 'xlim': (0, 100), 'ylim': (-0.1, 0)}
    main(**kwargs)
