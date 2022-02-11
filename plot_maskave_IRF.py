"""
Plot average IRFs computed with 3dmaskave.

Created 2/11/22 by Benjamin Velie.
"""

from os import PathLike
from pathlib import Path
from typing import List
import pandas
import matplotlib.pyplot as plt


def main():
    """
    Entrypoint of program.
    """
    irfs_paths = list(Path('processed/IRF_averages_in_rois/').glob('*.txt'))
    irfs_timeseries = [extract_timeseries(path) for path in irfs_paths]
    time_vector = [time*2 for time in range(len(irfs_timeseries[0]))]

    plot_irf(x_values=time_vector, timeseries_list=irfs_timeseries)


def plot_irf(x_values, timeseries_list) -> None:
    """
    Plot the IRFs
    """
    plt.rcParams.update({'font.size': 40})
    plt.rcParams.update({'axes.linewidth': 5})
    fig, ax = plt.subplots()
    plt.axhline(y=0, color='black', linestyle='--', linewidth=5)
    for timeseries in timeseries_list:
        ax.plot(x_values, timeseries, linewidth=10)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('BOLD Δ%')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(width=5, length=10)
    fig.set_size_inches(10, 9)
    plt.xticks(x_values[::2])
    plt.yticks([-.005, 0, .01])
    plt.tight_layout()
    plt.savefig('processed/plots/average_IRFs.png')


def extract_timeseries(maskave_path: PathLike) -> List[float]:
    """
    Extract a time series from a 3dmaskave output file.

    Args:
        maskave_path (PathLike): Path to 3dmask output file.

    Returns:
        List[float]: Vector containing time series in the file.
    """
    raw_data = pandas.read_csv(maskave_path, sep=" ", header=None)
    timeseries = raw_data[0].tolist()

    return timeseries


if __name__ == "__main__":
    main()
