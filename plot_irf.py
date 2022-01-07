"""
Make a nice plot of a representative occipital IRF.

Created 1/6/22 by Benjamin Velie.
"""

from os import PathLike, path
import matplotlib.pyplot as plt
import nibabel
from pathlib import Path
import numpy


def main():
    """
    Entrypoint of program.
    """
    path_to_irf = Path("processed/IRF_mean/IRF_mean+tlrc.HEAD")
    representative_voxel_timeseries = get_representative_time_series(path_to_irf)
    time_vector = [time*2 for time in range(len(representative_voxel_timeseries))]
    plot_irf(x_values=time_vector, y_values=representative_voxel_timeseries)


def plot_irf(x_values, y_values) -> None:
    """
    Plot the IRF.

    Args:
        x_values ([type]): [description]
        y_values ([type]): [description]
    """
    fig, ax = plt.subplots()
    ax.plot(x_values, y_values, linewidth=5, color='black')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('BOLD Δ%')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.set_size_inches(16, 9)
    plt.xticks(numpy.arange(min(x_values), max(x_values)+1, 2))
    plt.yticks((-.2, 0, .4))
    plt.rcParams.update({'font.size': 30})
    plt.rcParams['axes.linewidth'] = 5
    plt.tight_layout()
    plt.savefig('processed/plots/representative_IRF.png')


def get_representative_time_series(path_to_irf: PathLike) -> numpy.array:
    """
    Extract a time series from a representative voxel in the occipital pole.

    Args:
        path_to_irf (PathLike): Path to the IRF image.

    Returns:
        numpy.array: Vector containing time series of voxel from occipital pole.
    """
    irf_image = nibabel.load(path_to_irf)
    representative_voxel_coordinates = (29, 10, 33)
    representative_voxel_timeseries = irf_image.dataobj[representative_voxel_coordinates]

    return representative_voxel_timeseries


if __name__ == "__main__":
    main()
