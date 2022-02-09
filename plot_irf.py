"""
Make a nice plot of a representative occipital IRF.

Created 1/6/22 by Benjamin Velie.
"""

from os import PathLike, path
import matplotlib.pyplot as plt
import nibabel
import numpy


def main():
    """
    Entrypoint of program.
    """
    calcarine_timeseries = [0, -.07243, -.088155, -.171879, -.143361, .037872, .066945, .032208, .013069, 0]
    occipital_timeseries = [0, -.051867, .238449, .579762, .771096, .711357, .295275, -.034379, 0.182821, 0]
    time_vector = [time*2 for time in range(len(calcarine_timeseries))]

    plot_irf(x_values=time_vector, calcarine_values=calcarine_timeseries, occipital_values=occipital_timeseries)


def plot_irf(x_values, calcarine_values, occipital_values) -> None:
    """
    Plot the IRFs
    """
    plt.rcParams.update({'font.size': 40})
    plt.rcParams.update({'axes.linewidth': 5})
    fig, ax = plt.subplots()
    plt.axhline(y=0, color='black', linestyle='--', linewidth=5)
    ax.plot(x_values, calcarine_values, linewidth=10, color='blue')
    ax.plot(x_values, calcarine_values, color="blue", markevery=[4], marker="o", ls="", label="points", ms=20)
    ax.plot(x_values, occipital_values, linewidth=10, color='red')
    ax.plot(x_values, occipital_values, color="red", markevery=[4], marker="o", ls="", label="points", ms=20)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('BOLD Δ%')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(width=5, length=10)
    fig.set_size_inches(10, 9)
    plt.xticks(x_values[::2])
    plt.yticks([-.2, 0, .8])
    plt.tight_layout()
    plt.savefig('processed/plots/representative_IRFs.png')


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
