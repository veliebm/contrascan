"""
Calculate the distribution of a percentile from each image within a list of fMRI images.

Created 1/17/22 by Benjamin Velie.
veliebm@gmail.com
"""
from os import PathLike
from typing import List
import nibabel
import numpy
import pandas
from pathlib import Path


def main(in_permutations: List[PathLike], out_distribution: PathLike, percentile: float) -> None:
    """
    Get the distribution of a percentile from a list of permutation images.

    Args:
        in_permutations (List[PathLike]): List of paths to permutation images.
        out_distribution (PathLike): Table of percentiles from those images.
        percentile (PathLike): What percentile to extract from each image.
    """
    # Read images.
    permutation_images = [nibabel.load(path) for path in in_permutations]

    # Get data of images as list of numpy arrays.
    try:
        permutation_arrays = [image.dataobj[:, :, :, 0] for image in permutation_images]
    except IndexError:
        print("At least one image isn't 4d. Trying to read as 3d images.")
        permutation_arrays = [image.dataobj[:, :, :] for image in permutation_images]

    # Combine list of arrays into one big array.
    all_permutations_array = numpy.stack(permutation_arrays, 3)

    # Flatten array into 2 dimensions: 1 is the index of the image, 0 is the data of that image.
    num_voxels_in_each_image = all_permutations_array.shape[0] * all_permutations_array.shape[1] * all_permutations_array.shape[2]
    array_length = all_permutations_array.shape[3]
    flattened_array = numpy.reshape(all_permutations_array, (num_voxels_in_each_image, array_length))

    # Calculate percentile for each image in big array.
    percentile_array = numpy.percentile(flattened_array, percentile, axis=0)

    # Store arrays into a DataFrame.
    quantile_dataframe = pandas.DataFrame(data=percentile_array.T, index=in_permutations, columns=["value"])

    # Rank data.
    num_values = len(quantile_dataframe["value"])
    quantile_dataframe["percentile"] = (quantile_dataframe["value"].rank() - 1) / num_values * 100

    # Sort data.
    quantile_dataframe.sort_values("percentile", ascending=False, inplace=True)

    # Save results to disk.
    make_parent_dir(out_distribution)
    quantile_dataframe.to_csv(out_distribution)


def _test_module():
    """
    Test this module.
    """
    kwargs = {}
    main(**kwargs)


def make_parent_dir(path: PathLike) -> None:
    """
    Create the parent directory for a path.
    """
    Path(path).parent.mkdir(exist_ok=True, parents=True)
