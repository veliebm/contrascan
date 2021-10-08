"""
Get the maximum amplitude of each trial for each voxel in a functional image.

Created 10/8/2021 by Benjamin Velie
veliebm@ufl.edu
"""

from os import PathLike
from pathlib import Path
import numpy
import nibabel


def main(onsets_path: PathLike, func_path: PathLike, trials_path: PathLike) -> None:
    """
    Get the maximum amplitude of each trial for each voxel in a functional image.
    """
    # Read func image.
    func_image = nibabel.load(func_path)

    # Get volumes containing onsets.
    onsets = read_onsets(onsets_path)
    onset_volumes = get_onset_volumes(onsets)

    # Get array containing amplitudes of each trial for each voxel
    amplitudes = get_trial_amplitudes(func_image, onset_volumes)

    # Convert numpy array into a nibabel image.
    Path(trials_path).parent.mkdir(exist_ok=True, parents=True)
    affine = func_image.affine
    trials_image = nibabel.Nifti1Image(amplitudes, affine)
    trials_image.to_filename(trials_path)


def get_trial_amplitudes(func_image, onset_volumes: numpy.array) -> numpy.array:
    """
    Get array containing amplitudes of each trial for each voxel
    """
    # Initialize empty array to store our amplitudes in.
    volume_shape = func_image.dataobj[:,:,:,0].shape
    number_of_volumes = len(onset_volumes)
    shape = (*volume_shape, number_of_volumes)
    trial_amplitudes = numpy.zeros(shape=shape)

    # For each volume containing an onset, slice out that and the following 4 volumes.
    for i, volume in enumerate(onset_volumes):
        baseline = func_image.dataobj[:,:,:,volume]
        trial = func_image.dataobj[:,:,:,volume+1:volume+5]

        # For each voxel in the 4 volumes after the 1st, find the maximum amplitude in that voxel.
        trial_max = numpy.amax(trial, axis=3)

        # Subtract the baseline 1st volume from the maximum that we found.
        difference_from_baseline = trial_max - baseline

        # Store our calculation in the empty array we made at the top.
        trial_amplitudes[:,:,:,i] = difference_from_baseline
    
    return trial_amplitudes


def get_onset_volumes(onsets: numpy.array) -> numpy.array:
    """
    Calculate the volumes containing each onset.
    """
    # Round begnning and end of each trial to 2s. Divide by 2 to get start and end volumes to replace.
    volumes = onsets/2
    volumes = volumes.astype(int)

    return volumes


def read_onsets(onsets_path: PathLike) -> numpy.array:
    """
    Read a text file containing onsets. Return it as a list of floats.
    """
    with open(onsets_path, "r") as io:
        lines = io.readlines()
        onsets = numpy.array([float(line) for line in lines])

    return onsets


def _test_module() -> None:
    """
    Test this module.
    """
    kwargs = {'onsets_path': './processed/afniproc/sub-113/onsets.tsv', 'func_path': './processed/func_resample/sub-113_func_resampled+tlrc.HEAD', 'trials_path': './processed/func_trial_amplitudes/sub-113_trials.nii'}
    main(**kwargs)
