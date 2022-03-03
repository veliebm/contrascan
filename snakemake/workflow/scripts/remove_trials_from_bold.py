"""
Remove trials from a BOLD image.

Created 9/30/2022 By Benjamin Velie.
"""
import pandas
import nibabel
import numpy
import math


def main() -> None:
    """
    Entrypoint of the script.
    """
    # Load data.
    input_image = nibabel.load(snakemake.input.image)
    events = pandas.read_csv(snakemake.input.events, sep="\t")

    # Get volumes to remove as a 1D True/False array.
    volumes_to_remove = get_volumes_to_remove(events, input_image.shape[3])

    # Mask nibabel image using the True/False array.
    mask = numpy.tile(volumes_to_remove, (input_image.shape[0:3]))
    masked_array = numpy.ma.array(input_image.get_fdata(), mask=mask)
    non_trial_values = masked_array[~masked_array.mask]
    volume_size = math.prod(input_image.shape[0:3])
    non_trial_t_length = int(len(non_trial_values) / volume_size)
    trimmed_array = non_trial_values.reshape(input_image.shape[0:3] + tuple([non_trial_t_length]))

    # Save data.
    out_image = nibabel.Nifti1Image(trimmed_array, input_image.affine)
    out_image.to_filename(snakemake.output.image)
    numpy.savetxt(snakemake.output.quality_control, volumes_to_remove)


def get_volumes_to_remove(events: pandas.DataFrame, bold_length: int) -> numpy.array:
    """
    Returns a numpy array containing whether a trial was present for each volume.

    Args:
        events (pandas.DataFrame): Table of event information.
        bold_length (int): Length of the BOLD image in the Z direction.

    Returns:
        numpy.array: 1D array. Contains True for volumes in which a trial is present, else False.
    """
    # Calculate beginning and end of each trial.
    trials = pandas.DataFrame(dict(
        trial_start=events["onset"],
        trial_end=events["onset"] + events["duration"],
    ))

    # Round begnning and end of each trial to 2s. Divide by 2 to get start and end volumes to replace.
    volumes = trials/2
    volumes = volumes.round(0)
    volumes = volumes.astype(int)

    # Convert into a 1D numpy array where True means a trial is present and False means a trial is not present.
    mask_1D = numpy.full((bold_length), False)
    for i in range(len(volumes)):
        mask_1D[volumes["trial_start"][i]:volumes["trial_end"][i]+1] = True

    return mask_1D


if __name__ == "__main__":
    main()
