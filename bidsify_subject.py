#!/usr/bin/env python3
"""
Convert a single subject into BIDS format.

Created on 6/11/2021 by Benjamin Velie.
veliebm@gmail.com
"""
# Import external libraries and modules.
from os import PathLike
from shutil import copy
from typing import Dict, List
from pathlib import Path
import json
import pandas
import nibabel

# Import custom libraries and modules.
from vmrk import Vmrk
from dat import Dat

def main(sources: Dict, targets: Dict) -> Dict:
    """
    Converts a single subject into BIDS format.
    """
    copy_path(sources["eeg"], targets["eeg"])
    copy_path(sources["vmrk"], targets["vmrk"])
    copy_path(sources["vhdr"], targets["vhdr"])
    copy_path(sources["func"], targets["func"])
    copy_path(sources["anat"], targets["anat"])
    copy_path(sources["dat"], targets["dat"])

    write_func_tsv(targets["vmrk"], targets["dat"], targets["func"], targets["events"])
    write_func_json(targets["func"], targets["func_json"])

def copy_path(source: PathLike, destination: PathLike) -> None:
    """
    Copy a path to its new home.
    """
    destination = Path(destination).resolve()
    if not destination.exists():
        destination.parent.mkdir(parents=True, exist_ok=True)
        copy(source, destination)
        print(f"Copied {source}  ->  {destination}")

def write_func_tsv(path_to_vmrk: PathLike, path_to_dat: PathLike, path_to_func: PathLike, out_path: PathLike) -> None:
    """
    Given a complete dataframe of files, writes appropriate tsvs for func files.

    Parameters
    ----------
    file_dataframe : DataFrame
        DataFrame created by organize_files() containing metadata about each file
    """
    vmrk_data = Vmrk(path_to_vmrk)
    dat_data = Dat(path_to_dat)

    onsets = vmrk_data.onsets
    durations = dat_data.durations

    # Prep a dataframe to write to .tsv.
    tsv_tuples = [ ("onset", "duration", "trial_type") ]
    for i, onset in enumerate(onsets):
        duration = durations[i]
        tsv_tuples.append( (onset, duration, "gabor") )
    tsv_dataframe = pandas.DataFrame(tsv_tuples)

    # Write the .tsv
    tsv_dataframe.to_csv(out_path, sep="\t", header=False, index=False)

def write_func_json(path_to_func: PathLike, out_path: PathLike) -> None:
    """
    Writes appropriate jsons for functional files.

    Parameters
    ----------
    file_dataframe : DataFrame
        DataFrame created by organize_files() containing metadata about each file
    """
    # Acquire values we'll need in the json.
    repetition_time = 2
    task = "contrascan"
    volume_count = get_volume_count(path_to_func)
    slice_timings = _calculate_slice_timings(repetition_time, volume_count)

    # Prepare a dict to write to json.
    json_dict = {
        "RepetitionTime": repetition_time,
        "TaskName": task,
        "SliceTiming": slice_timings
    }

    # Write the json.
    with open(out_path, "w") as out_file:
        json.dump(json_dict, out_file, indent="\t")

def _calculate_slice_timings(repetition_time: float, volume_count: int) -> List[float]:
    """
    Returns the slice timings in a list.

    Slices must be interleaved in the + direction.

    Parameters
    ----------
    repetition_time : float
        Repetition time of each scan in the functional image.
    volume_count : int
        Number of volumes in the functional image.

    Returns
    -------
    list
        List of slice timings.
    """
    # Generate a slice order of length volume_count in interleaved order.
    slice_order = list(range(0, volume_count, 2)) + list(range(1, volume_count, 2))

    # Calculate the slice timing list from the slice order.
    slice_timings = [slice / volume_count * repetition_time for slice in slice_order]
    return slice_timings

def get_volume_count(path_to_func: PathLike) -> int:
    """
    Returns the number of volumes in a func image.
    """
    image = nibabel.load(path_to_func)
    volume_count = image.header["slice_end"] - image.header["slice_start"] + 1

    return volume_count
