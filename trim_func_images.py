#!/usr/bin/env python3
"""
Tools to truncate functional images to remove all data recorded before the first stimulus.
To use, call main().

Created on 7/6/2021 by Ben Velie.
"""

# Standard Python modules.
from os import PathLike
from pathlib import Path
import subprocess

def main(onsets_path: PathLike, func_path: PathLike, out_prefix: PathLike) -> None:
    """
    Truncate our images!
    """
    new_start_volume = _get_volume(_get_start_time(onsets_path))

    truncate_volumes(func_path, new_start_volume, out_prefix)


def _get_start_time(onsets_path: PathLike) -> float:
    """
    Returns the first onset time.
    """
    with open(onsets_path, "r") as f:
        onsets = f.readlines()
    
    return float(onsets[0])


def _get_volume(time_in_seconds: float) -> int:
    """
    Returns the volume containing the indicated time point.
    """
    return int((time_in_seconds - time_in_seconds % 2) / 2)


def truncate_volumes(func_path: PathLike, new_start_volume: int, out_prefix: PathLike) -> Path:
    """
    Removes volumes from the beginning of a 4D image. Returns the path to the .HEAD outfile.
    """
    command = [
        "3dTcat",
        f"{func_path}[{new_start_volume}..$]",
        "-prefix", out_prefix,
    ]
    subprocess.run(command)
