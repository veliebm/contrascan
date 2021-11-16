#!/usr/bin/env python3
"""
Do stuff involving resampling.

Created 7/19/2021 by Benjamin Velie.
veliebm@ufl.edu
"""
# Import external modules and libraries.
from os import PathLike
import subprocess
from pathlib import Path

def main(from_image: PathLike, to_image: PathLike, to_prefix: PathLike) -> None:
    """
    Resample an image to another image.

    Args:
        from_image (PathLike): Path to image you want to resample.
        to_image (PathLike): Path to image you want to resample to.
        to_prefix (PathLike): Where to write the resampled image to.
    """
    to_dir = Path(to_prefix).parent
    to_dir.mkdir(exist_ok=True, parents=True)

    command = f"""
        3dresample
        -master {to_image}
        -input {from_image}
        -prefix {to_prefix}
    """.split()
    subprocess.run(command)


def main2(from_image: PathLike, to_prefix: PathLike, dx: float, dy: float, dz: float) -> None:
    """
    Resample an image to a new resolution.

    Args:
        from_image (PathLike): Path to image you want to resample.
        to_prefix (PathLike): Where to write the resampled image to.
        dx (float): x resolution.
        dy (float): y resolution.
        dz (float): z resolution.
    """
    to_dir = Path(to_prefix).parent
    to_dir.mkdir(exist_ok=True, parents=True)

    command = f"""
        3dresample
        -dxyz {dx} {dy} {dz}
        -input {from_image}
        -prefix {to_prefix}
    """.split()
    subprocess.run(command)
