#!/usr/bin/env python3
"""
Scale an image.

Created 7/19/2021 by Benjamin Velie.
veliebm@ufl.edu
"""
from os import PathLike
import subprocess
from pathlib import Path


def main(image: PathLike, prefix: PathLike) -> None:
    """
    Scale an image.
    """
    out_dir = Path(prefix).parent
    out_dir.mkdir(exist_ok=True, parents=True)

    temp_prefix = prefix + "_temp"
    temp_image = temp_prefix + "+tlrc"

    tstat(image, temp_prefix)
    calc_percent_change(prefix, image, temp_image)


def tstat(image: PathLike, prefix: PathLike) -> None:
    """
    Get the mean of each voxel in a functional dataset.

    3dTstat info: https://afni.nimh.nih.gov/pub/dist/doc/htmldoc/programs/3dTstat_sphx.html#ahelp-3dtstat
    """
    args = f"""
        3dTstat
        -prefix {prefix}
        {image}
    """.split()
    subprocess.run(args)


def calc_percent_change(prefix: PathLike, image: PathLike, meaned_image: PathLike) -> None:
    """
    For each voxel in a smoothed func image, calculate the voxel as percent of the mean.

    3dcalc info: https://afni.nimh.nih.gov/pub/dist/doc/htmldoc/programs/3dcalc_sphx.html#ahelp-3dcalc
    """
    args = f"""
        3dcalc
        -float
        -a {image}
        -b {meaned_image}
        -expr ((a-b)/b)*100
        -prefix {prefix}
    """.split()
    subprocess.run(args)
