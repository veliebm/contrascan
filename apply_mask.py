#!/usr/bin/env python3
"""
Apply a mask to an image.

Created 8/9/2021 by Benjamin Velie.
"""

from os import PathLike
import subprocess
from pathlib import Path
from typing import List
import re


def main(in_image: PathLike, in_mask: PathLike, out_prefix: PathLike) -> None:
    """
    Apply a mask to an image.
    """
    Path(out_prefix).parent.mkdir(exist_ok=True, parents=True)
    args = f"""
        3dcalc
        -float
        -a {in_image}
        -b {in_mask}
        -expr a*step(b)
        -prefix {out_prefix}
        """.split()
    print(args)
    subprocess.run(args)


def mask_images(images: List[PathLike], mask: PathLike) -> None:
    """Apply a mask to each image within a list of images.

    Args:
        images (List[PathLike]): Images to mask.
        mask (PathLike): Path to mask to apply.
    """
    paths = {image: f"{get_prefix(image)}_masked.nii.gz" for image in images}
    for in_path, out_prefix in paths.items():
        main(in_path, mask, out_prefix)


def get_prefix(filename: str) -> str:
    """Returns the prefix of an AFNI file. (Everything before the final "+".)

    Args:
        filename (str): Path to an AFNI file.

    Returns:
        str: Prefix of the file. Includes parent dirs in path.
    """
    if is_afni(filename):
        return "".join(filename.split("+")[:-1])
    else:
        path = Path(filename)
        stem = path.stem
        return path.parent / stem


def is_afni(filename: PathLike) -> bool:
    """Returns true if the input path is an AFNI path.

    Args:
        filename (PathLike): Path to a file.

    Returns:
        bool: True if path is an AFNI path.
    """
    is_afni = False
    stem = Path(filename).stem
    if re.search(pattern=r"(\+tlrc|\+orig|\+acpc)$", string=stem):
        is_afni = True

    return is_afni
