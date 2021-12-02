"""
Concatenate fMRI images.

Created 12/2/2021 by Benjamin Velie.
veliebm@gmail.com
"""
from os import PathLike
from typing import List
import subprocess
from pathlib import Path


def main(images: List[PathLike], out_prefix: PathLike) -> None:
    """
    Concatenate fMRI images together.

    Goes in the order they appear in the list. Uses 3dTcat.

    Args:
        images (List[PathLike]): List of paths to images to concatenate.
        out_prefix (PathLike): AFNI prefix of the outfile.
    """
    make_parent_dir(out_prefix)
    command = ["3dTcat", *images, "-prefix", out_prefix]
    print(command)
    subprocess.run(command)


def make_parent_dir(path: PathLike) -> None:
    """
    Make parent directories for the target path.

    Args:
        path (PathLike): A path.
    """
    Path(path).parent.mkdir(exist_ok=True, parents=True)
