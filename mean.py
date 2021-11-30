"""
Calculate the mean of some images.

Created 11/29/2021 by Benjamin Velie.
"""

import subprocess
from os import PathLike
from pathlib import Path


def main(images: PathLike, out_prefix: PathLike) -> None:
    """
    Calculate the mean of some images.

    Args:
        images (PathLike): Paths to images.
        out_prefix (PathLike): Where to write the mean.
    """
    make_parent_dir(out_prefix)

    args = f"3dMean -prefix {out_prefix}".split()
    args += images

    print(args)
    subprocess.run(args)


def make_parent_dir(path: PathLike) -> None:
    """
    Create the parent directory for a path.
    """
    Path(path).parent.mkdir(exist_ok=True, parents=True)
