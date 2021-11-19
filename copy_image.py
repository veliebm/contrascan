"""
Wrapper for 3dcopy. Copy a dataset!

Created 11/16/2021 by Ben Velie.
veliebm@gmail.com
"""
from os import PathLike
import subprocess


def main(in_prefix: PathLike, out_prefix: PathLike) -> None:
    """
    Wrapper for 3dcopy. Copy an image to a new place.

    Args:
        in_prefix (PathLike): Old location of image.
        out_prefix (PathLike): New location of image.
    """
    command = f"""
        3dcopy
        {in_prefix}
        {out_prefix}
    """.split()
    print(command)
    subprocess.run(command, check=False)
