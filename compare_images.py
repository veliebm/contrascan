"""
Subtract an image from another image.

Created 11/18/2021 by Benjamin Velie.
"""
from os import PathLike
from pathlib import Path
import subprocess


def main(minuend: PathLike, subtrahend: PathLike, out_prefix: PathLike) -> None:
    """
    Subtract an image from another image.
    """
    make_parent_directory(out_prefix)
    args = f"""
        3dcalc
        -float
        -a {minuend}
        -b {subtrahend}
        -expr a-b
        -prefix {out_prefix}
        """.split()
    print(args)
    subprocess.run(args)


def make_parent_directory(path: PathLike) -> None:
    """
    Create the parent directory for a path.

    Args:
        path (PathLike): Path for which to create a parent dir.
    """
    Path(path).parent.mkdir(exist_ok=True, parents=True)
