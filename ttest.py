"""
t-test some images.

Created 7/28/21 by Benjamin Velie.
veliebm@ufl.edu
"""

from typing import List
from os import PathLike
import subprocess
from pathlib import Path


def main(images: List[PathLike], prefix: PathLike) -> None:
    """
    Runs 3dttest++ on 2 sets of images.
    """
    make_parent_directory(prefix)

    args = f"3dttest++ -zskip 100% -prefix {prefix}".split()
    args += ["-setA"] + [f"{path}" for path in images]

    subprocess.run(args, check=True)


def make_parent_directory(path: PathLike) -> None:
    """
    Create the parent directory for a path.
    """
    Path(path).parent.mkdir(exist_ok=True, parents=True)
