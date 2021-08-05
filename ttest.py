"""
t-test some images.

Created 7/28/21 by Benjamin Velie.
veliebm@ufl.edu
"""

from typing import List
from os import PathLike
import subprocess


def main(images: List[PathLike], prefix: PathLike) -> None:
    """
    Runs 3dttest++ on 2 sets of images.
    """
    args = f"3dttest++ -zskip 100% -prefix {prefix}".split()
    args += ["-setA"] + [f"{path}" for path in images]

    subprocess.run(args)
