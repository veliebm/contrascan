"""
t-test 2 groups of images against each other.

Created 7/22/21 by Benjamin Velie.
veliebm@ufl.edu
"""

from typing import List
from os import PathLike
import subprocess


def main(list_a: List[PathLike], list_b: List[PathLike], prefix: PathLike) -> None:
    """
    Runs 3dttest++ on 2 sets of images.
    """
    args = f"3dttest++ -paired -prefix {prefix}".split()
    args += ["-setA"] + [f"{path}" for path in list_a]
    args += ["-setB"] + [f"{path}" for path in list_b]

    subprocess.run(args)
