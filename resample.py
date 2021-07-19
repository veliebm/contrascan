#!/usr/bin/env python3
"""
Resample an image to another image.

Created 7/19/2021 by Benjamin Velie.
veliebm@ufl.edu
"""
# Import external modules and libraries.
from os import PathLike
import subprocess

def main(from_image: PathLike, to_image: PathLike, to_prefix: PathLike) -> None:
    """
    Aligns our IRFs to the MNI space used by the Kastner cortex mask.
    """
    command = f"""
        3dresample
        -master {to_image}
        -input {from_image}
        -prefix {to_prefix}
    """.split()
    subprocess.run(command)
