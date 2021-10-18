#!/usr/bin/env python3
"""
Use 3dmaskave to get the average signal within a 3d volume inside a specific ROI.

Created 10/18/2021 by Benjamin Velie
veliebm@gmail.com
"""

import subprocess
from os import PathLike
from pathlib import Path


def main(from_image: PathLike, to_text_file: PathLike, from_mask: PathLike) -> None:
    """
    Get the average signal within a 3d volume inside a specific ROI.
    """
    Path(to_text_file).parent.mkdir(exist_ok=True, parents=True)

    command = f"3dmaskave -mask '{from_mask}' '{from_image}' > '{to_text_file}'"
    print(command)
    subprocess.run(command, shell=True)
