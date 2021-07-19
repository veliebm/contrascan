#!/usr/bin/env python3
"""
Smooth an image.

Created 7/19/2021 by Benjamin Velie.
veliebm@ufl.edu
"""
# Import external modules and libraries.
from os import PathLike
import subprocess
from pathlib import Path

def main(image: PathLike, prefix: PathLike) -> None:
    """
    Smooth an image.

    3dmerge info: https://afni.nimh.nih.gov/pub/dist/doc/htmldoc/programs/3dmerge_sphx.html#ahelp-3dmerge
    """
    dir = Path(prefix).parent
    dir.mkdir(exist_ok=True, parents=True)
    
    args = f"""
        3dmerge
        -1blur_fwhm 4.0
        -doall
        -prefix {prefix}
        {image}
    """.split()
    subprocess.run(args)
