#!/usr/bin/env python3
"""
Apply a mask to an image.

Created 8/9/2021 by Benjamin Velie.
"""

from os import PathLike
import subprocess
from pathlib import Path


def main(in_image: PathLike, in_mask: PathLike, out_prefix: PathLike) -> None:
    """
    Apply a mask to an image.
    """
    Path(out_prefix).parent.mkdir(exist_ok=True, parents=True)
    args = f"""
        3dcalc
        -float
        -a {in_image}
        -b {in_mask}
        -expr a*b
        -prefix {out_prefix}
        """.split()
    print(args)
    subprocess.run(args)
