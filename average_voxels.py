#!/usr/bin/env python3
"""
Average together the voxels within an ROI. Works for either a 3d or 4d image.

3dmaskave help: https://afni.nimh.nih.gov/pub/dist/doc/htmldoc/programs/3dmaskave_sphx.html#ahelp-3dmaskave
"""

from os import PathLike
from pathlib import Path
import subprocess


def main(in_image: PathLike, out_file: PathLike, mask: PathLike) -> None:
    """
    Average together the voxels within an ROI.
    """
    args = f"""
    3dmaskave
    -mask {mask}
    {in_image}
    """.split()

    run_while_redirecting_stdout(args, out_file)


def run_while_redirecting_stdout(args, out_path: PathLike) -> None:
    """
    Run a program. Redirect stdout to a file.
    """
    Path(out_path).parent.mkdir(exist_ok=True, parents=True)
    with open(out_path, "w") as io:
        subprocess.run(args, stdout=io)
