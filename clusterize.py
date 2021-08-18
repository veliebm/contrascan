#!/usr/bin/env python3
"""
Make an image demarcating all clusters within an image.

Created 8/18/2021 by Benjamin Velie.
"""

from os import PathLike
from typing import Union
import subprocess
from pathlib import Path


def main(in_image: PathLike, out_prefix: PathLike, out_summary: PathLike, subbrick: Union[str, int]) -> None:
    """
    Make an image demarcating all clusters within an image.
    """
    args = f"""
    3dClusterize
    -inset {in_image}
    -ithr {subbrick}
    -pref_map {out_prefix}
    -bisided 0 0
    -NN 1
    """.split()

    run_while_redirecting_stdout(args, out_summary)


def run_while_redirecting_stdout(args, out_path: PathLike) -> None:
    """
    Run a program. Redirect stdout to a file.
    """
    Path(out_path).parent.mkdir(exist_ok=True, parents=True)
    with open(out_path, "w") as io:
        subprocess.run(args, stdout=io)
