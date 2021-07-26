"""
Write a JSON file.

Created 7/26/21 by Benjamin Velie.
veliebm@ufl.edu
"""

from os import PathLike
from typing import Any
import json
from pathlib import Path


def main(out_path: PathLike, data: Any):
    """
    Write a json file.
    """
    Path(out_path).parent.mkdir(exist_ok=True, parents=True)
    with open(out_path, "w") as f:
        json.dump(data, f, indent="\t")
