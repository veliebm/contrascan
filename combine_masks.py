#!/usr/bin/env python3
"""
Combine a list of masks together.

Created 8/9/2021 by Benjamin Velie.
"""

from os import PathLike
import subprocess
from pathlib import Path
from typing import List

ALPHABET = "abcdefghijklmnopqrstuvwxyz"


def main(in_masks: List[PathLike], out_prefix: PathLike) -> None:
    """
    Combine a list of masks together into one uber-mask.
    """
    Path(out_prefix).parent.mkdir(exist_ok=True, parents=True)

    variables = []
    for i, mask in enumerate(in_masks):
        letter = ALPHABET[i]
        variables += f"-{letter} {mask}".split()

    expression = _build_expression(len(in_masks))

    run_3dcalc(expression, variables, out_prefix)


def _build_expression(number_of_letters: int) -> str:
    """
    Returns the -expr part of 3dcalc.
    """
    
    letters = list(ALPHABET[:number_of_letters])
    expression_parts = "+".join(letters)
    expression = f"-expr step({expression_parts})".split()

    return expression


def run_3dcalc(expression: str, variables: List[str], out_prefix: PathLike) -> None:
    args = f"""
        3dcalc
        -float
        -prefix {out_prefix}
        """.split()
    args += expression
    args += variables
    print(args)
    subprocess.run(args)
