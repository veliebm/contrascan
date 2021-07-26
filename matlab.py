#!/usr/bin/env python3
"""
Run a MatLab script.

Created 6/29/21 by Ben Velie.
"""
from os import PathLike
import subprocess


def main(path_to_script: PathLike) -> None:
    """
    Run a MatLab script.
    """
    run_matlab_script(path_to_script)


def run_matlab_script(path_to_script: PathLike) -> None:
    """
    Run a matlab script.
    """
    run_matlab_function(f"run('{path_to_script}')")


def run_matlab_function(matlab_function: str) -> None:
    """
    Run a matlab function. It'll be really slow, but that's MatLab for ya :(
    """
    args = "matlab.exe -nosplash -nodesktop -minimize -r".split()
    args += [f"try, {matlab_function}, catch e, disp(getReport(e)), exit(1), end, exit(0);"]

    subprocess.run(args)
