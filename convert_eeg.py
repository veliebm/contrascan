#!/usr/bin/env python3
"""
Convert a subject from BrainVision format to EEGLAB format.

Created 6/29/21 by Ben Velie.
"""
from os import PathLike
import subprocess
from pathlib import Path

def main(brainvision_dir: PathLike, brainvision_name: str, converted_path: PathLike, setname: str) -> None:
    """
    Convert a subject from BrainVision format to EEGLAB format.
    """
    out_dir = Path(converted_path).parent
    out_dir.mkdir(exist_ok=True, parents=True)

    matlab_function = f"convert_eeg('{brainvision_dir}', '{brainvision_name}', '{converted_path}', '{setname}'); exit;"
    run_matlab_function(matlab_function, debug=True)

def run_matlab_function(matlab_function: str, debug=False) -> None:
    """
    Run a matlab function. It'll be really slow, but that's MatLab for ya :(
    """
    command = "matlab.exe -nosplash -nodesktop -minimize -r".split()
    command.append(matlab_function)

    if debug:
        print(f"Running {command}")

    subprocess.run(command)
