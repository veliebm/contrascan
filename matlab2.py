#!/usr/bin/env python3
"""
Write then run a MatLab script.

Created 9/14/21 by Ben Velie.
"""
from os import PathLike
import subprocess
import argparse
from pathlib import Path
import time


def main(script_contents: str, write_script_to: PathLike) -> None:
    """
    Run a MatLab script.
    """
    write_script(write_script_to, script_contents)
    run_matlab_script(write_script_to)


def write_script(path_to_script: PathLike, body: str) -> None:
    """
    Literally just dump a string to disk.
    """
    with open(path_to_script, "w") as io:
        io.write(body)


def run_matlab_script(path_to_script: PathLike) -> None:
    """
    Run a matlab script.
    """
    lock_file = f"{path_to_script}_LOCKFILE"
    run_matlab_function(f"run('{path_to_script}')", lock_file)


def run_matlab_function(matlab_function: str, lock_file: PathLike) -> None:
    """
    Run a matlab function.
    """
    # Make lock file.
    lock_file = Path(lock_file)
    with open(lock_file, "w") as io:
        pass

    # Run script.
    args = "matlab.exe -nosplash -nodesktop -minimize -r".split()
    args += [f"try, {matlab_function}, delete('{lock_file}'), catch e, disp(getReport(e)), pause, exit(1), end, exit(0);"]
    subprocess.run(args)

    # Wait for lock file to be deleted.
    while lock_file.exists():
        time.sleep(1)


if __name__ == "__main__":
    """
    This code only runs if this file is called as a script.
    """
    parser = argparse.ArgumentParser(description="Run a matlab script.")

    parser.add_argument("--script_contents", required=True, help="What we should write as the body of the script.")
    parser.add_argument("--write_script_to", required=True, help="Where to write the script to.")

    parsed_args = parser.parse_args()
    parsed_args_dict = vars(parsed_args)
    main(**parsed_args_dict)
