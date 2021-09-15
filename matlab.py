#!/usr/bin/env python3
"""
Run a MatLab script.

Created 6/29/21 by Ben Velie.
"""
from os import PathLike
import subprocess
import argparse
from pathlib import Path
import time


def main(run_script_at: PathLike) -> None:
    """
    Run a MatLab script.

    By default, Python doesn't wait for the script to finish executing. So we need to make a lock file, and then
    delete the lock file within the matlab script.

    main() will not finish executing until the lock file has been deleted.
    """
    lock_file = Path(f"{run_script_at}_LOCKFILE")
    
    with open(lock_file, "w") as io:
        pass
    run_matlab_script(run_script_at)

    while lock_file.exists():
        time.sleep(1)


def run_matlab_script(path_to_script: PathLike) -> None:
    """
    Run a matlab script.
    """
    run_matlab_function(f"run('{path_to_script}')")


def run_matlab_function(matlab_function: str) -> None:
    """
    Run a matlab function.
    """
    args = "matlab.exe -nosplash -nodesktop -minimize -r".split()
    args += [f"try, {matlab_function}, catch e, disp(getReport(e)), keyboard, exit(1), end, exit(0);"]

    subprocess.run(args)


if __name__ == "__main__":
    """
    This code only runs if this file is called as a script.
    """
    parser = argparse.ArgumentParser(description="Run a matlab script.")

    parser.add_argument("--run_script_at", required=True, help="Which script to run")

    parsed_args = parser.parse_args()
    parsed_args_dict = vars(parsed_args)
    main(**parsed_args_dict)
