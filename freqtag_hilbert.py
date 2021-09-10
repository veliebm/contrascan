#!/usr/bin/env python3
"""
Do a Hilbert transform!

Created 9/10/21 by Ben Velie.
"""
from os import PathLike
import subprocess
import sys
from pathlib import Path
import time


def main(write_script_to: PathLike, **kwargs) -> None:
    """
    THE NEXUS.
    """
    lock_file = Path(f"{write_script_to}_LOCKFILE")
    write_matlab_script(write_script_to, lock_file, **kwargs)
    run_matlab_script(write_script_to, lock_file)


def write_matlab_script(write_script_to, lock_file, **kwargs):
    """
    Write our MatLab script to disk.
    """
    Path(write_script_to).parent.mkdir(exist_ok=True, parents=True)   

    with open(write_script_to, "w") as io:
        io.write(f"""\
%% This script runs a Hilbert transform.
do_one()
delete('{lock_file}');

function do_one()
    %% Load dataset.
    dataset = load_dataset('{kwargs["segmented_eeg"]}');

    %% Read input variables.
    load('{kwargs["stimulus_start"]}')
    load('{kwargs["stimulus_end"]}')
    load('{kwargs["sampling_rate"]}')

    %% Run function.
    [powermat, phasemat, complexmat] = freqtag_HILB(mean(dataset(:, stimulus_start:stimulus_end, :), 3), {kwargs["frequency"]}, 8, {kwargs["electrode"]}, 0, sampling_rate);

    %% Save output variables.
    save('{kwargs["powermat"]}', 'powermat');
    save('{kwargs["phasemat"]}', 'phasemat');
    save('{kwargs["complexmat"]}', 'complexmat');
end

function [dataset] = load_dataset(path)
    % Load a dataset.
    [parent_dir, stem, suffix] = fileparts(path);
    eeglab;
    EEG = pop_loadset('filename', strcat(stem, suffix), 'filepath', parent_dir);
    EEG = eeg_checkset(EEG);
    dataset = double(EEG.data);
end
""")


def run_matlab_script(path_to_script: PathLike, lock_file: PathLike) -> None:
    """
    Run a matlab script.

    By default Python doesn't wait for the script to finish executing. We need to make a lock file and then
    delete the lock file within the matlab script.

    This function won't finish executing until the lock file has been deleted.
    """   
    with open(lock_file, "w") as io:
        pass
    run_matlab_function(f"run('{path_to_script}')")

    while lock_file.exists():
        time.sleep(1)


def run_matlab_function(matlab_function: str) -> None:
    """
    Run a matlab function. It'll be really slow, but that's MatLab for ya :(
    """
    args = "matlab.exe -nosplash -nodesktop -minimize -r".split()
    args += [f"try, {matlab_function}, catch e, disp(getReport(e)), pause, exit(1), end, exit(0);"]

    subprocess.run(args)


if __name__ == "__main__":
    """
    This code only runs if this file is called as a script.
    """
    kwargs = dict(arg.split("=") for arg in sys.argv[1:])
    main(**kwargs)
