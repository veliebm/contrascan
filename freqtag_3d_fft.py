#!/usr/bin/env python3
"""
Run FFT for each individual trial.

Created 9/9/21 by Ben Velie.
"""
from os import PathLike
import subprocess
import argparse
from pathlib import Path
import time


def main(write_script_to: PathLike, **kwargs) -> None:
    """
    Write and run a MatLab script that calculates our freqtag parameters.
    """
    lock_file = Path(f"{write_script_to}_LOCKFILE")
    write_matlab_script(write_script_to, lock_file, **kwargs)
    run_matlab_script(write_script_to, lock_file)


def write_matlab_script(write_script_to, lock_file, **kwargs):
    """
    Write our MatLab script to disk.
    """
    Path(write_script_to).parent.mkdir(exist_ok=True, parents=True)   

    read_segmented_eeg_from = Path(kwargs["read_segmented_eeg_from"])

    with open(write_script_to, "w") as io:
        io.write(f"""\
%% This script runs a fast fourier transform on our data.
do_one()
delete('{lock_file}');

function do_one()
    %% Load dataset.
    dataset = load_dataset('{read_segmented_eeg_from.name}', '{read_segmented_eeg_from.parent}');

    %% Read input variables.
    load('{kwargs["read_stimulus_start_from"]}')
    load('{kwargs["read_stimulus_end_from"]}')
    load('{kwargs["read_sampling_rate_from"]}')

    %% Run function.
    [spec] = freqtag3D_FFT(dataset, stimulus_start:stimulus_end, sampling_rate);

    %% Save output variables.
    save('{kwargs["write_spec_to"]}', 'spec');
end

function [dataset] = load_dataset(file_name, directory)
    % Load a dataset.
    eeglab;
    EEG = pop_loadset('filename', file_name, 'filepath', directory);
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
    parser = argparse.ArgumentParser(description="Run an FFT on segmented EEG data.")

    parser.add_argument("--write_script_to", required=True)
    parser.add_argument("--write_spec_to", required=True)
    parser.add_argument("--read_segmented_eeg_from", required=True)
    parser.add_argument("--read_sampling_rate_from", required=True)
    parser.add_argument("--read_stimulus_start_from", required=True)
    parser.add_argument("--read_stimulus_end_from", required=True)

    parsed_args = parser.parse_args()
    parsed_args_dict = vars(parsed_args)
    main(**parsed_args_dict)
