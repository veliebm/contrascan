#!/usr/bin/env python3
"""
Calculate the parameters for the freqtag steps.

Created 9/7/21 by Ben Velie.
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
    with open(write_script_to, "w") as io:
        io.write(f"""\
%% This script calculates the parameters we'll need for the freqtag pipeline.
do_one()
delete('{lock_file}');

function do_one()
    %% Oz sensor number.
    oz_id = 20;

    %% Get epoch duration.
    num_time_points = 2650;
    baseline = 4/5*500;
    stimulus_start = baseline + 500;
    stimulus_end = 2483;
    epoch_duration = stimulus_end - stimulus_start
    
    
    %% Define all frequencies contained in the data, see section 3.1 for a guide on how to delimit the frequency axis.
    frequency_resolution = 1/(epoch_duration*2);
    sampling_rate = 500;
    num_frequencies_all = sampling_rate / 2;
    faxisall = 0:frequency_resolution:num_frequencies_all;

    
    %% Eliminate unnecessary frequencies.
    faxis = faxisall(12:2:100);

    %% Save variables.
    save('{kwargs["write_faxis_to"]}', 'faxis');
    save('{kwargs["write_faxisall_to"]}', 'faxisall');
    save('{kwargs["write_stimulus_start_to"]}', 'stimulus_start');
    save('{kwargs["write_stimulus_end_to"]}', 'stimulus_end');
    save('{kwargs["write_epoch_duration_to"]}', 'epoch_duration');
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
    parser = argparse.ArgumentParser(description="Calculate the basic parameters of the freqtag pipeline.")

    parser.add_argument("--write_script_to", required=True, help="Where to write the generated MatLab script.")
    parser.add_argument("--write_faxis_to", required=True, help="Where to write faxis.")
    parser.add_argument("--write_faxisall_to", required=True, help="Where to write faxis_all.")
    parser.add_argument("--write_stimulus_start_to", required=True, help="Where to write the stimulus start time.")
    parser.add_argument("--write_stimulus_end_to", required=True, help="Where to write the stimulus end time.")
    parser.add_argument("--write_epoch_duration_to", required=True, help="Where to write epoch duration.")

    parsed_args = parser.parse_args()
    parsed_args_dict = vars(parsed_args)
    main(**parsed_args_dict)
