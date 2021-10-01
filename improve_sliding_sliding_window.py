"""
Doit task to overwrite the trial parts of the sliding sliding window with more accurate values.

We have good trial estimates from the sliding window. The sliding sliding window uses poor trial estimates.
If we overwrite the bad trial estimates with our good trial estimates, we expect that the sliding sliding
window analysis will produce better results.

Created 9/30/2021 By Benjamin Velie, veliebm@ufl.edu
"""

# Import external Python stuff.
from os import PathLike
import pandas

# Import custom CSEA stuff.
from read_matlab import get_amplitudes


def main(better_sliding_window: PathLike, sliding_sliding_window: PathLike, events: PathLike, improved_sliding_sliding_window: PathLike):
    """
    Overwrite the trial parts of the sliding sliding window with more accurate values.

    We have good trial estimates from the sliding window. The sliding sliding window uses poor trial estimates.
    If we overwrite the bad trial estimates with our good trial estimates, we expect that the sliding sliding
    window analysis will produce better results.
    """
    # Read onsets and durations from events file. Adjust onsets to begin 1s before the first onset.
    events_table = get_events(events)

    # Read sliding window data.
    sliding_window_series = get_amplitudes(better_sliding_window)

    # Calculate the start times, end times, and slid win values for each trial.
    trials_table = calculate_trials(events_table, sliding_window_series)

    # Get exactly which volumes to replace with which values.
    volume_replacements = calculate_volumes_to_replace(trials_table)

    # Read the sliding sliding window amplitudes.
    sliding_sliding_window_series = get_amplitudes(sliding_sliding_window)

    # Fill trial blocks in sliding sliding window with sliding window results.


    # Write the improved amplitudes to a .m file.


def calculate_volumes_to_replace(trials_table: pandas.DataFrame) -> pandas.Series:
    """
    Get exactly which volumes to replace with which values.
    """
    volumes = dict()
    for i in trials_table.index:
        start = trials_table["trial_start"][i]
        end = trials_table["trial_end"][i]
        amplitude = trials_table["amplitude"][i]
        volumes_to_replace = list(range(start, end))
        for volume in volumes_to_replace:
            volumes[volume] = amplitude

    return volumes


def calculate_trials(events_table: pandas.DataFrame, sliding_window_series: pandas.Series) -> pandas.DataFrame:
    """
    Calculate the start times, end times, and slid win values for each trial.
    """
    # Calculate beginning and end of each trial.
    trials = pandas.DataFrame(dict(
        trial_start = events_table["onset"],
        trial_end = events_table["onset"] + events_table["duration"],
    ))
    
    # Round begnning and end of each trial to 2s. Divide by 2 to get start and end volumes to replace.
    volumes = trials/2
    volumes = volumes.round(0)
    volumes = volumes.astype(int)

    # Return DataFrame containing a column of start volumes, end volumes, and sliding window values for those intervals.
    volumes["amplitude"] = sliding_window_series
    return volumes


def get_events(path_to_events: PathLike) -> pandas.DataFrame:
    """
    Read an events file then adjust onsets to 1s before the first stimulus. Returns a dataframe of events.
    """
    events_table = _read_events(path_to_events)

    first_onset = events_table["onset"][0]
    adjust_times_by = first_onset - 1

    adjusted_events_table = events_table.copy(deep=True)
    adjusted_events_table["onset"] = events_table["onset"] - adjust_times_by

    return adjusted_events_table


def _read_events(path_to_events: PathLike) -> pandas.DataFrame:
    """
    Return the events file as a DataFrame.
    """
    events_table = pandas.read_csv(path_to_events, sep="\t")
    return events_table


def _test_task():
    """
    Runs our code with dummy data.
    """
    main(
        "./processed/freqtag_better_sliding_window/sub-122_variable-trialpow_channel-20.mat",
        "./processed/eeg_sliding_sliding_window/sub-104_frequency-12_oz_amplitudes.mat",
        "./processed/bids/sub-122/func/sub-122_task-contrascan_events.tsv",
        "./processed/eeg_sliding_sliding_window_improved/sub-122_improved_amplitudes.m",
    )
