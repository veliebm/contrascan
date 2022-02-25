"""
Remove trials from a BOLD image.

Created 9/30/2022 By Benjamin Velie.
"""

# Import external Python stuff.
from os import PathLike
from pathlib import Path
from typing import Dict
import pandas
import scipy.io

# Import custom CSEA stuff.
from read_matlab import get_amplitudes


def main(better_sliding_window: PathLike, sliding_sliding_window: PathLike, events: PathLike, improved_sliding_sliding_window: PathLike):
    """
    Entrypoint of the script.
    """
    # Read onsets and durations from events file. EEG data begins 1s before first onset and ends 10s after final onset.
    # BOLD data has been trimmed such that the first onset is located SOMEWHERE inside 1st BOLD volume.
    # So we should still adjust the onsets to begin from that point.
    # We should then remove all TRs from our BOLD image that contain trials.
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
    improved_amplitudes = overwrite_sliding_sliding_window(sliding_sliding_window_series, volume_replacements)


def overwrite_sliding_sliding_window(sliding_sliding_window_series: pandas.Series, volume_replacements: Dict) -> pandas.Series:
    """
    Fill trial blocks in sliding sliding window with sliding window results.
    """
    new_series = sliding_sliding_window_series.copy(deep=True)
    for volume, amplitude in volume_replacements.items():
        new_series[volume] = amplitude
    
    return new_series


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


if __name__ == "__main__":
    main()
