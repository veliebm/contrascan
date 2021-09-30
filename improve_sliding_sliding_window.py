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

    # Calculate which volumes to replace.
    volumes_to_replace = calculate_volumes(events_table, sliding_window_series)

    # Downsample onsets into 2s intervals coded based on which trial they belong to. Intervals with no trial contain null or 0.

    # Fill trial blocks with sliding window results.

    # Read the sliding sliding window amplitudes.

    # Splice the downsampled better sliding window results over the sliding sliding window amplitudes.

    # Write the improved amplitudes to a .m file.


def calculate_volumes(events_table: pandas.DataFrame, sliding_window_series: pandas.Series) -> pandas.DataFrame:
    """
    Downsample onsets into 2s intervals coded based on which trial they belong to. Intervals with no trial contain null or 0.
    """
    # Calculate beginning and end of each trial.

    # Round begnning and end of each trial to 2s.

    # Divide by 2 to get start and end volumes to replace.

    # Return DataFrame containing a column of start volumes, end volumes, and sliding window values for those intervals.
    breakpoint()


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
        "./processed/eeg_sliding_sliding_window/sub-122_frequency-12_amplitudes.fftamp.mat",
        "./processed/bids/sub-122/func/sub-122_task-contrascan_events.tsv",
        "./processed/eeg_sliding_sliding_window_improved/sub-122_improved_amplitudes.m",
    )
