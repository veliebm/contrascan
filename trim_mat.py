#!/usr/bin/env python3
"""
Module to trim a vector within a mat file by a specified number of units.

Created by Ben Velie, 11/9/2021.
veliebm@gmail.com
"""

# Import standard libraries.
from os import PathLike

# Import custom CSEA libraries.
from read_matlab import get_amplitudes, save_series_to_mat


def main(from_mat: PathLike, to_trimmed_mat: PathLike, start_index: int = None, end_index: int = None) -> None:
    """
    Module to trim a vector within a mat file by a specified number of units.

    Args:
        from_mat (PathLike): Path to mat file containing a vector.
        to_trimmed_mat (PathLike): Where to write the trimmed mat to.
        start_index (int, optional): Index to trim from from the start of the vector. Defaults to None.
        end_index (int, optional): Index to trim from from the end of the vector. Defaults to None.

    Raises:
        RuntimeError: If neither start_index nor end_endix are specified.
    """
    if start_index is None and end_index is None:
        raise RuntimeError(f"You must specify either a start index or an end index.")

    vector = get_amplitudes(from_mat)
    if start_index is not None:
        vector = vector[start_index:]
    if end_index is not None:
        vector = vector[:end_index]

    save_series_to_mat(vector, "trimmed", to_trimmed_mat)


def _test_module() -> None:
    """
    Test this module.
    """
    kwargs = {'from_mat': './processed/canonical_bold_response/sub-125_canonical.mat', 'to_trimmed_mat': './processed/canonical_bold_response_trimmed/sub-125_canonical.mat', 'start_index': 4, 'end_index': -2}
    main(**kwargs)
