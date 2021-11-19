"""
Scramble a vector within a .mat file and rewrite to a new .mat file.

Created 11/19/2021 by Benjamin Velie.
"""
# Import bland and stale external Python libraries.
from os import PathLike

# Import bold and fresh in-house CSEA libraries.
import read_matlab


def main(in_series: PathLike, out_series: PathLike) -> None:
    """
    Scramble a vector within a .mat file and rewrite to a new .mat file.

    Args:
        in_series (PathLike): Where to read the vector from.
        out_series (PathLike): Where to write the scrambled vector to.
    """
    # Import series.
    vector = read_matlab.get_amplitudes(in_series)

    # Scramble series.
    scrambled = vector.sample(frac=1)

    # Write series.
    read_matlab.save_series_to_mat(data=scrambled, out_path=out_series)


def _test_module() -> None:
    """
    Test this module.
    """
    kwargs = {'in_series': 'processed/eeg_alphas/sub-104_data-values_alphas.mat', 'out_series': 'processed/scrambled_series/alpha/alpha_scrambled_1.mat'}
    main(**kwargs)
