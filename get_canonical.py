"""
Calculate the canonical BOLD response we expect given our contrast function and the stimulus onsets.

Created 11/8/2021 by Benjamin Velie.
"""
from os import PathLike
from pathlib import Path
import textwrap
from typing import List
import matlab2


def main(from_onsets: PathLike, to_trs: PathLike, write_script_to: PathLike) -> None:
    """Calculate the canonical BOLD response we expect given our contrast function and the stimulus onsets.

    Args:
        from_onsets (PathLike): Path to the afniproc onsets file.
        to_trs (PathLike): Where to write our BOLD response to. Should be a .mat file.
        write_script_to (PathLike): Where to write the temporary MatLab script to.
    """
    onsets = extract_onsets(from_onsets)
    formatted_onsets = _format_onsets_to_matlab_code(onsets)

    execute_script(formatted_onsets, to_trs, write_script_to)


def execute_script(formatted_onsets: str, to_trs: PathLike, write_script_to: PathLike) -> None:
    """Write and execute a lil matlab script to calc the canonical BOLD response.

    Args:
        formatted_onsets (str): A list of onsets formatted as a MatLab matrix.
        to_trs (PathLike): Where to write the TRs to.
        write_script_to (PathLike): Where to write the temp script to.
    """
    make_parent_dir(to_trs)
    script = textwrap.dedent(f"""\
        onsets = {formatted_onsets};
        canonical = contrascan_canonical(onsets);
        save('{to_trs}', 'canonical');
        """)
    matlab2.main(script_contents=script, write_script_to=write_script_to)


def _format_onsets_to_matlab_code(onsets: List[float]) -> str:
    """Formats a list of onsets into a MatLab matrix.

    Args:
        onsets (List[float]): List of onsets.

    Returns:
        str: Matlab matrix.
    """
    onsets_as_strings = [str(onset) for onset in onsets]
    matrix = f"[{' '.join(onsets_as_strings)}]"

    return matrix


def extract_onsets(onsets_path: PathLike) -> List[float]:
    """
    Extract onsets from a path.
    """
    with open(onsets_path, "r") as f:
        onsets = f.readlines()

    return [float(onset) for onset in onsets]


def make_parent_dir(path: PathLike) -> None:
    """Make parent directories for the target path.

    Args:
        path (PathLike): A path.
    """
    Path(path).parent.mkdir(exist_ok=True, parents=True)


def _test_module() -> None:
    """
    Test this module.
    """
    kwargs = {'from_onsets': './processed/afniproc/sub-104/onsets.tsv', 'to_trs': './processed/canonical_bold_response/sub-104_canonical.mat', 'write_script_to': './TEMP_sub104_eeg_canonical.m'}
    main(**kwargs)
