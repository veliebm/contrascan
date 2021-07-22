#!/usr/bin/env python3
"""
Deconvolve an fMRIPrep result.
"""
# Import external modules and libraries.
from os import PathLike
from pathlib import Path
import subprocess
from typing import List
import pandas

def main(
    regressors_tsv: PathLike,
    regressors_dir: PathLike,
    deconvolved_prefix: PathLike,
    IRF_prefix: PathLike,
    func_path: PathLike,
    events_tsv: PathLike, 
    events_dir: PathLike,
    ):
    """
    Deconvolve an fMRIPrep result.
    """
    regressors = "csf csf_derivative1 csf_power2 csf_derivative1_power2 white_matter white_matter_derivative1 white_matter_derivative1_power2 white_matter_power2 csf_wm trans_x trans_x_derivative1 trans_x_power2 trans_x_derivative1_power2 trans_y trans_y_derivative1 trans_y_derivative1_power2 trans_y_power2 trans_z trans_z_derivative1 trans_z_derivative1_power2 trans_z_power2 rot_x rot_x_derivative1 rot_x_power2 rot_x_derivative1_power2 rot_y rot_y_derivative1 rot_y_power2 rot_y_derivative1_power2 rot_z rot_z_derivative1 rot_z_power2 rot_z_derivative1_power2".split()
    onsets_path = Path(events_dir) / "onset.txt"

    prepare_to_deconvolve_fmriprep(regressors_tsv, events_tsv, regressors_dir, events_dir)
    deconvolve(regressors, regressors_dir, func_path, deconvolved_prefix, IRF_prefix, onsets_path)


def prepare_to_deconvolve_fmriprep(regressors_tsv: PathLike, events_tsv: PathLike, regressors_dir: PathLike, events_dir: PathLike):
    """
    Does what it says on the tin.
    """
    Path(regressors_dir).mkdir(exist_ok=True, parents=True)
    Path(events_dir).mkdir(exist_ok=True, parents=True)

    _split_columns_into_text_files(events_tsv, events_dir)
    _split_columns_into_text_files(regressors_tsv, regressors_dir)


def deconvolve(regressors: List[str], regressors_dir: PathLike, func_path: PathLike, deconvolved_prefix: PathLike, IRF_prefix: PathLike, onsets_path: PathLike):
    """
    Deconvolve a functional image.

    3dDeconvolve info: https://afni.nimh.nih.gov/pub/dist/doc/htmldoc/programs/3dDeconvolve_sphx.html#ahelp-3ddeconvolve
    """
    regressors_dir = Path(regressors_dir)

    # Total amount of regressors to include in the analysis.
    amount_of_regressors = 1 + len(regressors)

    # Create list of arguments.
    args = f"""
        3dDeconvolve
        -input {func_path}
        -GOFORIT 4
        -polort A
        -fout
        -bucket {deconvolved_prefix}
        -num_stimts {amount_of_regressors}
        -stim_times 1 {onsets_path} CSPLINzero(0,18,10)
        -stim_label 1 all
        -iresp 1 {IRF_prefix}
    """.split()

    # Add individual stim files to the string.
    for i, regressor_name in enumerate(regressors):
        stim_number = i + 2
        stim_file_info = f"-stim_file {stim_number} {regressors_dir/regressor_name}.txt -stim_base {stim_number}"
        stim_label_info = f"-stim_label {stim_number} {regressor_name}"
        args += stim_file_info.split() + stim_label_info.split()

    subprocess.run(args)


def _split_columns_into_text_files(tsv_path: PathLike, output_dir: PathLike):
    """
    Converts a tsv file into a collection of text files.

    Each column name becomes the name of a text file. Each value in that column is then
    placed into the text file. Don't worry - this won't hurt your .tsv file, which will lay
    happily in its original location.

    Parameters
    ----------
    tsv_path : str or Path
        Path to the .tsv file to break up.
    output_dir : str or Path
        Directory to write columns of the .tsv file to.
    """
    # Prepare our paths.
    tsv_path = Path(tsv_path).absolute()
    output_dir = Path(output_dir).absolute()
    output_dir.mkdir(exist_ok=True, parents=True)
    print(f"Storing the columns of {tsv_path.name} as text files in directory {output_dir}")

    # Read the .tsv file into a DataFrame and fill n/a values with zero.
    tsv_info = pandas.read_table(
        tsv_path,
        sep="\t",
        na_values="n/a"
    ).fillna(value=0)

    # Write each column of the dataframe as a text file.
    for column_name in tsv_info:
        column_path = output_dir / f"{column_name}.txt"
        tsv_info[column_name].to_csv(column_path, sep=' ', index=False, header=False)
