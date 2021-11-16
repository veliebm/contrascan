"""
Prepare an image for being a figure.

NOT SUITABLE FOR PIPELINE USE. SUITABLE FOR INTERACTIVE USE.

Created 11/16/2021 by Benjamin Velie.
veliebm@ufl.edu
"""
from os import PathLike
import apply_matrix
import apply_mask
from pathlib import Path


def main(
    in_image: PathLike,
    in_matrix: PathLike = "data/misc/kastner_cortex_masks/MNI152_T1_1mm_TTN27_mat.aff12.1D",
    mask_to: PathLike = "data/misc/kastner_cortex_masks/TT_N27_2.5mm.nii",
) -> None:
    """
    Align a 2.5mm image from Kastner to the TT_N27 brain then mask it.

    Args:
        in_image (PathLike): Image to process.
        apply_matrix (PathLike, optional): Path to 1D file to apply. Defaults to "data/misc/kastner_cortex_masks/MNI152_T1_1mm_TTN27_mat.aff12.1D".
        mask_to (PathLike, optional): Path to 2.5mm mask to apply. Defaults to "data/misc/kastner_cortex_masks/TT_N27_2.5mm.nii".
    """
    prefix = get_prefix(in_image)

    aligned_image = prefix.parent / (prefix.name + "_aligned.nii.gz")
    apply_matrix.main(in_image=in_image, in_matrix=in_matrix, out_prefix=aligned_image)

    masked_image = prefix.parent / (prefix.name + "_aligned_masked.nii.gz")
    apply_mask.main(in_image=aligned_image, in_mask=mask_to, out_prefix=masked_image)


def get_prefix(path: PathLike) -> str:
    """
    Returns the prefix of an AFNI file. (Everything before the final "+".)

    If not AFNI then return stem.
    """
    pathlib_path = Path(path)
    parent_dir = pathlib_path.parent
    prefix = "".join(pathlib_path.name.split(".")[:-1])

    if "+" in prefix:
        prefix = "".join(prefix.split("+")[:-1])

    return parent_dir / prefix
