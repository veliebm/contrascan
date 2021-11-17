"""
Prepare an image for being a figure.

NOT SUITABLE FOR PIPELINE USE. SUITABLE FOR INTERACTIVE USE.

Created 11/16/2021 by Benjamin Velie.
veliebm@ufl.edu
"""
# Import external code..
from os import PathLike
from pathlib import Path

# Import custom CSEA code.
import apply_matrix
import apply_mask
import copy_image


def prep_TTN27(
    in_image: PathLike,
    in_matrix: PathLike = "data/misc/kastner_cortex_masks/MNI152_T1_1mm_TTN27_mat.aff12.1D",
    mask_to: PathLike = "data/misc/kastner_cortex_masks/TT_N27_2.5mm.nii",
    get_underlay_from: PathLike = "data/misc/kastner_cortex_masks/TT_N27.nii.gz",
) -> None:
    """
    Align a 2.5mm image from Kastner to the TT_N27 brain then mask it. Copy an underlay to the dir.

    Args:
        in_image (PathLike): Image to process.
        in_matrix (PathLike, optional): Path to 1D file to apply. Defaults to "data/misc/kastner_cortex_masks/MNI152_T1_1mm_TTN27_mat.aff12.1D".
        mask_to (PathLike, optional): Path to 2.5mm mask to apply. Defaults to "data/misc/kastner_cortex_masks/TT_N27_2.5mm.nii".
        get_underlay_from (PathLike, optional): Where to copy an underlay from. Defaults to "data/misc/kastner_cortex_masks/TT_N27.nii.gz".
    """
    prefix = get_prefix(in_image)

    copy_underlay_to = prefix.parent / Path(get_underlay_from).name
    copy_image.main(get_underlay_from, copy_underlay_to)

    aligned_image = prefix.parent / (prefix.name + "_aligned.nii.gz")
    apply_matrix.main(in_image=in_image, in_matrix=in_matrix, out_prefix=aligned_image)

    masked_image = prefix.parent / (prefix.name + "_aligned_masked.nii.gz")
    apply_mask.main(in_image=aligned_image, in_mask=mask_to, out_prefix=masked_image)


def prep_MNI(
    in_image: PathLike,
    mask_to: PathLike = "data/misc/kastner_cortex_masks/MNI152_T1_2.5mm_full_mask.nii.gz",
    get_underlay_from: PathLike = "data/misc/kastner_cortex_masks/MNI152_T1_1mm.nii.gz",
) -> None:
    """
    Mask an image to the brain used by the Kastner atlas. Copy the atlas to the dir if it isn't already present.

    Args:
        in_image (PathLike): Image to process.
        mask_to (PathLike, optional): Path to 2.5mm mask to apply. Defaults to "data/misc/kastner_cortex_masks/MNI152_T1_2.5mm_full_mask.nii.gz".
        get_underlay_from (PathLike, optional): Where to copy an underlay from. Defaults to "data/misc/kastner_cortex_masks/MNI152_T1_1mm.nii.gz".
    """
    prefix = get_prefix(in_image)

    copy_underlay_to = prefix.parent / Path(get_underlay_from).name
    copy_image.main(get_underlay_from, copy_underlay_to)

    masked_image = prefix.parent / (prefix.name + "_aligned_masked.nii.gz")
    apply_mask.main(in_image=in_image, in_mask=mask_to, out_prefix=masked_image)


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
