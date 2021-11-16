"""
Apply a transformation matrix produced by align_epi_anat.py.

Created 11/16/2021 by Ben Velie.
veliebm@gmail.com
"""
from os import PathLike
import subprocess


def main(in_image: PathLike, in_matrix: PathLike, out_prefix: PathLike) -> None:
    """
    Apply a 1D transformation matrix to an image.

    Args:
        in_image (PathLike): Path to image to transform.
        in_matrix (PathLike): Path to matrix to use.
        out_prefix (PathLike): Where to write transformed image.
    """
    command = f"""
        3dAllineate
        -cubic
        -1Dmatrix_apply {in_matrix}
        -prefix {out_prefix}
        {in_image}
    """.split()
    print(command)
    subprocess.run(command, check=True)
