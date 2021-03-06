#!/usr/bin/env python3
"""
Using two template images, align images from one template space to another.
"""
import subprocess
from os import PathLike
from pathlib import Path
from typing import List


def main(from_images: List[PathLike], from_template: PathLike, to_template: PathLike, to_dir: PathLike, suffix: str = "_aligned") -> None:
    """
    Using two template images, align images from one template space to another.

    Args:
        from_images (List[PathLike]): Paths to images you want to align.
        from_template (PathLike): Template the input images are aligned to.
        to_template (PathLike): What to align the images to.
        to_dir (PathLike): Where you want the results to be written to.
        suffix (str, optional): Suffix to append to filenames. Defaults to "_aligned".
    """
    from_images = [Path(from_image).resolve() for from_image in from_images]
    from_template = Path(from_template).resolve()
    to_template = Path(to_template).resolve()
    to_dir = Path(to_dir).resolve()

    to_dir.mkdir(parents=True, exist_ok=True)

    _align_using_templates(from_images, from_template, to_template, to_dir, suffix)


def _align_using_templates(from_images: List[PathLike], from_template: PathLike, to_template: PathLike, cwd: PathLike, suffix: str) -> None:
    """
    Using two template images, align images from one template space to another. Requires that input paths be full paths.

    Args:
        from_images (List[PathLike]): Paths to images you want to align.
        from_template (PathLike): Template the input images are aligned to.
        to_template (PathLike): What to align the images to.
        cwd (PathLike): Where you want the results to be written to.
        suffix (str): Suffix to append to filenames.
    """
    print(f"Aligning images to {to_template}")
    command = f"""
        align_epi_anat.py
        -dset1 {from_template}
        -dset2 {to_template}
        -dset1to2
        -dset1_strip None
        -dset2_strip None
        -suffix {suffix}
        -child_dset1
        """.split()
    command += from_images
    print(command)
    subprocess.run(command, cwd=cwd)
