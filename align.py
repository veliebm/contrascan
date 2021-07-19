#!/usr/bin/env python3
"""
Align 2 datasets together.
"""
import subprocess
from os import PathLike
from pathlib import Path
from typing import List

def main(from_images: List[PathLike], from_template: PathLike, to_template: PathLike, to_dir: PathLike) -> None:
    """
    Aligns 2 datasets together.
    """
    from_images = [Path(from_image).resolve() for from_image in from_images]
    from_template = Path(from_template).resolve()
    to_template = Path(to_template).resolve()
    to_dir = Path(to_dir).resolve()

    to_dir.mkdir(parents=True, exist_ok=True)

    align_using_templates(from_images, from_template, to_template, to_dir)

def align_using_templates(from_images: List[PathLike], from_template: PathLike, to_template: PathLike, cwd: PathLike) -> Path:
    """
    Using two template images, align an image from one template space to another.
    """
    print(f"Aligning images to {to_template}")
    command = f"""
        align_epi_anat.py
        -dset1 {from_template}
        -dset2 {to_template}
        -dset1to2
        -dset1_strip None
        -dset2_strip None
        -suffix _aligned
        -child_dset1
        """.split()
    command += from_images
    print(command)
    subprocess.run(command, cwd=cwd)
