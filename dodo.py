#!/usr/bin/env python3
"""
Do-it script to execute the entire pipeline using the doit tool:
http://pydoit.org

All the filenames are defined in config.py
"""
# Import external modules and libraries.
from pathlib import Path
from typing import Dict

# Import internal modules and libraries.
from config import fname, SUBJECTS, n_jobs

# Import tasks
import create_bids_root
import bidsify_subject
import afniproc
import convert_eeg
import trim_func_images
import align
import resample
import smooth
import scale

# Configuration for the "doit" tool.
DOIT_CONFIG = dict(
    # While running scripts, output everything the script is printing to the
    # screen.
    verbosity=2,

    # Tell doit to use all processors I've said this machine has.
    num_process=n_jobs,

    # When the user executes "doit list", list the tasks in the order they are
    # defined in this file, instead of alphabetically.
    sort='definition',
)


def task_check() -> Dict:
    """Check the system dependencies."""
    return dict(
        file_dep=['check_system.py'],
        targets=[fname.system_check],
        actions=['python check_system.py']
    )


def task_create_bids_root() -> Dict:
    """
    Create the root of our bids dataset. We'll finish BIDSifiying the data when we add our individual subjects to the dataset.
    """
    return dict(
        actions=[(create_bids_root.main, [fname.bids_dir])],
        task_dep=["check"],
        file_dep=["create_bids_root.py"],
        targets=[fname.bids_description],
    )


def task_bidsify_subject() -> Dict:
    """
    Convert all subjects to BIDS format.

    We'll need this for fMRIPrep. Also, when we submit our dataset to the NIH, they'll want it in BIDS format.
    """
    for subject in SUBJECTS:
        sources = dict(
            eeg=fname.raw_eeg(subject=subject),
            vhdr=fname.raw_vhdr(subject=subject),
            vmrk=fname.raw_vmrk(subject=subject),
            dat=fname.raw_dat(subject=subject),
            func=fname.raw_func(subject=subject),
            anat=fname.raw_anat(subject=subject),
        )
        targets = dict(
            anat=fname.bids_anat(subject=subject),
            eeg=fname.bids_eeg(subject=subject),
            vmrk=fname.bids_vmrk(subject=subject),
            vhdr=fname.bids_vhdr(subject=subject),
            func_json=fname.bids_func_json(subject=subject),
            func=fname.bids_func(subject=subject),
            events=fname.bids_events(subject=subject),
            dat=fname.bids_dat(subject=subject),
        )

        yield dict(
            name=subject,
            actions=[(bidsify_subject.main, [sources, targets])],
            file_dep=list(sources.values()),
            targets=list(targets.values()),
        )


def task_afniproc() -> Dict:
    """
    Run afni_proc.py to preprocess and deconvolve our subjects.

    This is the base of our pipeline. We'll use the outputs of afni_proc.py for our more advanced analyses.
    """
    for subject in SUBJECTS:
        sources = dict(
            vmrk=fname.bids_vmrk(subject=subject),
            func=fname.bids_func(subject=subject),
            anat=fname.bids_anat(subject=subject),
        )
        targets = dict(
            log=fname.afniproc_log(subject=subject),
            command=fname.afniproc_command(subject=subject),
            deconvolved=fname.afniproc_deconvolved(subject=subject),
            irf=fname.afniproc_irf(subject=subject),
            anat=fname.afniproc_anat(subject=subject),
        )

        kwargs = dict(
            vmrk_path=sources["vmrk"],
            func_path=sources["func"],
            anat_path=sources["anat"],
            out_dir=fname.afniproc_subject_dir(subject=subject),
            subject_id=subject,
            remove_first_trs=1,
        )

        yield dict(
            name=subject,
            actions=[(afniproc.main, [], kwargs)],
            file_dep=list(sources.values()),
            targets=list(targets.values()),
        )


def task_align_afniproc_irfs() -> Dict:
    """
    Align our afniproc IRFs to the space of the Kastner cortex masks so we may compare them with fMRIPrep's IRFs.

    Note that the aligned IRFs appear in the afniproc dir, not the alignment template dir.
    I'm aware that this is quirky, but it's convenient.

    Also, we don't separate each subject into a sub-task for a very good reason. It's more efficient to calculate the alignment
    for one subject, then apply the transformation to all the others.
    """
    sources = [fname.afniproc_irf(subject=subject) for subject in SUBJECTS]
    sources += [
        fname.atlas_template,
        fname.afniproc_template,
    ]

    targets = [fname.afniproc_aligned_irf(subject=subject) for subject in SUBJECTS]

    kwargs = dict(
        from_images=[fname.afniproc_irf(subject=subject) for subject in SUBJECTS],
        from_template=fname.afniproc_template,
        to_template=fname.atlas_template,
        to_dir=fname.afniproc_alignment_dir,
    )

    return dict(
        actions=[(align.main, [], kwargs)],
        file_dep=sources,
        targets=targets,
    )


def task_resample_afniproc_irfs() -> Dict:
    """
    Resample our afniproc IRFs to the space of the Kastner cortex masks so we may compare them with fMRIPrep's IRFs.
    """
    for subject in SUBJECTS:
        sources = [
            fname.atlas_template,
            fname.afniproc_aligned_irf(subject=subject)
        ]
        targets = [
            fname.afniproc_resampled_irf(subject=subject)
        ]

        kwargs = dict(
            from_image=fname.afniproc_aligned_irf(subject=subject),
            to_image=fname.atlas_template,
            to_prefix=get_prefix(fname.afniproc_resampled_irf(subject=subject)),
        )

        yield dict(
            name=subject,
            actions=[(resample.main, [], kwargs)],
            file_dep=sources,
            targets=targets,
        )


def task_smooth_fmriprep() -> Dict:
    """
    Smooth our fMRIPrep images.
    """
    for subject in SUBJECTS:
        sources = [
            fname.fmriprep_func(subject=subject)
        ]

        targets = [
            fname.fmriprep_smoothed(subject=subject)
        ]

        kwargs = dict(
            image=fname.fmriprep_func(subject=subject),
            prefix=get_prefix(fname.fmriprep_smoothed(subject=subject)),
        )

        yield dict(
            name=subject,
            actions=[(smooth.main, [], kwargs)],
            file_dep=sources,
            targets=targets,
        )


def task_scale_fmriprep() -> Dict:
    """
    Scale our smoothed fMRIPrep images.
    """
    for subject in SUBJECTS:
        sources = [
            fname.fmriprep_smoothed(subject=subject)
        ]
        targets = [
            fname.fmriprep_scaled(subject=subject)
        ]

        kwargs = dict(
            image=fname.fmriprep_smoothed(subject=subject),
            prefix=get_prefix(fname.fmriprep_scaled(subject=subject)),
        )

        yield dict(
            name=subject,
            actions=[(scale.main, [], kwargs)],
            file_dep=sources,
            targets=targets,
        )


def task_trim_func_images() -> Dict:
    """
    Truncate our functional images to begin when the first stimulus was presented.

    Necessary for correlating the functional images with our EEG data.
    We got the stimulus presentation times from the afniproc results.
    """
    for subject in SUBJECTS:
        sources = dict(
            onsets=fname.afniproc_onsets(subject=subject),
            func=fname.afniproc_func(subject=subject),
        )
        targets = dict(
            trimmed_func=fname.trimmed_func(subject=subject)
        )

        kwargs = dict(
            onsets_path=sources["onsets"],
            func_path=sources["func"],
            out_prefix="".join(targets["trimmed_func"].split("+")[:-1]),
        )

        yield dict(
            name=subject,
            actions=[(trim_func_images.main, [], kwargs)],
            file_dep=list(sources.values()),
            targets=list(targets.values()),
        )


def task_convert_eeg() -> Dict:
    """
    Convert BrainVision EEG files into EEGLAB EEG files, which are easier to edit.
    """
    for subject in SUBJECTS:
        sources = dict(
            brainvision_path=fname.brainvision_eeg(subject=subject),
        )
        targets = dict(
            converted_path=fname.converted_eeg(subject=subject),
        )
        kwargs = dict(
            brainvision_dir=fname.brainvision_dir,
            brainvision_name = Path(sources["brainvision_path"]).name,
            converted_path = targets["converted_path"],
            setname = f"sub-{subject}",
        )

        yield dict(
            name=subject,
            actions=[(convert_eeg.main, [], kwargs)],
            file_dep=list(sources.values()),
            targets=list(targets.values()),
        )


def _print_dict(dictionary: Dict, name: str=None) -> None:
    """
    Prints a dictionary and whether its paths exist.

    Useful for debugging.
    """
    if name:
        print(name)
    for key, value in dictionary.items():
        path = Path(value)
        exists = path.exists()
        exists_str = "exists" if exists else "doesn't exist"
        print(f"{key}  :  {value}  :  {exists_str}")
    print()


def get_prefix(filename: str) -> str:
    """
    Returns the prefix of an AFNI file. (Everything before the final "+".)
    """
    return "".join(filename.split("+")[:-1])


def get_view(filename: str) -> str:
    """
    Returns the view of an AFNI file. (Everything after the final "+".)
    """
    return filename.split("+")[-1]


def add_suffix(filename: str, suffix: str) -> str:
    """
    Adds a suffix to the end of an AFNI filename.
    """
    prefix = get_prefix(filename)
    view = get_view(filename)
    return prefix + suffix + "+" + view
