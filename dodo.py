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


# This example task executes a single analysis script for each subject, giving
# the subject as a command line parameter to the script.
#def task_example_step():
#    """Step 00: An example analysis step that is executed for each subject."""
#    # Run the example script for each subject in a sub-task.
#    for subject in SUBJECTS:
#        yield dict(
#            # This task should come after `task_check`
#            task_dep=['check'],
#
#            # A name for the sub-task: set to the name of the subject
#            name=subject,
#
#            # If any of these files change, the script needs to be re-run. Make
#            # sure that the script itself is part of this list!
#            file_dep=[fname.input(subject=subject), '00_example_step.py'],
#
#            # The files produced by the script
#            targets=[fname.output(subject=subject)],
#
#            # How the script needs to be called. Here we indicate it should
#            # have one command line parameter: the name of the subject.
#            actions=['python 00_example_step.py %s' % subject],
#        )
#
#
## Here is another example task that averages across subjects.
#def task_example_summary():
#    """Step 01: Average across subjects."""
#    return dict(
#        task_dep=['example_step'],  # This task should come after `task_example_step`
#        file_dep=[fname.output(subject=s) for s in SUBJECTS] + ['01_grand_average.py'],
#        targets=[fname.grand_average],
#        actions=['python 01_grand_average.py'],
#    )
