#!/usr/bin/env python3
"""
Do-it script to execute the entire pipeline using the doit tool:
http://pydoit.org

All the filenames are defined in config.py
"""
# Import external modules and libraries.
from os import PathLike
from pathlib import Path
from typing import Dict, Iterable, List, Tuple
import textwrap

# Import internal modules and libraries.
from config import EXPANDED_START_VOLUMES, FREQUENCIES, fname, SUBJECTS, COMPONENTS_TO_REMOVE, START_VOLUMES, PERMUTATIONS

# Import actions for tasks to use.
import create_bids_root
import bidsify_subject
import afniproc
import trim_func_images
import align
import resample
import write_json
import correlate_whole_brain
import ttest
import combine_masks
import apply_mask
import clusterize
import average_voxels
import correlate_regions
import matlab2
import improve_sliding_sliding_window
import correlate_timeseries
import correlate_tables
import func_get_trials
import maskave
import ttest_averages
import fisher_transform
import cohens
import get_canonical
import trim_mat
import compare_images
import scramble_series
import mean
import mean_means
import test_permutations
import cat_images
import calc_percentiles
import get_maxes_and_mins
import plot_distribution


# Configuration for the pydoit tool.
DOIT_CONFIG = dict(
    # While running scripts, output everything the script is printing to the screen.
    verbosity=2,

    # When the user executes "doit list", list the tasks in the order they are defined in this file, instead of alphabetically.
    sort='definition',
)


# General tasks.
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
            name=f"subject--{subject}",
            actions=[(bidsify_subject.main, [sources, targets])],
            file_dep=list(sources.values()),
            targets=list(targets.values()),
        )


# fMRI preprocessing.
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
            preprocessed_func=fname.afniproc_preprocessed_func(subject=subject)
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
            name=f"subject--{subject}",
            actions=[(afniproc.main, [], kwargs)],
            file_dep=list(sources.values()),
            targets=list(targets.values()),
        )
def task_align_func_images() -> Dict:
    """
    Align our preprocessed functional images to the space of the Kastner cortex masks.

    Note that the aligned images appear in the afniproc dir, not the alignment template dir.
    I'm aware that this is quirky, but it's convenient.

    Also, we don't separate each subject into a sub-task for a very good reason. It's more efficient to calculate the alignment
    for one subject, then apply the transformation to all the others.
    """
    sources = [fname.afniproc_preprocessed_func(subject=subject) for subject in SUBJECTS]
    sources += [
        fname.atlas_template,
        fname.afniproc_template,
    ]

    targets = [fname.aligned_func(subject=subject) for subject in SUBJECTS]

    kwargs = dict(
        from_images=[fname.afniproc_preprocessed_func(subject=subject) for subject in SUBJECTS],
        from_template=fname.afniproc_template,
        to_template=fname.atlas_template,
        to_dir=fname.alignment_dir,
    )

    return dict(
        actions=[(align.main, [], kwargs)],
        file_dep=sources,
        targets=targets,
    )
def task_resample_template() -> Dict:
    """
    The Kastner template is WAY too fine-grain. Smooth smooth smooth!
    """
    sources = [fname.atlas_template]
    targets = [fname.resampled_template]

    kwargs = dict(
        from_image=fname.atlas_template,
        to_prefix=get_prefix(fname.resampled_template),
        dx="2.5",
        dy="2.5",
        dz="2.5",
    )

    return dict(
        actions=[(resample.main2, [], kwargs)],
        file_dep=sources,
        targets=targets,
    )
def task_resample_func_images() -> Dict:
    """
    Resample our afniproc func images to the space of the Kastner cortex masks.
    """
    for subject in SUBJECTS:
        sources = [
            fname.resampled_template,
            fname.aligned_func(subject=subject)
        ]
        targets = [
            fname.resampled_func(subject=subject)
        ]

        kwargs = dict(
            from_image=fname.aligned_func(subject=subject),
            to_image=fname.resampled_template,
            to_prefix=get_prefix(fname.resampled_func(subject=subject)),
        )

        yield dict(
            name=f"subject--{subject}",
            actions=[(resample.main, [], kwargs)],
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
        sources = [
            fname.afniproc_onsets(subject=subject),
            fname.resampled_func(subject=subject),
        ]
        targets = [
            fname.trimmed_func(subject=subject),
        ]

        kwargs = dict(
            onsets_path=fname.afniproc_onsets(subject=subject),
            func_path=fname.resampled_func(subject=subject),
            out_prefix=get_prefix(fname.trimmed_func(subject=subject)),
        )

        yield dict(
            name=f"subject--{subject}",
            actions=[(trim_func_images.main, [], kwargs)],
            file_dep=sources,
            targets=targets,
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
            fname.resampled_template,
            fname.afniproc_aligned_irf(subject=subject)
        ]
        targets = [
            fname.afniproc_resampled_irf(subject=subject)
        ]

        kwargs = dict(
            from_image=fname.afniproc_aligned_irf(subject=subject),
            to_image=fname.resampled_template,
            to_prefix=get_prefix(fname.afniproc_resampled_irf(subject=subject)),
        )

        yield dict(
            name=f"subject--{subject}",
            actions=[(resample.main, [], kwargs)],
            file_dep=sources,
            targets=targets,
        )
def task_trim_func_images_again() -> Dict:
    """
    Trim func images one last time so we can correlate the BOLD response with EEG signal.
    """
    for subject in SUBJECTS:

        sources = [fname.trimmed_func(subject=subject)]

        for volumes_to_remove in EXPANDED_START_VOLUMES:

            targets = [
                fname.final_func(subject=subject, start_volume=volumes_to_remove),
            ]

            kwargs = dict(
                new_start_volume=volumes_to_remove,
                func_path=fname.trimmed_func(subject=subject),
                out_prefix=get_prefix(fname.final_func(subject=subject, start_volume=volumes_to_remove)),
            )

            yield dict(
                name=f"sub--{subject}, startvolume--{volumes_to_remove}",
                actions=[(trim_func_images.main2, [], kwargs)],
                file_dep=sources,
                targets=targets,
            )
def task_get_trial_func_amplitudes() -> Dict:
    """
    Get the biggest amplitude of the BOLD signal for each trial.

    For each S  2 stimulus, get its TR. Then look at the next 4 TRs after that. Select the TR containing the largest value.
    Then subtract from it the signal of the TR containing the S  2.
    """
    def create_task(sources: Dict[str, PathLike], targets: Dict[str, PathLike], name: str) -> dict:
        """
        Allows this task to easily be generalizable.

        Parameters
        ----------
        name : str
            Name of task.
        sources : Dict
            Contains paths to inputs of the task.
        targets : Dict
            Contains paths to outputs of the task.
        """

        kwargs = {**sources, **targets}

        return dict(
            name=name,
            actions=[(func_get_trials.main, [], kwargs)],
            file_dep=list(sources.values()),
            targets=list(targets.values()),
        )

    for subject in SUBJECTS:
        yield create_task(
            sources=dict(
                onsets_path=fname.afniproc_onsets(subject=subject),
                func_path=fname.resampled_func(subject=subject),
            ),
            targets=dict(
                trials_path=fname.func_trial_amplitudes(subject=subject)
            ),
            name=f"subject--{subject}",
        )


# fMRI-only analysis.
def task_resample_deconvolutions() -> Dict:
    """
    Resample the deconvolution results so we can t-test them against each other.
    """
    for subject in SUBJECTS:
        sources = [
            fname.resampled_template,
            fname.afniproc_deconvolved(subject=subject),
        ]
        targets = [
            fname.afniproc_deconvolved_resampled(subject=subject),
        ]

        kwargs = dict(
            from_image=fname.afniproc_deconvolved(subject=subject),
            to_image=fname.resampled_template,
            to_prefix=get_prefix(fname.afniproc_deconvolved_resampled(subject=subject)),
        )

        yield dict(
            name=f"subject--{subject}",
            actions=[(resample.main, [], kwargs)],
            file_dep=sources,
            targets=targets,
        )
def task_ttest_deconvolutions() -> Dict:
    """
    Here we t-test our 3dDeconvolve results. Onward!
    """
    for subbrick in range(18):
        sources = [fname.afniproc_deconvolved_resampled(subject=subject) for subject in SUBJECTS]
        targets = [fname.afniproc_ttest_result(subbrick=subbrick)]

        kwargs = dict(
            images=[fname.afniproc_deconvolved_resampled(subject=subject) + f"[{subbrick}]" for subject in SUBJECTS],
            prefix=get_prefix(fname.afniproc_ttest_result(subbrick=subbrick)),
        )
        yield dict(
            name=f"subbrick--{subbrick}",
            actions=[(ttest.main, [], kwargs)],
            file_dep=sources,
            targets=targets,
        )
def task_get_IRF_mean() -> Dict:
    """
    Calculate the mean IRF.
    """
    def create_task(images: PathLike, out_prefix: PathLike, name: str) -> Dict:
        """
        Calculate the mean of some images.

        Args:
            images (PathLike): Path to images to mean together.
            out_prefix (PathLike): Where to write the outfile to.
            name (str): Name of the task.

        Returns:
            Dict: [description]
        """
        sources = dict(images=images)
        targets = dict(out_prefix=out_prefix)
        kwargs = {**sources, **targets}

        return dict(
            name=name,
            actions=[(mean.main, [], kwargs)],
            file_dep=list(sources.values())[0],
            targets=list(targets.values()),
        )

    yield create_task(
        images=[fname.afniproc_resampled_irf(subject=subject) for subject in SUBJECTS],
        out_prefix=fname.IRF_mean,
        name="IRF"
    )

    yield create_task(
        images=[fname.afniproc_deconvolved(subject=subject) for subject in SUBJECTS],
        out_prefix=fname.deconvolve_mean,
        name="deconvolve"
    )


# Getting masks.
def task_get_occipital_mask() -> Dict:
    """
    Combine 4 Kastner masks to get 1 occipital pole mask. Now that's a steal!

    Specifically, combines Kastner V2d and V3d L/R.
    """
    sources = []
    for hemisphere in "lr":
        for roi_number in [4, 6]:
            sources.append(fname.kastner_mask(roi_number=roi_number, hemisphere=hemisphere))

    targets = [fname.occipital_pole_mask]

    kwargs = dict(
        in_masks=sources,
        out_prefix=get_prefix(fname.occipital_pole_mask),
    )
    return dict(
        actions=[(combine_masks.main, [], kwargs)],
        file_dep=sources,
        targets=targets,
    )
def task_get_calcarine_mask() -> Dict:
    """
    Combine 2 Kastner masks to get 1 calcarine mask. Now that's a steal!

    Specifically, combines Kastner V1v L/R.
    """
    sources = []
    for hemisphere in "lr":
        sources.append(fname.kastner_mask(roi_number=1, hemisphere=hemisphere))

    targets = [fname.calcarine_mask]

    kwargs = dict(
        in_masks=sources,
        out_prefix=get_prefix(fname.calcarine_mask),
    )
    return dict(
        actions=[(combine_masks.main, [], kwargs)],
        file_dep=sources,
        targets=targets,
    )
def task_resample_masks() -> Dict:
    """
    Resample our masks so we may apply them to our correlation results.
    """
    masks = [
        dict(name="calcarine", in_path=fname.calcarine_mask, out_path=fname.resampled_mask(mask="calcarine")),
        dict(name="occipital", in_path=fname.occipital_pole_mask, out_path=fname.resampled_mask(mask="occipital")),
    ]

    for mask in masks:

        sources = [
            mask["in_path"],
            fname.resampled_template,
        ]
        targets = [mask["out_path"]]

        kwargs = dict(
            from_image=mask["in_path"],
            to_image=fname.resampled_template,
            to_prefix=get_prefix(mask["out_path"]),
        )

        yield dict(
            name=f"mask--{mask['name']}",
            actions=[(resample.main, [], kwargs)],
            file_dep=sources,
            targets=targets,
        )
def task_apply_masks_to_irfs() -> Dict:
    """
    FOR each subject, USING her IRF, USING sub-brick 3, FOR each ROI, APPLY ROI mask to IRF.

    THEORY
    ------
    For each subject, we want to compare the most positive blob within the occipital pole with the most negative blob within the calcarine. Get these blobs from the IRFs. Presumable, sub-brick 3 or sub-brick 4.

    What does it mean to be the most positive blob? Easy. Within the ROI, find all contiguous blobs. I didn't ask whether touching edges count, but I bet they do. Within each contiguous blob, sum the signal of each voxel together. Pick the blob with the largest sum. Or in the case of the calcarine, the lowest sum. Huzzah!

    Then, use each voxel in that blob in the original correlation test. Maybe make a mask for that blob, then apply the mask to the original correlation t-tests. That's a good idea!
    """
    masks = [
        dict(name="calcarine", in_path=fname.resampled_mask(mask="calcarine")),
        dict(name="occipital", in_path=fname.resampled_mask(mask="occipital")),
    ]
    for subject in SUBJECTS:
        for mask in masks:
            sources = [
                fname.afniproc_resampled_irf(subject=subject),
                mask["in_path"],

            ]
            targets = [fname.masked_irf(subject=subject, mask=mask["name"])]

            kwargs = dict(
                in_image=fname.afniproc_resampled_irf(subject=subject),
                in_mask=mask["in_path"],
                out_prefix=get_prefix(fname.masked_irf(subject=subject, mask=mask["name"])),
            )

            yield dict(
                name=f"sub--{subject}, mask--{mask['name']}",
                actions=[(apply_mask.main, [], kwargs)],
                file_dep=sources,
                targets=targets,
            )
def task_clusterize_irfs() -> Dict:
    """
    FOR each subject, FOR each masked IRF, GET all clusters.
    """
    regions = "calcarine occipital".split()
    for subject in SUBJECTS:
        for region in regions:

            # Get sources.
            masked_irf = fname.masked_irf(subject=subject, mask=region)
            sources = [masked_irf]

            # Get targets.
            clusters = fname.clusters(subject=subject, mask=region)
            clusters_summary = fname.clusters_summary(subject=subject, mask=region)
            targets = [clusters, clusters_summary]

            # Get args to run.
            kwargs = dict(
                in_image=masked_irf,
                out_prefix=get_prefix(clusters),
                out_summary=clusters_summary,
                subbrick=3,
            )

            yield dict(
                name=f"subject--{subject}, region--{region}",
                actions=[(clusterize.main, [], kwargs)],
                file_dep=sources,
                targets=targets,
            )
def task_make_micromasks() -> Dict:
    """
    Convert our clusters into bona fide masks! Wow!

    FOR each image of clusters, GET strongest cluster.
    """
    regions = "calcarine occipital".split()
    for subject in SUBJECTS:
        for region in regions:

            # Get sources.
            clusters_image = fname.clusters(subject=subject, mask=region)
            clusters_summary = fname.clusters_summary(subject=subject, mask=region)
            sources = [clusters_image, clusters_summary]

            # Get targets.
            micromask = fname.micromask(subject=subject, mask=region)
            targets = [micromask]

            # Get action.
            action = f"""
                python3 micromasks.py
                --clusters_image {clusters_image}
                --clusters_1d_file {clusters_summary}
                --out_prefix {get_prefix(micromask)}
                {"--get_negative_cluster" if region == "calcarine" else ""}
            """.split()

            yield dict(
                name=f"subject--{subject}, region--{region}",
                actions=[action],
                file_dep=sources,
                targets=targets,
            )
def task_apply_micromasks_to_trimmed_trimmed_funcs() -> Dict:
    """
    Does what it says on the tin.
    """
    regions = "calcarine occipital".split()
    for subject in SUBJECTS:
        for region in regions:
            for start_volume in START_VOLUMES:
                trimmed_trimmed_func = fname.final_func(subject=subject, start_volume=start_volume)
                micromask = fname.micromask(subject=subject, mask=region)
                micromasked_func = fname.micromasked_func(subject=subject, mask=region, start_volume=start_volume)

                sources = [
                    trimmed_trimmed_func,
                    micromask,
                ]
                targets = [micromasked_func]

                kwargs = dict(
                    in_image=trimmed_trimmed_func,
                    in_mask=micromask,
                    out_prefix=get_prefix(micromasked_func),
                )

                yield dict(
                    name=f"subject--{subject}, mask--{region}, start_volume--{start_volume}",
                    actions=[(apply_mask.main, [], kwargs)],
                    file_dep=sources,
                    targets=targets,
                )
def task_average_microregion_voxels() -> Dict:
    """
    Within each microregion of the trimmed trimmed funcs, average together all voxels of that region into their own timeseries.
    """
    regions = "calcarine occipital".split()
    for subject in SUBJECTS:
        for region in regions:
            for start_volume in START_VOLUMES:
                micromasked_func = fname.micromasked_func(subject=subject, mask=region, start_volume=start_volume)
                averages = fname.microregion_average(subject=subject, mask=region, start_volume=start_volume)

                sources = [micromasked_func]
                targets = [averages]

                kwargs = dict(
                    in_image=micromasked_func,
                    out_file=averages,
                    mask="SELF"
                )

                yield dict(
                    name=f"subject--{subject}, mask--{region}, start_volume--{start_volume}",
                    actions=[(average_voxels.main, [], kwargs)],
                    file_dep=sources,
                    targets=targets,
                )


# EEG preprocessing.
def task_convert_eeg() -> Dict:
    """
    Convert brainvision files to EEGLAB files.
    """
    Path(fname.converteeg_dir).mkdir(exist_ok=True, parents=True)

    for subject in SUBJECTS:

        # Get sources.
        sources = dict(brainvision=fname.brainvision_eeg(subject=subject))
        sources_list = list(sources.values())

        # Get targets.
        targets = dict(converted_eeg=fname.converted_eeg(subject=subject), script=fname.converteeg_script(subject=subject))
        targets_list = list(targets.values())

        # Make MatLab script.
        script = textwrap.dedent(f"""\
        % Convert a subject from BrainVision format to EEGLAB format.
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
        EEG = eeg_checkset( EEG );
        EEG = pop_loadbv('{Path(sources["brainvision"]).parent}', '{Path(sources["brainvision"]).name}');
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname','sub-{subject}','savenew','{targets["converted_eeg"]}','gui','off');
        """)

        # Make action to run script.
        action = f"python3 matlab2.py".split()
        action += ["--script_contents", script]
        action += ["--write_script_to", targets["script"]]

        # Go!
        yield dict(
            name=f"subject--{subject}",
            actions=[action],
            file_dep=sources_list,
            targets=targets_list,
        )
def task_prepare_to_preprocess_eeg() -> Dict:
    """
    Write a JSON file that will be read later by a MatLab script I wrote to preprocess EEG files.
    """
    sources = [fname.converted_eeg(subject=subject) for subject in SUBJECTS]
    targets = [fname.preprocesseeg_json]

    data = []
    for subject in SUBJECTS:
        data.append(dict(
            subject=subject,
            components_to_remove=COMPONENTS_TO_REMOVE[subject],
            in_file=Path(fname.converted_eeg(subject=subject)).name,
            in_directory=fname.converteeg_dir,
            out_directory=fname.preprocesseeg_dir,
        ))

    kwargs = dict(
        out_path=fname.preprocesseeg_json,
        data=data,
    )

    return dict(
        actions=[(write_json.main, [], kwargs)],
        file_dep=sources,
        targets=targets,
    )
def task_preprocess_eeg() -> Dict:
    """
    Preprocess our EEG files.

    Because MatLab is ~quirky~, we can't do multithreading on this one.
    We must do all subjects in one glorious blaze!
    Sigh.
    """
    sources = [fname.preprocesseeg_json]
    targets = [fname.preprocessed_eeg(subject=subject) for subject in SUBJECTS]

    path_to_script = "preprocess_eegs.m"
    action = f"""
        python3 matlab.py
        --run_script_at {path_to_script}
    """.split()

    return dict(
        actions=[action],
        file_dep=sources,
        targets=targets,
    )
def task_prepare_to_segment_eeg() -> Dict:
    """
    Write a JSON file that will be read later by a MatLab script I wrote to segment our preprocessed EEG files.
    """
    sources = [fname.preprocessed_eeg(subject=subject) for subject in SUBJECTS]
    targets = [fname.segmenteeg_json]

    data = []
    for subject in SUBJECTS:
        data.append(dict(
            in_name=Path(fname.preprocessed_eeg(subject=subject)).name,
            in_directory=fname.preprocesseeg_dir,
            out_directory=fname.segmenteeg_dir,
        ))

    kwargs = dict(
        out_path=fname.segmenteeg_json,
        data=data,
    )

    return dict(
        actions=[(write_json.main, [], kwargs)],
        file_dep=sources,
        targets=targets,
    )
def task_segment_eeg() -> Dict:
    """
    Segment our preprocessed EEG files.

    Because we are using MatLab, we can't do multithreading on this one.
    We must do all subjects in one glorious blaze!
    """
    sources = [fname.segmenteeg_json]
    targets = [fname.segmented_eeg(subject=subject) for subject in SUBJECTS]

    path_to_script = "segment_eegs.m"
    action = f"""
        python3 matlab.py
        --run_script_at {path_to_script}
    """.split()

    return dict(
        actions=[action],
        file_dep=sources,
        targets=targets,
    )
def task_prepare_to_trim_eeg() -> Dict:
    """
    Write a JSON file that will be read later by a MatLab script I wrote to trim our
    preprocessed EEG files to the time at which the fMRI turned on.
    """
    sources = [fname.preprocessed_eeg(subject=subject) for subject in SUBJECTS]
    targets = [fname.trimeeg_json]

    data = []
    for subject in SUBJECTS:
        data.append(dict(
            subject=subject,
            in_filename=Path(fname.preprocessed_eeg(subject=subject)).name,
            in_dir=fname.preprocesseeg_dir,
            out_dir=fname.trimeeg_dir,
        ))

    kwargs = dict(
        out_path=fname.trimeeg_json,
        data=data,
    )

    return dict(
        actions=[(write_json.main, [], kwargs)],
        file_dep=sources,
        targets=targets,
    )
def task_trim_eeg() -> Dict:
    """
    Make our EEG data start at first stimulus and end 10s after final stimulus.
    """
    sources = [fname.trimeeg_json]
    targets = [fname.trimmed_eeg(subject=subject) for subject in SUBJECTS]

    path_to_script = "trim_eegs.m"
    action = f"""
        python3 matlab.py
        --run_script_at {path_to_script}
    """.split()

    return dict(
        actions=[action],
        file_dep=sources,
        targets=targets,
    )
def task_eeg_get_flicker_frequencies() -> Dict:
    """
    Get the flicker frequency of each trial. Save as .mat files.
    """
    Path(fname.eeg_flicker_frequencies_dir).mkdir(exist_ok=True, parents=True)
    for subject in SUBJECTS:

        # Get sources.
        sources = dict(
            dat=fname.bids_dat(subject=subject)
        )
        sources_list = list(sources.values())

        # Get targets.
        targets = dict(
            frequencies=fname.eeg_flicker_frequencies(subject=subject, variable="frequencies"),
            durations=fname.eeg_flicker_frequencies(subject=subject, variable="durations"),
            write_script_to=fname.eeg_flicker_frequencies_script(subject=subject),
        )
        targets_list = list(targets.values())

        # Make the script to run.
        script = textwrap.dedent(f"""\
            %% Get the flicker frequency of each trial. Save as a .mat file.
            do_one()

            function do_one()
                %% Read input variables.
                durations = get_durations('{sources["dat"]}')

                %% Run functions.
                frequencies = [];
                for i = 1:numel(durations)
                    duration = durations(i)
                    frequency = 1./(duration/50)
                    frequencies = [frequencies, frequency]
                end

                %% Save output variables.
                save('{targets["frequencies"]}', 'frequencies');
                save('{targets["durations"]}', 'durations');
            end

            function [durations] = get_durations(dat_path)
                datmat = importdata(dat_path)
                durations = datmat(1:end, 5);
            end
            """)

        # Make action to run script.
        action = f"python3 matlab2.py".split()
        action += ["--script_contents", script]
        action += ["--write_script_to", targets["write_script_to"]]

        # Go!
        yield dict(
            name=f"subject--{subject}",
            actions=[action],
            file_dep=sources_list,
            targets=targets_list,
        )


# Alpha analysis.
def task_eeg_get_alphas() -> Dict:
    """
    Get the amount of alpha present in each TR.
    """
    Path(fname.eeg_alpha_dir).mkdir(exist_ok=True, parents=True)
    for subject in SUBJECTS:

        # Get sources.
        sources = dict(
            eeg=fname.trimmed_eeg(subject=subject)
        )
        sources_list = list(sources.values())

        # Get targets.
        targets = dict(
            values=fname.eeg_alpha(subject=subject, data="values"),
            SNRs=fname.eeg_alpha(subject=subject, data="SNRs"),
            average_power=fname.average_power(subject=subject),
            write_script_to=fname.eeg_alpha_script(subject=subject),
        )
        targets_list = list(targets.values())

        # Make the script to run.
        script = textwrap.dedent(f"""\
            %% This script calculates the alpha value for each TR.

            %% Load dataset.
            dataset = load_dataset('{sources["eeg"]}');

            %% Run functions.
            values = [];
            pows = [];
            SNRs = [];
            for i = 1000:1000:numel(dataset(1,:))
                sub_dataset = dataset(:,i-999:i);

                [pow, phase, freqs] = FFT_spectrum(sub_dataset, 500);

                alpha = mean(mean(pow([20 31 19 7 8 9 10 ], 18:26)));

                [SNRdb, SNRratio] = freqtag_simpleSNR(pow, [11:16 28:33]);
                SNR = mean(mean(SNRdb([20 31 19 7 8 9 10 ], 18:26)));

                values = [values; alpha];
                SNRs = [SNRs; SNR];
                pows = cat(3, pows, pow);
            end
            average_power =  mean(pows, 3);

            %% Save output variables.
            save('{targets["values"]}', 'values');
            save('{targets["SNRs"]}', 'SNRs');
            save('{targets["average_power"]}', 'average_power');

            function [dataset] = load_dataset(path)
                % Load a dataset.
                [parent_dir, stem, suffix] = fileparts(path);
                eeglab;
                EEG = pop_loadset('filename', strcat(stem, suffix), 'filepath', parent_dir);
                EEG = eeg_checkset(EEG);
                dataset = double(EEG.data);
            end
            """)

        # Make action to run script.
        action = f"python3 matlab2.py".split()
        action += ["--script_contents", script]
        action += ["--write_script_to", targets["write_script_to"]]

        # Go!
        yield dict(
            name=f"subject--{subject}",
            actions=[action],
            file_dep=sources_list,
            targets=targets_list,
        )
def task_eeg_get_trial_by_trial_alpha() -> Dict:
    """
    Get the amount of alpha present in each trial.
    """
    def create_task(sources: Dict[str, PathLike], targets: Dict[str, PathLike], name: str) -> dict:
        """
        Allows this task to easily be generalizable.

        Parameters
        ----------
        out_dir : PathLike
            Dir where all outputs will be written.
        name : str
            Name of task.
        sources : Dict
            Contains paths to inputs of the task.
                eeg : PathLike
                    The segmented EEG data for this subject.
        targets : Dict
            Contains paths to outputs of the task.
                values : PathLike
                    The amplitudes for the trials, stored in a vector in a .m file.
                SNRs : PathLike
                    The signal to noise ratio for our amplitudes.
                average_power : PathLike
                    Average power over all trials.
        """
        out_dir = Path(targets["values"]).parent
        out_dir.mkdir(exist_ok=True, parents=True)

        script = textwrap.dedent(f"""\
            %% This script calculates the alpha value for each trial.

            %% Load dataset.
            dataset = load_dataset('{sources["eeg"]}');

            %% Run functions.
            values = [];
            pows = [];
            raw_pows = [];
            baseline_pows = [];
            SNRs = [];
            for i = 1:numel(dataset(1,1,:))
                sub_dataset = dataset(:,:,i);

                [pow_baseline, phase_baseline, freqs_baseline] = FFT_spectrum3D_singtrial(sub_dataset, 1:400, 500);
                [pow_post, phase_post, freqs_post] = FFT_spectrum3D_singtrial(sub_dataset, 2083:2482, 500);
                pow_difference = pow_post - pow_baseline;

                amplitude = mean(mean(pow_difference([20 31 19 7 8 9 10 ], 8:11)));

                [SNRdb, SNRratio] = freqtag_simpleSNR(pow_difference, [5:7 13:15]);
                SNR = mean(mean(SNRdb([20 31 19 7 8 9 10 ], 8:11)));

                values = [values; amplitude];
                SNRs = [SNRs; SNR];
                pows = cat(3, pows, pow_difference);
                raw_pows = cat(3, raw_pows, pow_post);
                baseline_pows = cat(3, baseline_pows, pow_baseline);
            end
            averagepower =  mean(pows, 3);
            averagerawpower = mean(raw_pows, 3);
            averagebaselinepower = mean(pow_baseline, 3);

            %% Save output variables.
            save('{targets["values"]}', 'values');
            save('{targets["SNRs"]}', 'SNRs');
            save('{targets["averagepower"]}', 'averagepower');
            save('{targets["averagerawpower"]}', 'averagerawpower');
            save('{targets["averagebaselinepower"]}', 'averagebaselinepower');

            function [dataset] = load_dataset(path)
                % Load a dataset.
                [parent_dir, stem, suffix] = fileparts(path);
                eeglab;
                EEG = pop_loadset('filename', strcat(stem, suffix), 'filepath', parent_dir);
                EEG = eeg_checkset(EEG);
                dataset = double(EEG.data);
            end
            """)

        return dict(
            name=name,
            actions=[(matlab2.main, [], dict(script_contents=script, write_script_to=targets["write_script_to"]))],
            file_dep=list(sources.values()),
            targets=list(targets.values()),
        )

    for subject in SUBJECTS:
        yield create_task(
            sources=dict(
                eeg=fname.segmented_eeg(subject=subject)
            ),
            targets=dict(
                values=fname.eeg_trial_alpha(subject=subject, data="values"),
                SNRs=fname.eeg_trial_alpha(subject=subject, data="SNRs"),
                averagepower=fname.eeg_trial_alpha(subject=subject, data="averagepower"),
                averagerawpower=fname.eeg_trial_alpha(subject=subject, data="averagerawpower"),
                averagebaselinepower=fname.eeg_trial_alpha(subject=subject, data="averagebaselinepower"),
                write_script_to=fname.eeg_trial_alpha_script(subject=subject),
            ),
            name=f"subject--{subject}",
        )
def task_correlate_alpha_and_snr() -> Dict:
    """
    Correlate the signal and SNR of our alpha data.
    """
    def create_task(sources: Dict[str, PathLike], targets: Dict[str, PathLike], other_kwargs: Dict, name: str) -> dict:
        """
        Allows this task to easily be generalizable.

        Parameters
        ----------
        sources : Dict
            data_1 : PathLike
                Path to a 1xN or Nx1 MatLab .mat file containing a list of numbers.
            data_2 : PathLike
                Path to a 1xN or Nx1 MatLab .mat file containing a list of numbers.
        targets : Dict
            save_spearman_to : PathLike
                Where to output our correlation results to.
            save_scatter_to : PathLike
                Where to output our scatter plot to.
            save_table_to : PathLike
                Where to output a table containing both streams of data.
        other_kwargs : Dict
            data_1_name : str
                What to name the data of data_1.
            data_2_name : str
                What to name the data of data_2.
        """
        kwargs = {**sources, **targets, **other_kwargs}

        return dict(
            name=name,
            actions=[(correlate_timeseries.main, [], kwargs)],
            file_dep=list(sources.values()),
            targets=list(targets.values()),
        )

    for subject in SUBJECTS:
        yield create_task(
            sources=dict(
                data_1=fname.eeg_alpha(subject=subject, data="values"),
                data_2=fname.eeg_alpha(subject=subject, data="SNRs"),
            ),
            targets=dict(
                save_spearman_to=fname.eeg_correlation_alpha_snr_results(subject=subject),
                save_table_to=fname.eeg_correlation_alpha_snr_table(subject=subject),
                save_scatter_to=fname.eeg_correlation_alpha_snr_scatter(subject=subject),
            ),
            other_kwargs=dict(
                data_1_name="alpha amplitude",
                data_2_name="alpha SNR",
            ),
            name=f"sub-{subject}",
        )
def task_correlate_alpha_and_snr_across_subjects() -> Dict:
    """
    Correlate the signal and SNR of our alpha data across all subjects.
    """
    def create_task(sources: Dict[str, PathLike], targets: Dict[str, PathLike]) -> dict:
        """
        Allows this task to be easily generalizable.

        Parameters
        ----------
        sources : Dict
            load_tables_from : List[PathLike]
                Where to get each of our microROI+amplitude tables
        targets : Dict
            save_scatter_to : PathLike
                Where to save our scatter plot.
            save_table_to : PathLike
                Where to save the table of data we make from our tiny tables.
            save_spearman_to : PathLike
                Where to save our Spearman results.
        """
        kwargs = {**sources, **targets}

        return dict(
            actions=[(correlate_tables.main, [], kwargs)],
            file_dep=list(sources.values())[0],
            targets=list(targets.values()),
        )

    return create_task(
        sources=dict(
            load_tables_from=[fname.eeg_correlation_alpha_snr_table(subject=subject) for subject in SUBJECTS],
        ),
        targets=dict(
            save_scatter_to=fname.eeg_correlation_across_subjects_alpha_snr_scatter,
            save_table_to=fname.eeg_correlation_across_subjects_alpha_snr_table,
            save_spearman_to=fname.eeg_correlation_across_subjects_alpha_snr_results,
        ),
    )


# Sliding window analysis.
def task_eeg_sliding_sliding_window() -> Dict:
    """
    Run a sliding sliding window analysis on our trimmed, preprocessed EEG files.

    Recall that stead2singtrials.m starts 1s before the first "S  2" stimulus and ends 10s after the final "S  2" stimulus.
    """
    Path(fname.eeg_sliding_sliding_window_dir).mkdir(exist_ok=True, parents=True)
    for subject in SUBJECTS:
        for frequency in FREQUENCIES:

            # Get sources.
            sources = dict(
                eeg=fname.preprocessed_eeg(subject=subject),
            )
            sources_list = list(sources.values())

            # Get targets.
            targets = dict(
                fftamp=fname.eeg_sliding_sliding_window_amplitudes(subject=subject, frequency=frequency),
                SNR=fname.eeg_sliding_sliding_window_SNR(subject=subject, frequency=frequency),
                oz_SNR=fname.eeg_sliding_sliding_window_oz_SNR(subject=subject, frequency=frequency),
                oz_fftamp=fname.eeg_sliding_sliding_window_oz_amplitudes(subject=subject, frequency=frequency),
                write_script_to=fname.eeg_sliding_sliding_window_script(subject=subject, frequency=frequency),
            )
            targets_list = list(targets.values())
            amplitudes_path = Path(targets["fftamp"])
            amplitudes_stem = str(amplitudes_path.parent / amplitudes_path.name.split(".")[0])

            # Make the script to run.
            script = textwrap.dedent(f"""\
                % Gets a SLIDING sliding window average for a subject.

                %% Main script.

                eeglab;
                [in_filename, in_dir] = split_path('{sources["eeg"]}')
                [winmat, SNR] = stead2singtrialsCont(in_filename, in_dir, 0, 1:1000, 1:1000, {frequency}, 600, 500, 1, '{amplitudes_stem}');
                load('{targets["fftamp"]}')
                oz_fftamp = fftamp(20,:)'
                oz_SNR = SNR(20,:)'

                % Save variables.
                save('{targets["oz_fftamp"]}', 'oz_fftamp');
                save('{targets["SNR"]}', 'SNR');
                save('{targets["oz_SNR"]}', 'oz_SNR');
                """)

            # Make action to run script.
            action = f"python3 matlab2.py".split()
            action += ["--script_contents", script]
            action += ["--write_script_to", targets["write_script_to"]]

            # Go!
            yield dict(
                name=f"subject--{subject}, frequency--{frequency}",
                actions=[action],
                file_dep=sources_list,
                targets=targets_list,
            )
def task_improve_sliding_sliding_window() -> Dict:
    """
    Overwrite the trial parts of the sliding sliding window with more accurate values.

    We have good trial estimates from the sliding window. The sliding sliding window uses poor trial estimates.
    If we overwrite the bad trial estimates with our good trial estimates, we expect that the sliding sliding
    window analysis will produce better results.

    Recall that the sliding sliding window amplitudes start 1s before the first "S  2" stimulus and end 10s after the final "S  2" stimulus.
    Each amplitude is for a single 2s TR.

    Recall that the events file is trimmed to begin when the fMRI turns on. We should adjust the timings for the events such that they are relative to
    1s before the first "S  2" stimulus. Then we should downsample the events into 2s intervals.
    """
    def create_task(sources: Dict[str, PathLike], targets: Dict[str, PathLike], name: str) -> dict:
        """
        Allows this task to easily be generalizable.

        Parameters
        ----------
        name : str
            Name of the task.
        sources : dict
            Contains paths to source files for the task.
                better_sliding_window : PathLike
                    Path to the better sliding window results.
                sliding_sliding_window : PathLike
                    Path to the TR sliding window results.
                events : PathLike
                    Path to the BIDS events file containing the onset and duration of each trial.
        targets : dict
            Contains paths to target files for the task.
                improved_sliding_sliding_window : PathLike
                    Where to output improved sliding sliding window results.
        """
        kwargs = dict(**sources, **targets)

        return dict(
            name=name,
            actions=[(improve_sliding_sliding_window.main, [], kwargs)],
            file_dep=list(sources.values()),
            targets=list(targets.values()),
        )

    for subject in SUBJECTS:
        yield create_task(
            sources=dict(
                better_sliding_window=fname.freqtag_better_sliding_window_channel(subject=subject, channel=20, variable="trialpow"),
                sliding_sliding_window=fname.eeg_sliding_sliding_window_oz_amplitudes(subject=subject, frequency=FREQUENCIES[0]),
                events=fname.bids_events(subject=subject),
            ),
            targets=dict(
                improved_sliding_sliding_window=fname.eeg_sliding_sliding_window_improved(subject=subject, variable="amplitudes"),
            ),
            name=f"subject--{subject}, variable--amplitudes"
        )
        yield create_task(
            sources=dict(
                better_sliding_window=fname.freqtag_better_sliding_window_channel(subject=subject, channel=20, variable="trialSNR"),
                sliding_sliding_window=fname.eeg_sliding_sliding_window_oz_SNR(subject=subject, frequency=FREQUENCIES[0]),
                events=fname.bids_events(subject=subject),
            ),
            targets=dict(
                improved_sliding_sliding_window=fname.eeg_sliding_sliding_window_improved(subject=subject, variable="SNRs"),
            ),
            name=f"subject--{subject}, variable--SNRs"
        )


# Freqtag pipeline.
def task_freqtag_calculate_parameters() -> Dict:
    """
    Calculate the parameters of the freqtag pipeline.
    """
    Path(fname.freqtag_parameters_dir).mkdir(exist_ok=True, parents=True)

    # Get sources.
    sources = ["freqtag_calculate_parameters.py"]

    # Get targets.
    script = fname.freqtag_parameters_script
    faxis = fname.freqtag_faxis
    faxisall = fname.freqtag_faxisall
    stimulus_start = fname.freqtag_stimulus_start
    stimulus_end = fname.freqtag_stimulus_end
    epoch_duration = fname.freqtag_epoch_duration
    sampling_rate = fname.freqtag_sampling_rate
    targets = [
        script,
        faxis,
        faxisall,
        stimulus_start,
        stimulus_end,
        epoch_duration,
        sampling_rate
    ]

    # Get args.
    action = f"""\
        python3 freqtag_calculate_parameters.py
        --write_script_to {script}
        --write_faxis_to {faxis}
        --write_faxisall_to {faxisall}
        --write_stimulus_start_to {stimulus_start}
        --write_stimulus_end_to {stimulus_end}
        --write_epoch_duration_to {epoch_duration}
        --write_sampling_rate_to {sampling_rate}
    """.split()

    # Go!
    return dict(
        actions=[action],
        file_dep=sources,
        targets=targets,
    )
def task_freqtag_fft() -> Dict:
    """
    Make an average of our trials then do a Fast Fourier Transform on it.
    """
    Path(fname.freqtag_fft_dir).mkdir(exist_ok=True, parents=True)

    for subject in SUBJECTS:
        # Get sources.
        segmented_eeg = fname.segmented_eeg(subject=subject)
        sampling_rate = fname.freqtag_sampling_rate
        stimulus_start = fname.freqtag_stimulus_start
        stimulus_end = fname.freqtag_stimulus_end
        sources = [
            segmented_eeg,
            sampling_rate,
            stimulus_start,
            stimulus_end
        ]

        # Get targets.
        script = fname.freqtag_fft_script(subject=subject)
        pow = fname.freqtag_fft_pow(subject=subject)
        phase = fname.freqtag_fft_phase(subject=subject)
        freqs = fname.freqtag_fft_freqs(subject=subject)
        targets = [
            script,
            pow,
            phase,
            freqs,
        ]

        # Get args.
        action = f"""\
            python3 freqtag_fft.py
            --write_script_to {script}
            --write_pow_to {pow}
            --write_phase_to {phase}
            --write_freqs_to {freqs}
            --read_segmented_eeg_from {segmented_eeg}
            --read_sampling_rate_from {sampling_rate}
            --read_stimulus_start_from {stimulus_start}
            --read_stimulus_end_from {stimulus_end}
        """.split()

        # Go!
        yield dict(
            name=f"subject--{subject}",
            actions=[action],
            file_dep=sources,
            targets=targets,
        )
def task_freqtag_3d_fft() -> Dict:
    """
    Do an FFT of each and every trial whooooooooooooooooah man
    """
    Path(fname.freqtag_3d_fft_dir).mkdir(exist_ok=True, parents=True)

    for subject in SUBJECTS:
        # Get sources.
        segmented_eeg = fname.segmented_eeg(subject=subject)
        sampling_rate = fname.freqtag_sampling_rate
        stimulus_start = fname.freqtag_stimulus_start
        stimulus_end = fname.freqtag_stimulus_end
        sources = [
            segmented_eeg,
            sampling_rate,
            stimulus_start,
            stimulus_end
        ]

        # Get targets.
        script = fname.freqtag_3d_fft_script(subject=subject)
        spec = fname.freqtag_3d_fft_spec(subject=subject)
        targets = [
            script,
            spec,
        ]

        # Get args.
        action = f"""\
            python3 freqtag_3d_fft.py
            --write_script_to {script}
            --write_spec_to {spec}
            --read_segmented_eeg_from {segmented_eeg}
            --read_sampling_rate_from {sampling_rate}
            --read_stimulus_start_from {stimulus_start}
            --read_stimulus_end_from {stimulus_end}
        """.split()

        # Go!
        yield dict(
            name=f"subject--{subject}",
            actions=[action],
            file_dep=sources,
            targets=targets,
        )
def task_freqtag_sliding_window() -> Dict:
    """
    Do a sliding window average of our dataset.
    """
    Path(fname.freqtag_sliding_window_dir).mkdir(exist_ok=True, parents=True)

    for subject in SUBJECTS:
        for frequency in FREQUENCIES:
            # Get sources.
            sources_dict = dict(
                segmented_eeg=fname.segmented_eeg(subject=subject),
                sampling_rate=fname.freqtag_sampling_rate,
                stimulus_start=fname.freqtag_stimulus_start,
                stimulus_end=fname.freqtag_stimulus_end,
            )
            sources = list(sources_dict.values())

            # Get targets.
            targets_dict = dict(
                write_script_to=fname.freqtag_sliding_window_script(subject=subject, frequency=frequency),
                winmat3d=fname.freqtag_sliding_window_winmat3d(subject=subject, frequency=frequency),
                phasestabmat=fname.freqtag_sliding_window_phasestabmat(subject=subject, frequency=frequency),
                trialSNR=fname.freqtag_sliding_window_trialSNR(subject=subject, frequency=frequency),
                trialpow=fname.freqtag_sliding_window_trialpow(subject=subject, frequency=frequency),
                outfile=fname.freqtag_sliding_window_outfile(subject=subject, frequency=frequency),
                meanwinmat=fname.freqtag_sliding_window_meanwinmat(subject=subject, frequency=frequency),
            )
            targets = list(targets_dict.values())

            # Get args.
            action = f"""\
                python3 freqtag_sliding_window.py

                frequency={frequency}

                {encode_dict_in_bash(targets_dict)}
                {encode_dict_in_bash(sources_dict)}

            """.split()

            # Go!
            yield dict(
                name=f"subject--{subject}, frequency--{frequency}",
                actions=[action],
                file_dep=sources,
                targets=targets,
            )
def task_freqtag_fft_sliding_window() -> Dict:
    """
    Run an FFT on our sliding window results.
    """
    Path(fname.freqtag_fft_sliding_window_dir).mkdir(exist_ok=True, parents=True)

    for subject in SUBJECTS:
        for frequency in FREQUENCIES:
            # Get sources.
            sources_dict = dict(
                meanwinmat=fname.freqtag_sliding_window_meanwinmat(subject=subject, frequency=frequency),
            )
            sources = list(sources_dict.values())

            # Get targets.
            targets_dict = dict(
                write_script_to=fname.freqtag_fft_sliding_window_script(subject=subject, frequency=frequency),
                pow=fname.freqtag_fft_sliding_window_pow(subject=subject, frequency=frequency),
                phase=fname.freqtag_fft_sliding_window_phase(subject=subject, frequency=frequency),
                freqs=fname.freqtag_fft_sliding_window_freqs(subject=subject, frequency=frequency),
            )
            targets = list(targets_dict.values())

            # Get args.
            action = f"""\
                python3 freqtag_fft_sliding_window.py

                frequency={frequency}

                {encode_dict_in_bash(targets_dict)}
                {encode_dict_in_bash(sources_dict)}

            """.split()

            # Go!
            yield dict(
                name=f"subject--{subject}, frequency--{frequency}",
                actions=[action],
                file_dep=sources,
                targets=targets,
            )
def task_freqtag_better_sliding_window() -> Dict:
    """
    Run a sliding window using exactly correct frequencies rather than approximating them to 12Hz or 24Hz.
    """
    for subject in SUBJECTS:

        # Get sources.
        sources = dict(
            segmented_eeg=fname.segmented_eeg(subject=subject),
            sampling_rate=fname.freqtag_sampling_rate,
            stimulus_start=fname.freqtag_stimulus_start,
            stimulus_end=fname.freqtag_stimulus_end,
            frequencies=fname.eeg_flicker_frequencies(subject=subject, variable="frequencies")
        )
        sources_list = list(sources.values())

        # Get targets.
        targets = dict(
            write_script_to=fname.freqtag_better_sliding_window_script(subject=subject),
            winmat3d=fname.freqtag_better_sliding_window_winmat3d(subject=subject),
            phasestabmat=fname.freqtag_better_sliding_window_phasestabmat(subject=subject),
            trialSNR=fname.freqtag_better_sliding_window_trialSNR(subject=subject),
            trialpow=fname.freqtag_better_sliding_window_trialpow(subject=subject),
            meanwinmat=fname.freqtag_better_sliding_window_meanwinmat(subject=subject),
            oz_trialpow=fname.freqtag_better_sliding_window_channel(subject=subject, channel=20, variable="trialpow"),
            oz_trialSNR=fname.freqtag_better_sliding_window_channel(subject=subject, channel=20, variable="trialSNR"),
        )
        targets_list = list(targets.values())
        Path(targets["meanwinmat"]).parent.mkdir(exist_ok=True, parents=True)

        # Make the script to run.
        script = textwrap.dedent(f"""\
            %% This script runs a sliding window analysis. It's special in that it uses a different frequency for each trial in the analysis.
            do_one()

            function do_one()
                %% Load dataset.
                dataset = load_dataset('{sources["segmented_eeg"]}');

                %% Read input variables.
                load('{sources["stimulus_start"]}')
                load('{sources["stimulus_end"]}')
                load('{sources["sampling_rate"]}')
                load('{sources["frequencies"]}')

                %% Run functions.
                trialpow = [];
                winmat3d = [];
                phasestabmat = [];
                trialSNR = [];
                for i = 1:numel(dataset(1,1,:))
                    frequency = frequencies(i)
                    [minitrialpow,miniwinmat3d,miniphasestabmat,minitrialSNR] = freqtag_slidewin(dataset(:,:,i), 0, stimulus_start:stimulus_end, stimulus_start:stimulus_end, frequency, 600, sampling_rate, 'TEMP');
                    trialpow = [trialpow, minitrialpow];
                    winmat3d = [winmat3d, miniwinmat3d];
                    phasestabmat = [phasestabmat, miniphasestabmat];
                    trialSNR = [trialSNR, minitrialSNR];
                end
                meanwinmat = mean(winmat3d, 3);

                %% Save output variables.
                save('{targets["trialpow"]}', 'trialpow');
                oz_trialpow = trialpow(20,:);
                save('{targets["oz_trialpow"]}', 'oz_trialpow');
                oz_trialSNR = trialSNR(20,:);
                save('{targets["oz_trialSNR"]}', 'oz_trialSNR');
                save('{targets["winmat3d"]}', 'winmat3d');
                save('{targets["phasestabmat"]}', 'phasestabmat');
                save('{targets["trialSNR"]}', 'trialSNR');
                save('{targets["meanwinmat"]}', 'meanwinmat');
            end

            function [dataset] = load_dataset(path)
                % Load a dataset.
                [parent_dir, stem, suffix] = fileparts(path);
                eeglab;
                EEG = pop_loadset('filename', strcat(stem, suffix), 'filepath', parent_dir);
                EEG = eeg_checkset(EEG);
                dataset = double(EEG.data);
            end
            """)

        # Make action to run script.
        action = f"python3 matlab2.py".split()
        action += ["--script_contents", script]
        action += ["--write_script_to", targets["write_script_to"]]

        # Go!
        yield dict(
            name=f"subject--{subject}",
            actions=[action],
            file_dep=sources_list,
            targets=targets_list,
        )
def task_freqtag_better_fft_sliding_window() -> Dict:
    """
    Run an FFT on our better sliding window analysis.
    """
    Path(fname.freqtag_better_fft_sliding_window_dir).mkdir(exist_ok=True, parents=True)

    for subject in SUBJECTS:

        # Get sources.
        sources = dict(
            meanwinmat=fname.freqtag_better_sliding_window_meanwinmat(subject=subject),
        )
        sources_list = list(sources.values())

        # Get targets.
        targets = dict(
            write_script_to=fname.freqtag_better_fft_sliding_window_script(subject=subject),
            pow=fname.freqtag_better_fft_sliding_window_pow(subject=subject),
            phase=fname.freqtag_better_fft_sliding_window_phase(subject=subject),
            freqs=fname.freqtag_better_fft_sliding_window_freqs(subject=subject),
        )
        targets_list = list(targets.values())

        Path(targets["pow"]).parent.mkdir(exist_ok=True, parents=True)

        # Make MatLab script.
        script = textwrap.dedent(f"""\
        %% This script runs a fast fourier transform on our sliding window average.
        do_one()

        function do_one();
            %% Read input variables.
            load('{sources["meanwinmat"]}');

            %% Run functions.
            [pow, phase, freqs] = freqtag_FFT(meanwinmat, 600);

            %% Save output variables.
            save('{targets["pow"]}', 'pow');
            save('{targets["phase"]}', 'phase');
            save('{targets["freqs"]}', 'freqs');
        end
        """)

        # Make action to run script.
        action = f"python3 matlab2.py".split()
        action += ["--script_contents", script]
        action += ["--write_script_to", targets["write_script_to"]]

        # Go!
        yield dict(
            name=f"subject--{subject}",
            actions=[action],
            file_dep=sources_list,
            targets=targets_list,
        )
def task_freqtag_hilbert() -> Dict:
    """
    Run hilbert transform!
    """
    Path(fname.freqtag_hilbert_dir).mkdir(exist_ok=True, parents=True)

    for subject in SUBJECTS:
        for frequency in FREQUENCIES:
            # Get sources.
            sources_dict = dict(
                segmented_eeg=fname.segmented_eeg(subject=subject),
                sampling_rate=fname.freqtag_sampling_rate,
                stimulus_start=fname.freqtag_stimulus_start,
                stimulus_end=fname.freqtag_stimulus_end,
            )
            sources = list(sources_dict.values())

            # Get targets.
            targets_dict = dict(
                write_script_to=fname.freqtag_hilbert_script(subject=subject, frequency=frequency),
                powermat=fname.freqtag_hilbert_powermat(subject=subject, frequency=frequency),
                phasemat=fname.freqtag_hilbert_phasemat(subject=subject, frequency=frequency),
                complexmat=fname.freqtag_hilbert_complexmat(subject=subject, frequency=frequency),
            )
            targets = list(targets_dict.values())

            # Get args.
            action = f"""\
                python3 freqtag_hilbert.py

                frequency={frequency}
                electrode=20

                {encode_dict_in_bash(targets_dict)}
                {encode_dict_in_bash(sources_dict)}

            """.split()

            # Go!
            yield dict(
                name=f"subject--{subject}, frequency--{frequency}",
                actions=[action],
                file_dep=sources,
                targets=targets,
            )


# fMRI/EEG correlation tasks.
def task_correlate_whole_brain() -> Dict:
    """
    This is it! Huzzah! Correlate EEG and fMRI data across the whole brain. EEG data must be a 1xN or Nx1 .mat file.

    1st subbrick will be the correlation coefficient.
    2nd subbrick will be its p value.
    """
    def create_task(eeg_data: PathLike, in_image: PathLike, out_image: PathLike, name: str) -> dict:
        """
        Allows this task to easily be generalizable.

        Parameters
        ----------
        eeg_data : PathLike
            Path to a list of numbers in a 1xN or Nx1 MatLab .mat file
        in_image : PathLike
            Path to the fMRI image you want to correlate with the EEG data.
        out_image : PathLike
            Where you want to write your out image to.
        name : str
            Name of the task.
        """
        sources = dict(
            in_image_path=in_image,
            in_eeg_path=eeg_data,
        )

        targets = dict(
            out_image_path=out_image,
        )

        kwargs = {**sources, **targets}

        return dict(
            name=name,
            actions=[(correlate_whole_brain.main, [], kwargs)],
            file_dep=list(sources.values()),
            targets=list(targets.values()),
        )

    for subject in SUBJECTS:
        yield create_task(
                eeg_data=fname.eeg_trial_alpha(subject=subject, data="values"),
                out_image=fname.correlation_whole_brain_trials(subject=subject, variable="alpha"),
                in_image=fname.func_trial_amplitudes(subject=subject),
                name=f"trial alphas with trial funcs, sub--{subject}"
        )
        yield create_task(
                eeg_data=fname.freqtag_better_sliding_window_channel(subject=subject, channel=20, variable="trialpow"),
                out_image=fname.correlation_whole_brain_trials(subject=subject, variable="slidewinamp"),
                in_image=fname.func_trial_amplitudes(subject=subject),
                name=f"trial sliding window amplitudes with trial funcs, sub--{subject}"
        )
        yield create_task(
                eeg_data=fname.freqtag_better_sliding_window_channel(subject=subject, channel=20, variable="trialSNR"),
                out_image=fname.correlation_whole_brain_trials(subject=subject, variable="slidewinSNR"),
                in_image=fname.func_trial_amplitudes(subject=subject),
                name=f"trial sliding window SNR with trial funcs, sub--{subject}"
        )
        yield create_task(
                eeg_data=fname.canonical_trimmed(subject=subject),
                out_image=fname.correlation_whole_brain_canonical(subject=subject),
                in_image=fname.resampled_func(subject=subject),
                name=f"canonical BOLD, sub--{subject}"
        )
        for start_volume in EXPANDED_START_VOLUMES:
            for variable in "amplitudes SNRs".split():
                yield create_task(
                    eeg_data=fname.eeg_sliding_sliding_window_improved(subject=subject, variable=variable),
                    out_image=fname.correlation_whole_brain_improved(subject=subject, start_volume=start_volume, variable=variable),
                    in_image=fname.final_func(subject=subject, start_volume=start_volume),
                    name=f"improved sliding sliding window, sub--{subject}, startvolume--{start_volume}, variable--{variable}"
                )
        for start_volume in START_VOLUMES:
            for alpha_data in "values SNRs".split():
                yield create_task(
                    eeg_data=fname.eeg_alpha(subject=subject, data=alpha_data),
                    out_image=fname.correlation_whole_brain_alpha(subject=subject, start_volume=start_volume, data=alpha_data),
                    in_image=fname.final_func(subject=subject, start_volume=start_volume),
                    name=f"alpha, {alpha_data}, sub--{subject}, startvolume--{start_volume}"
                )
        for permutation in PERMUTATIONS:
            start_volume = 4
            analysis = "alpha"
            for variable in "amplitude SNR".split():
                yield create_task(
                    eeg_data=fname.scrambled_series(subject=subject, start_volume=start_volume, variable=variable, analysis=analysis, permutation=permutation),
                    out_image=fname.correlation_whole_brain_permutation(subject=subject, start_volume=start_volume, variable=variable, analysis=analysis, permutation=permutation),
                    in_image=fname.final_func(subject=subject, start_volume=start_volume),
                    name=f"permutation correlation, sub--{subject}, startvolume--{start_volume}, variable--{variable}, analysis--{analysis}, permutation--{permutation}"
                )
            variable = "amplitude"
            start_volume = "na"
            yield create_task(
                eeg_data=fname.scrambled_series(subject=subject, start_volume=start_volume, variable=variable, analysis=analysis, permutation=permutation),
                out_image=fname.correlation_whole_brain_permutation(subject=subject, start_volume=start_volume, variable=variable, analysis=analysis, permutation=permutation),
                in_image=fname.func_trial_amplitudes(subject=subject),
                name=f"permutation correlation, sub--{subject}, startvolume--{start_volume}, variable--{variable}, analysis--{analysis}, permutation--{permutation}"
            )

            start_volume = 5
            analysis = "ssvep"
            for variable in "amplitude SNR".split():
                yield create_task(
                    eeg_data=fname.scrambled_series(subject=subject, start_volume=start_volume, variable=variable, analysis=analysis, permutation=permutation),
                    out_image=fname.correlation_whole_brain_permutation(subject=subject, start_volume=start_volume, variable=variable, analysis=analysis, permutation=permutation),
                    in_image=fname.final_func(subject=subject, start_volume=start_volume),
                    name=f"permutation correlation, sub--{subject}, startvolume--{start_volume}, variable--{variable}, analysis--{analysis}, permutation--{permutation}"
                )
            start_volume = "na"
            variable = "amplitude"
            yield create_task(
                eeg_data=fname.scrambled_series(subject=subject, start_volume=start_volume, variable=variable, analysis=analysis, permutation=permutation),
                out_image=fname.correlation_whole_brain_permutation(subject=subject, start_volume=start_volume, variable=variable, analysis=analysis, permutation=permutation),
                in_image=fname.func_trial_amplitudes(subject=subject),
                name=f"permutation correlation trials, sub--{subject}, startvolume--{start_volume}, variable--{variable}, analysis--{analysis}, permutation--{permutation}"
            )
def task_fisher_transform_whole_brain() -> Dict:
    """
    Apply the Fisher Z transform to our whole brain correlation so we can do more valid t-tests.
    """
    def create_task(in_correlation_path: PathLike, out_fisher_path: PathLike, name: str) -> dict:
        """
        Allows this task to easily be generalizable.
        """
        sources = dict(
            in_correlation_path=in_correlation_path,
        )

        targets = dict(
            out_fisher_path=out_fisher_path,
        )

        kwargs = {**sources, **targets}

        return dict(
            name=name,
            actions=[(fisher_transform.main, [], kwargs)],
            file_dep=list(sources.values()),
            targets=list(targets.values()),
        )

    for subject in SUBJECTS:
        for start_volume in EXPANDED_START_VOLUMES:
            for variable in "amplitudes SNRs".split():
                analysis = "slidslidwin_improved"
                yield create_task(
                    in_correlation_path=fname.correlation_whole_brain_improved(subject=subject, start_volume=start_volume, variable=variable),
                    out_fisher_path=fname.correlation_fisher(subject=subject, start_volume=start_volume, variable=variable, analysis=analysis),
                    name=f"sub--{subject}, startvolume--{start_volume}, variable--{variable}, analysis--{analysis}",
                )
        for start_volume in START_VOLUMES:
            for variable in "values SNRs".split():
                analysis = "alpha"
                yield create_task(
                    in_correlation_path=fname.correlation_whole_brain_alpha(subject=subject, start_volume=start_volume, data=variable),
                    out_fisher_path=fname.correlation_fisher(subject=subject, start_volume=start_volume, variable=variable, analysis=analysis),
                    name=f"sub--{subject}, startvolume--{start_volume}, variable--{variable}, analysis--{analysis}",
                )
        for start_volume in "na".split():
            for variable in "alpha slidewinamp".split():
                analysis = "trials"
                yield create_task(
                    in_correlation_path=fname.correlation_whole_brain_trials(subject=subject, variable=variable),
                    out_fisher_path=fname.correlation_fisher(subject=subject, start_volume=start_volume, variable=variable, analysis=analysis),
                    name=f"sub--{subject}, startvolume--{start_volume}, variable--{variable}, analysis--{analysis}",
                )
def task_ttest_whole_brain_fishers() -> Dict:
    """
    t-test the fisher transformed correlations.
    """
    def create_task(images: PathLike, out_path: PathLike, name: str) -> dict:
        """
        Allows this task to easily be generalizable.

        Parameters
        ----------
        images : PathLike
            List of paths to images to ttest against each other.
        out_path : PathLike
            Where to write our ttest results to.
        name : str
            Name of the task.
        """
        sources = dict(
            images=images,
        )

        targets = dict(
            ttest=out_path,
        )

        kwargs = dict(
            images=[f"{path}[0]" for path in sources["images"]],
            prefix=targets["ttest"]
        )

        return dict(
            name=name,
            actions=[(ttest.main, [], kwargs)],
            file_dep=list(sources.values())[0],
            targets=list(targets.values()),
        )

    for start_volume in EXPANDED_START_VOLUMES:
        for variable in "amplitudes SNRs".split():
            analysis = "slidslidwin_improved"
            yield create_task(
                images=[fname.correlation_fisher(subject=subject, start_volume=start_volume, variable=variable, analysis=analysis) for subject in SUBJECTS],
                out_path=fname.correlations_ttest_fisher(start_volume=start_volume, variable=variable, analysis=analysis),
                name=f"startvolume--{start_volume}, variable--{variable}, analysis--{analysis}",
            )
    for start_volume in START_VOLUMES:
        for variable in "values SNRs".split():
            analysis = "alpha"
            yield create_task(
                images=[fname.correlation_fisher(subject=subject, start_volume=start_volume, variable=variable, analysis=analysis) for subject in SUBJECTS],
                out_path=fname.correlations_ttest_fisher(start_volume=start_volume, variable=variable, analysis=analysis),
                name=f"startvolume--{start_volume}, variable--{variable}, analysis--{analysis}",
            )
    for start_volume in "na".split():
        for variable in "alpha slidewinamp".split():
            analysis = "trials"
            yield create_task(
                images=[fname.correlation_fisher(subject=subject, start_volume=start_volume, variable=variable, analysis=analysis) for subject in SUBJECTS],
                out_path=fname.correlations_ttest_fisher(start_volume=start_volume, variable=variable, analysis=analysis),
                name=f"startvolume--{start_volume}, variable--{variable}, analysis--{analysis}",
            )
def task_cohens_d_whole_brain() -> Dict:
    """
    Calculate Cohen's d for each voxel in our ttests.
    """
    def create_task(in_ttest: PathLike, out_cohen: PathLike, name: str) -> dict:
        """
        Allows this task to easily be generalizable.
        """
        sources = dict(
            in_ttest=in_ttest,
        )

        targets = dict(
            out_cohen=out_cohen,
        )

        kwargs = {**sources, **targets, "n_samples": len(SUBJECTS)}

        return dict(
            name=name,
            actions=[(cohens.main, [], kwargs)],
            file_dep=list(sources.values()),
            targets=list(targets.values()),
        )

    for start_volume in EXPANDED_START_VOLUMES:
        for variable in "amplitudes SNRs".split():
            analysis = "slidslidwin_improved"
            yield create_task(
                in_ttest=fname.correlations_ttest_fisher(start_volume=start_volume, variable=variable, analysis=analysis),
                out_cohen=fname.correlations_cohens_fisher(start_volume=start_volume, variable=variable, analysis=analysis),
                name=f"startvolume--{start_volume}, variable--{variable}, analysis--{analysis}",
            )
    for start_volume in START_VOLUMES:
        for variable in "values SNRs".split():
            analysis = "alpha"
            yield create_task(
                in_ttest=fname.correlations_ttest_fisher(start_volume=start_volume, variable=variable, analysis=analysis),
                out_cohen=fname.correlations_cohens_fisher(start_volume=start_volume, variable=variable, analysis=analysis),
                name=f"startvolume--{start_volume}, variable--{variable}, analysis--{analysis}",
            )
    for start_volume in "na".split():
        for variable in "alpha slidewinamp".split():
            analysis = "trials"
            yield create_task(
                in_ttest=fname.correlations_ttest_fisher(start_volume=start_volume, variable=variable, analysis=analysis),
                out_cohen=fname.correlations_cohens_fisher(start_volume=start_volume, variable=variable, analysis=analysis),
                name=f"startvolume--{start_volume}, variable--{variable}, analysis--{analysis}",
            )
def task_ttest_whole_brain_correlations() -> Dict:
    """
    ttest the correlations we calculated.
    """
    def create_task(images: List[PathLike], out_path: PathLike, name: str) -> dict:
        """
        Allows this task to easily be generalizable.

        Parameters
        ----------
        images : List[PathLike]
            List of paths to images to ttest against each other.
        out_path : PathLike
            Where to write our ttest results to.
        name : str
            Name of the task.
        """
        sources = dict(
            images=images,
        )

        targets = dict(
            ttest=out_path,
        )

        kwargs = dict(
            images=[f"{path}[0]" for path in sources["images"]],
            prefix=get_prefix(targets["ttest"])
        )

        return dict(
            name=name,
            actions=[(ttest.main, [], kwargs)],
            file_dep=list(sources.values())[0],
            targets=list(targets.values()),
        )

    yield create_task(
        images=[fname.correlation_whole_brain_canonical(subject=subject) for subject in SUBJECTS],
        out_path=fname.correlations_whole_brain_canonical_ttest,
        name=f"canonical BOLD",
    )
    for variable in "alpha slidewinamp slidewinSNR".split():
        yield create_task(
            images=[fname.correlation_whole_brain_trials(subject=subject, variable=variable) for subject in SUBJECTS],
            out_path=fname.correlations_whole_brain_trials_ttest(variable=variable),
            name=f"trials, variable--{variable}",
        )
    for start_volume in EXPANDED_START_VOLUMES:
        for variable in "amplitudes SNRs".split():
            yield create_task(
                images=[fname.correlation_whole_brain_improved(subject=subject, start_volume=start_volume, variable=variable) for subject in SUBJECTS],
                out_path=fname.correlations_improved_whole_brain_ttest(start_volume=start_volume, variable=variable),
                name=f"improved sliding sliding window, startvolume--{start_volume}, variable--{variable}",
            )
    for start_volume in START_VOLUMES:
        for alpha_data in "values SNRs".split():
            yield create_task(
                images=[fname.correlation_whole_brain_alpha(subject=subject, start_volume=start_volume, data=alpha_data) for subject in SUBJECTS],
                out_path=fname.correlations_whole_brain_alpha_ttest(start_volume=start_volume, data=alpha_data),
                name=f"alphas, data--{alpha_data}, startvolume--{start_volume}",
            )
    for permutation in PERMUTATIONS:
        start_volume = 4
        variable = "amplitude"
        analysis = "alpha"
        yield create_task(
            images=[fname.correlation_whole_brain_permutation(subject=subject, start_volume=start_volume, variable=variable, analysis=analysis, permutation=permutation) for subject in SUBJECTS],
            out_path=fname.correlations_whole_brain_permutations_ttest(start_volume=start_volume, variable=variable, analysis=analysis, permutation=permutation),
            name=f"permutation ttest, startvolume--{start_volume}, variable--{variable}, analysis--{analysis}, permutation--{permutation}",
        )

        start_volume = 5
        analysis = "ssvep"
        yield create_task(
            images=[fname.correlation_whole_brain_permutation(subject=subject, start_volume=start_volume, variable=variable, analysis=analysis, permutation=permutation) for subject in SUBJECTS],
            out_path=fname.correlations_whole_brain_permutations_ttest(start_volume=start_volume, variable=variable, analysis=analysis, permutation=permutation),
            name=f"permutation ttest, startvolume--{start_volume}, variable--{variable}, analysis--{analysis}, permutation--{permutation}",
        )
def task_correlate_eeg_with_average_microregion_timeseries() -> Dict:
    """
    Correlate the time series of each microregion with EEG data.

    Inputs
    ------
    eeg_sliding_sliding_window_oz_amplitudes
        Amplitudes of the Oz electrode calculated with a sliding sliding window using frequency estimates of 12Hz and 24Hz.
    eeg_sliding_sliding_window_oz_SNR
        Signal to noise ratio of the Oz electrode calculated with a sliding sliding window using frequency estimates of 12Hz and 24Hz.
    """
    def create_task(sources: Dict[str, PathLike], targets: Dict[str, PathLike], name: str) -> dict:
        """
        Allows this task to easily be generalizable.

        Parameters
        ----------
        sources : Dict
            load_microROI_from : PathLike
                Path to a text file containing the average fMRI time series of a region.
            load_amplitudes_from : PathLike
                Path to a 1xN or Nx1 MatLab .mat file containing a list of numbers.
        targets : Dict
            save_spearman_to : PathLike
                Where to output our correlation results to.
            save_scatter_to : PathLike
                Where to output our scatter plot to.
            save_table_to : PathLike
                Where to output a table containing the EEG data and average fMRI time series.
        """
        kwargs = {**sources, **targets}

        return dict(
            name=name,
            actions=[(correlate_regions.main, [], kwargs)],
            file_dep=list(sources.values()),
            targets=list(targets.values()),
        )

    for subject in SUBJECTS:
        for region in "calcarine occipital".split():
            for start_volume in START_VOLUMES:
                for variable in "amplitudes SNRs".split():
                    yield create_task(
                        sources=dict(
                            load_microROI_from=fname.microregion_average(subject=subject, mask=region, start_volume=start_volume),
                            load_amplitudes_from=fname.eeg_sliding_sliding_window_improved(subject=subject, variable=variable),
                        ),
                        targets=dict(
                            save_spearman_to=fname.microregions_correlation_better(subject=subject, mask=region, start_volume=start_volume, variable=variable),
                            save_table_to=fname.microregions_and_amplitudes_better(subject=subject, mask=region, start_volume=start_volume, variable=variable),
                            save_scatter_to=fname.microregions_correlation_scatter_plot_better(subject=subject, mask=region, start_volume=start_volume, variable=variable),
                        ),
                        name=f"improved sliding sliding window, variable--{variable}, subject--{subject}, mask--{region}, start_volume--{start_volume}"
                    )
                for alpha_data in "values SNRs".split():
                    yield create_task(
                        sources=dict(
                            load_microROI_from=fname.microregion_average(subject=subject, mask=region, start_volume=start_volume),
                            load_amplitudes_from=fname.eeg_alpha(subject=subject, data=alpha_data),
                        ),
                        targets=dict(
                            save_spearman_to=fname.microregions_correlation_alpha_results(subject=subject, mask=region, start_volume=start_volume, data=alpha_data),
                            save_table_to=fname.microregions_and_alpha_amplitudes(subject=subject, mask=region, start_volume=start_volume, data=alpha_data),
                            save_scatter_to=fname.microregions_correlation_alpha_scatter_plot(subject=subject, mask=region, start_volume=start_volume, data=alpha_data),
                        ),
                        name=f"alpha, data--{alpha_data}, subject--{subject}, mask--{region}, start_volume--{start_volume}"
                    )
def task_correlate_eeg_with_average_microregion_timeseries_across_subjects() -> Dict:
    """
    How big is the correlation across ALL subjects for each microROI?
    """
    def create_task(sources: Dict[str, PathLike], targets: Dict[str, PathLike], name: str) -> dict:
        """
        Allows this task to easily be generalizable.

        Parameters
        ----------
        sources : Dict
            load_tables_from : List[PathLike]
                Where to get each of our microROI+amplitude tables
        targets : Dict
            save_scatter_to : PathLike
                Where to save our scatter plot.
            save_table_to : PathLike
                Where to save the table of data we make from our tiny tables.
            save_spearman_to : PathLike
                Where to save our Spearman results.
        """
        kwargs = {**sources, **targets}

        return dict(
            name=name,
            actions=[(correlate_tables.main, [], kwargs)],
            file_dep=list(sources.values())[0],
            targets=list(targets.values()),
        )

    for region in "calcarine occipital".split():
        for start_volume in START_VOLUMES:
            for variable in "amplitudes SNRs".split():
                yield create_task(
                    sources=dict(
                        load_tables_from=[fname.microregions_and_amplitudes_better(subject=subject, mask=region, start_volume=start_volume, variable=variable) for subject in SUBJECTS],
                    ),
                    targets=dict(
                        save_scatter_to=fname.correlation_across_subjects_scatter_better(mask=region, variable=variable, start_volume=start_volume),
                        save_table_to=fname.correlation_across_subjects_table_better(mask=region, variable=variable, start_volume=start_volume),
                        save_spearman_to=fname.correlation_across_subjects_better(mask=region, variable=variable, start_volume=start_volume),
                    ),
                    name=f"better sliding sliding window, variable--{variable}, mask--{region}, start_volume--{start_volume}",
                )
            for alpha_data in "values SNRs".split():
                yield create_task(
                    sources=dict(
                        load_tables_from=[fname.microregions_and_alpha_amplitudes(subject=subject, mask=region, start_volume=start_volume, data=alpha_data) for subject in SUBJECTS],
                    ),
                    targets=dict(
                        save_scatter_to=fname.correlation_across_subjects_alpha_scatter(mask=region, data=alpha_data, start_volume=start_volume),
                        save_table_to=fname.correlation_across_subjects_alpha_table(mask=region, data=alpha_data, start_volume=start_volume),
                        save_spearman_to=fname.correlation_across_subjects_alpha(mask=region, data=alpha_data, start_volume=start_volume),
                    ),
                    name=f"alpha, data--{alpha_data}, mask--{region}, start_volume--{start_volume}",
                )
def task_mask_and_average_correlations() -> Dict:
    """
    Get the average correlation within a specific ROI.

    Uses AFNI's 3dmaskave: https://afni.nimh.nih.gov/pub/dist/doc/htmldoc/programs/3dmaskave_sphx.html#ahelp-3dmaskave
    """
    def create_task(from_image: PathLike, to_text_file: PathLike, from_mask: PathLike, name: str):
        """
        Allows this task to easily be generalizable.
        """
        sources = dict(from_image=from_image, from_mask=from_mask)
        targets = dict(to_text_file=to_text_file)
        kwargs = {**sources, **targets}

        return dict(
            name=name,
            actions=[(maskave.main, [], kwargs)],
            file_dep=list(sources.values()),
            targets=list(targets.values()),
        )

    for subject in SUBJECTS:
        for mask in "calcarine occipital".split():
            for start_volume in START_VOLUMES:
                for data in "values SNRs".split():
                    yield create_task(
                        from_image=fname.correlation_whole_brain_alpha(subject=subject, start_volume=start_volume, data=data),
                        to_text_file=fname.correlation_alpha_average(subject=subject, start_volume=start_volume, mask=mask, data=data),
                        from_mask=fname.micromask(subject=subject, mask=mask),
                        name=f"sub--{subject}, startvolume--{start_volume}, mask--{mask}, data--{data}",
                    )
            for start_volume in EXPANDED_START_VOLUMES:
                for variable in "amplitudes SNRs".split():
                    yield create_task(
                        from_image=fname.correlation_whole_brain_improved(subject=subject, start_volume=start_volume, variable=variable),
                        to_text_file=fname.correlation_improved_sliding_sliding_window_average(subject=subject, start_volume=start_volume, mask=mask, variable=variable),
                        from_mask=fname.micromask(subject=subject, mask=mask),
                        name=f"sub--{subject}, startvolume--{start_volume}, mask--{mask}, variable--{variable}",
                    )
            for permutation in PERMUTATIONS:
                start_volume = 4
                variable = "amplitude"
                analysis = "alpha"
                yield create_task(
                    from_image=fname.correlation_whole_brain_permutation(subject=subject, start_volume=start_volume, variable=variable, analysis=analysis, permutation=permutation),
                    to_text_file=fname.correlation_permutation_average(subject=subject, start_volume=start_volume, mask=mask, variable=variable, analysis=analysis, permutation=permutation),
                    from_mask=fname.micromask(subject=subject, mask=mask),
                    name=f"sub--{subject}, startvolume--{start_volume}, mask--{mask}, variable--{variable}, analysis--{analysis}, permutation--{permutation}",
                )

                start_volume = 5
                analysis = "ssvep"
                yield create_task(
                    from_image=fname.correlation_whole_brain_permutation(subject=subject, start_volume=start_volume, variable=variable, analysis=analysis, permutation=permutation),
                    to_text_file=fname.correlation_permutation_average(subject=subject, start_volume=start_volume, mask=mask, variable=variable, analysis=analysis, permutation=permutation),
                    from_mask=fname.micromask(subject=subject, mask=mask),
                    name=f"sub--{subject}, startvolume--{start_volume}, mask--{mask}, variable--{variable}, analysis--{analysis}, permutation--{permutation}",
                )
def task_average_averages() -> Dict:
    """
    Average the ROI correlation averages. (Do across subjects.)

    Returns:
        Dict: Dummy pydoit dict.
    """
    def create_task(from_text_files: List[PathLike], to_table: PathLike, to_result: PathLike, name: str) -> Dict:
        """
        Average ROI correlations from a list of ROI correlation results.

        Args:
            from_text_files (List[PathLike]): List of ROI correlation results for individual subjects.
            to_table (PathLike): Big table to help us check that we got the right results.
            to_result (PathLike): Average correlation across subjects.
            name (str): What to name pydoit subtask.

        Returns:
            Dict: Pydoit dummy dict.
        """
        sources = dict(from_text_files=from_text_files)
        targets = dict(to_table=to_table, to_result=to_result)
        kwargs = {**sources, **targets}

        return dict(
            name=name,
            actions=[(mean_means.main, [], kwargs)],
            file_dep=list(sources.values())[0],
            targets=list(targets.values()),
        )

    for mask in "calcarine occipital".split():
        for start_volume in START_VOLUMES:
            analysis = "alpha"
            variable = "amplitude"
            yield create_task(
                from_text_files=[fname.correlation_alpha_average(subject=subject, start_volume=start_volume, mask=mask, data="values") for subject in SUBJECTS],
                to_table=fname.average_averages_correlation(variable=variable, start_volume=start_volume, mask=mask, analysis=analysis, outfile="table"),
                to_result=fname.average_averages_correlation(variable=variable, start_volume=start_volume, mask=mask, analysis=analysis, outfile="results"),
                name=f"startvolume--{start_volume}, mask--{mask}, variable--{variable}, analysis--{analysis}",
            )
        for start_volume in EXPANDED_START_VOLUMES:
            variable = "amplitude"
            analysis = "ssvep"
            yield create_task(
                from_text_files=[fname.correlation_improved_sliding_sliding_window_average(subject=subject, start_volume=start_volume, mask=mask, variable="amplitudes") for subject in SUBJECTS],
                to_table=fname.average_averages_correlation(variable=variable, start_volume=start_volume, mask=mask, analysis=analysis, outfile="table"),
                to_result=fname.average_averages_correlation(variable=variable, start_volume=start_volume, mask=mask, analysis=analysis, outfile="results"),
                name=f"startvolume--{start_volume}, mask--{mask}, variable--{variable}, analysis--{analysis}",
            )
        for permutation in PERMUTATIONS:
            start_volume = 4
            variable = "amplitude"
            analysis = "alpha"
            yield create_task(
                from_text_files=[fname.correlation_permutation_average(subject=subject, analysis=analysis, start_volume=start_volume, mask=mask, variable=variable, permutation=permutation) for subject in SUBJECTS],
                to_table=fname.average_averages_correlation_permutation(permutation=permutation, variable=variable, start_volume=start_volume, mask=mask, analysis=analysis, outfile="table"),
                to_result=fname.average_averages_correlation_permutation(permutation=permutation, variable=variable, start_volume=start_volume, mask=mask, analysis=analysis, outfile="results"),
                name=f"startvolume--{start_volume}, mask--{mask}, variable--{variable}, analysis--{analysis}, permutation--{permutation}",
            )

            start_volume = 5
            analysis = "ssvep"
            yield create_task(
                from_text_files=[fname.correlation_permutation_average(subject=subject, analysis=analysis, start_volume=start_volume, mask=mask, variable=variable, permutation=permutation) for subject in SUBJECTS],
                to_table=fname.average_averages_correlation_permutation(permutation=permutation, variable=variable, start_volume=start_volume, mask=mask, analysis=analysis, outfile="table"),
                to_result=fname.average_averages_correlation_permutation(permutation=permutation, variable=variable, start_volume=start_volume, mask=mask, analysis=analysis, outfile="results"),
                name=f"startvolume--{start_volume}, mask--{mask}, variable--{variable}, analysis--{analysis}, permutation--{permutation}",
            )
def task_ttest_averages() -> Dict:
    """
    ttest the correlation averages.
    """
    def create_task(from_text_files: List[PathLike], to_table: PathLike, to_statistics: PathLike, name: str) -> Dict:
        """
        Allows this task to easily be generalizable.
        """
        sources = dict(from_text_files=from_text_files)
        targets = dict(to_table=to_table, to_statistics=to_statistics)
        kwargs = {**sources, **targets}

        return dict(
            name=name,
            actions=[(ttest_averages.main, [], kwargs)],
            file_dep=list(sources.values())[0],
            targets=list(targets.values()),
        )

    for mask in "calcarine occipital".split():
        for start_volume in START_VOLUMES:
            for data in "values SNRs".split():
                yield create_task(
                    from_text_files=[fname.correlation_alpha_average(subject=subject, start_volume=start_volume, mask=mask, data=data) for subject in SUBJECTS],
                    to_table=fname.ttest_averages(variable=data, start_volume=start_volume, mask=mask, analysis="alpha", outfile="table"),
                    to_statistics=fname.ttest_averages(variable=data, start_volume=start_volume, mask=mask, analysis="alpha", outfile="results"),
                    name=f"startvolume--{start_volume}, mask--{mask}, variable--{data}, analysis--alpha",
                )
        for start_volume in EXPANDED_START_VOLUMES:
            for variable in "amplitudes SNRs".split():
                yield create_task(
                    from_text_files=[fname.correlation_improved_sliding_sliding_window_average(subject=subject, start_volume=start_volume, mask=mask, variable=variable) for subject in SUBJECTS],
                    to_table=fname.ttest_averages(variable=variable, start_volume=start_volume, mask=mask, analysis="improvedslidslidwin", outfile="table"),
                    to_statistics=fname.ttest_averages(variable=variable, start_volume=start_volume, mask=mask, analysis="improvedslidslidwin", outfile="results"),
                    name=f"startvolume--{start_volume}, mask--{mask}, variable--{variable}",
                )
def task_get_canonical_bold() -> Dict:
    """
    For each subject, calculate the canonical BOLD expected based on stimulus times alone.
    """
    def create_task(from_onsets: PathLike, to_trs: PathLike, write_script_to: PathLike, name: str) -> Dict:
        """Make this task generalizable."""
        sources = dict(from_onsets=from_onsets)
        targets = dict(to_trs=to_trs)
        kwargs = {**sources, **targets, "write_script_to": write_script_to}

        return dict(
            name=name,
            actions=[(get_canonical.main, [], kwargs)],
            file_dep=list(sources.values()),
            targets=list(targets.values()),
        )

    for subject in SUBJECTS:
        yield create_task(
            from_onsets=fname.afniproc_onsets(subject=subject),
            to_trs=fname.canonical(subject=subject),
            write_script_to=fname.canonical_script(subject=subject),
            name=f"subject--{subject}"
        )
def task_trim_canonical_bold() -> Dict:
    """
    Trim canonical bolds to mirror the afni_proc.py outputs.
    """
    def create_task(from_mat: PathLike, to_trimmed_mat: PathLike, name: str, start_index: int = None, end_index: int = None) -> Dict:
        """
        Make this task generalizable.
        """
        sources = dict(from_mat=from_mat)
        targets = dict(to_trimmed_mat=to_trimmed_mat)
        kwargs = {**sources, **targets, "start_index": start_index, "end_index": end_index}

        return dict(
            name=name,
            actions=[(trim_mat.main, [], kwargs)],
            file_dep=list(sources.values()),
            targets=list(targets.values()),
        )

    for subject in SUBJECTS:
        yield create_task(
            from_mat=fname.canonical(subject=subject),
            to_trimmed_mat=fname.canonical_trimmed(subject=subject),
            end_index=-2,
            name=f"subject--{subject}"
        )
def task_subtract_canonical_bold() -> Dict:
    """
    Subtract the mean canonical BOLD correlation from the mean alpha and mean ssVEP correlations.
    """
    def create_task(minuend: PathLike, out_prefix: PathLike, name: str, subtrahend: PathLike = fname.correlations_whole_brain_canonical_ttest) -> Dict:
        """
        Make this task generalizable.
        """
        sources = dict(minuend=minuend, subtrahend=subtrahend)
        targets = dict(out_prefix=out_prefix)
        kwargs = {
            "minuend": minuend + "[0]",
            "subtrahend": subtrahend + "[0]",
            "out_prefix": out_prefix,
        }

        return dict(
            name=name,
            actions=[(compare_images.main, [], kwargs)],
            file_dep=list(sources.values()),
            targets=list(targets.values()),
        )

    for start_volume in START_VOLUMES:
        old_and_new_variables = {"values": "amplitude", "SNRs": "SNR"}
        for old_variable, new_variable in old_and_new_variables.items():
            yield create_task(
                minuend=fname.correlations_whole_brain_alpha_ttest(start_volume=start_volume, data=old_variable),
                out_prefix=fname.compared_to_canonical(start_volume=start_volume, variable=new_variable, analysis="alpha"),
                name=f"analysis--alpha, variable--{new_variable}, start_volume--{start_volume}",
            )
    for start_volume in EXPANDED_START_VOLUMES:
        old_and_new_variables = {"amplitudes": "amplitude", "SNRs": "SNR"}
        for old_variable, new_variable in old_and_new_variables.items():
            yield create_task(
                minuend=fname.correlations_improved_whole_brain_ttest(start_volume=start_volume, variable=old_variable),
                out_prefix=fname.compared_to_canonical(start_volume=start_volume, variable=new_variable, analysis="ssvep"),
                name=f"analysis--ssvep, variable--{new_variable}, start_volume--{start_volume}",
            )
    for permutation in PERMUTATIONS:
        start_volume = 4
        variable = "amplitude"
        analysis = "alpha"
        yield create_task(
            minuend=fname.correlations_whole_brain_permutations_ttest(start_volume=start_volume, variable=variable, analysis=analysis, permutation=permutation),
            out_prefix=fname.compared_permutations_to_canonical(start_volume=start_volume, variable=variable, analysis=analysis, permutation=permutation),
            name=f"permutations, analysis--{analysis}, variable--{variable}, start_volume--{start_volume}, permutation--{permutation}",
        )

        start_volume = 5
        analysis = "ssvep"
        yield create_task(
            minuend=fname.correlations_whole_brain_permutations_ttest(start_volume=start_volume, variable=variable, analysis=analysis, permutation=permutation),
            out_prefix=fname.compared_permutations_to_canonical(start_volume=start_volume, variable=variable, analysis=analysis, permutation=permutation),
            name=f"permutations, analysis--{analysis}, variable--{variable}, start_volume--{start_volume}, permutation--{permutation}",
        )
def task_calculate_variance() -> Dict:
    """
    Square our correlation images to get them as variance.
    """
    def create_task(in_correlation_image: PathLike, out_variance_image: PathLike, name: str) -> Dict:
        """
        Square our correlation images to get them as variance.

        Assumes correlations are located in 1st subbrick of correlation image.

        Args:
            in_correlation_image (PathLike): Path to a correlation image.
            out_variance_image (PathLike): Where to write the squared image.
            name (str): What to name the task.

        Returns:
            Dict: Dummy pydoit dict.
        """
        out_parent_dir = Path(out_variance_image).parent
        return dict(
            name=name,
            file_dep=[in_correlation_image],
            targets=[out_variance_image],
            actions=[f"mkdir -p '{out_parent_dir}'",
                    f"3dcalc -float -a '{in_correlation_image}[0]' -expr 'a^2' -prefix '{out_variance_image}'"]
        )

    # Canonical BOLD.
    start_volume = "na"
    variable = "na"
    analysis = "baseline"
    yield create_task(
        in_correlation_image=fname.correlations_whole_brain_canonical_ttest,
        out_variance_image=fname.variance_whole_brain(start_volume=start_volume, variable=variable, analysis=analysis),
        name=f"analysis--{analysis}, variable--{variable}, start_volume--{start_volume}",
    )

    # Trial by trial alpha.
    analysis = "alpha"
    start_volume = "na"
    variable = "amplitude"
    yield create_task(
        in_correlation_image=fname.correlations_whole_brain_trials_ttest(variable=analysis),
        out_variance_image=fname.variance_whole_brain(start_volume=start_volume, variable=variable, analysis=analysis),
        name=f"analysis--{analysis}, variable--{variable}, start_volume--{start_volume}",
    )

    # Continuous alpha.
    old_and_new_variables = {"values": "amplitude", "SNRs": "SNR"}
    for start_volume in START_VOLUMES:
        for old_variable, new_variable in old_and_new_variables.items():
            yield create_task(
                in_correlation_image=fname.correlations_whole_brain_alpha_ttest(start_volume=start_volume, data=old_variable),
                out_variance_image=fname.variance_whole_brain(start_volume=start_volume, variable=new_variable, analysis=analysis),
                name=f"analysis--{analysis}, variable--{new_variable}, start_volume--{start_volume}",
            )

    # Continuous ssVEP.
    analysis = "ssvep"
    old_and_new_variables = {"amplitudes": "amplitude", "SNRs": "SNR"}
    for start_volume in EXPANDED_START_VOLUMES:
        for old_variable, new_variable in old_and_new_variables.items():
            yield create_task(
                in_correlation_image=fname.correlations_improved_whole_brain_ttest(start_volume=start_volume, variable=old_variable),
                out_variance_image=fname.variance_whole_brain(start_volume=start_volume, variable=new_variable, analysis=analysis),
                name=f"analysis--{analysis}, variable--{new_variable}, start_volume--{start_volume}",
            )
    
    # Trial by trial ssVEP.
    start_volume = "na"
    old_and_new_variables = {"slidewinamp": "amplitude", "slidewinSNR": "SNR"}
    for old_variable, new_variable in old_and_new_variables.items():
        yield create_task(
            in_correlation_image=fname.correlations_whole_brain_trials_ttest(variable=old_variable),
            out_variance_image=fname.variance_whole_brain(start_volume=start_volume, variable=new_variable, analysis=analysis),
            name=f"analysis--{analysis}, variable--{new_variable}, start_volume--{start_volume}",
        )
def task_subtract_canonical_variance():
    """
    Subtract the canonical variance from the actual variance.
    """
    def create_task(in_actual_variance: PathLike, in_baseline_variance: PathLike, out_corrected_variance: PathLike, name: str) -> Dict:
        """
        Subtract baseline variance from actual variance.

        Args:
            in_actual_variance (PathLike): Path to image containing actual variance.
            in_baseline_variance (PathLike): Path to image containing baseline variance.
            out_corrected_variance (PathLike): Where to write corrected variance image.
            name (str): What to name the task.

        Returns:
            Dict: Dummy pydoit dict.
        """
        out_parent_dir = Path(out_corrected_variance).parent
        return dict(
            name=name,
            file_dep=[in_actual_variance, in_baseline_variance],
            targets=[out_corrected_variance],
            actions=[f"mkdir -p '{out_parent_dir}'",
                    f"3dcalc -float -a '{in_actual_variance}' -b '{in_baseline_variance}' -expr 'a-b' -prefix '{out_corrected_variance}'"]
        )

    # Continuous alpha.
    variables = "amplitude SNR".split()
    analysis = "alpha"
    for start_volume in START_VOLUMES:
        for variable in variables:
            yield create_task(
                in_actual_variance=fname.variance_whole_brain(start_volume=start_volume, variable=variable, analysis=analysis),
                in_baseline_variance=fname.variance_whole_brain(start_volume="na", variable="na", analysis="baseline"),
                out_corrected_variance=fname.variance_whole_brain_baselined(start_volume=start_volume, variable=variable, analysis=analysis),
                name=f"analysis--{analysis}, variable--{variable}, start_volume--{start_volume}",
            )

    # Continuous ssVEP.
    analysis = "ssvep"
    for start_volume in EXPANDED_START_VOLUMES:
        for variable in variables:
            yield create_task(
                in_actual_variance=fname.variance_whole_brain(start_volume=start_volume, variable=variable, analysis=analysis),
                in_baseline_variance=fname.variance_whole_brain(start_volume="na", variable="na", analysis="baseline"),
                out_corrected_variance=fname.variance_whole_brain_baselined(start_volume=start_volume, variable=variable, analysis=analysis),
                name=f"analysis--{analysis}, variable--{variable}, start_volume--{start_volume}",
            )


# Permutation testing.
def task_scramble_data() -> Dict:
    """
    Scramble our ssVEP and alpha time series so we can use them for permutation thresholding.
    """
    def create_task(in_series: PathLike, out_series: PathLike, name: str) -> Dict:
        """
        Make this task generalizable.

        Args:
            in_series (PathLike): Path to series to scramble.
            out_series (PathLike): Where to write scrambled series.
            name (str): Name of the sub-task.

        Returns:
            Dict: Placeholder dict used by pydoit.
        """
        sources = dict(in_series=in_series)
        targets = dict(out_series=out_series)
        kwargs = {**sources, **targets}

        return dict(
            name=name,
            actions=[(scramble_series.main, [], kwargs)],
            file_dep=list(sources.values()),
            targets=list(targets.values()),
        )

    for subject in SUBJECTS:
        for permutation in PERMUTATIONS:
            analysis = "alpha"
            start_volume = 4
            old_and_new_variables = {"values": "amplitude", "SNRs": "SNR"}
            for old_variable, variable in old_and_new_variables.items():
                yield create_task(
                    in_series=fname.eeg_alpha(subject=subject, data=old_variable),
                    out_series=fname.scrambled_series(subject=subject, start_volume=start_volume, variable=variable, analysis=analysis, permutation=permutation),
                    name=f"subject--{subject}, start_volume--{start_volume}, variable--{variable}, analysis--{analysis}, permutation--{permutation}",
                )
            start_volume = "na"
            yield create_task(
                in_series=fname.eeg_trial_alpha(subject=subject, data="values"),
                out_series=fname.scrambled_series(subject=subject, start_volume=start_volume, variable=variable, analysis=analysis, permutation=permutation),
                name=f"subject--{subject}, start_volume--{start_volume}, variable--{variable}, analysis--{analysis}, permutation--{permutation}",
            )

            analysis = "ssvep"
            start_volume = 5
            old_and_new_variables = {"amplitudes": "amplitude", "SNRs": "SNR"}
            for old_variable, variable in old_and_new_variables.items():
                yield create_task(
                    in_series=fname.eeg_sliding_sliding_window_improved(subject=subject, variable="amplitudes"),
                    out_series=fname.scrambled_series(subject=subject, start_volume=start_volume, variable=variable, analysis=analysis, permutation=permutation),
                    name=f"subject--{subject}, start_volume--{start_volume}, variable--{variable}, analysis--{analysis}, permutation--{permutation}",
                )
            start_volume = "na"
            yield create_task(
                in_series=fname.freqtag_better_sliding_window_channel(subject=subject, channel=20, variable="trialpow"),
                out_series=fname.scrambled_series(subject=subject, start_volume=start_volume, variable=variable, analysis=analysis, permutation=permutation),
                name=f"subject--{subject}, start_volume--{start_volume}, variable--{variable}, analysis--{analysis}, permutation--{permutation}",
            )
def task_get_maxes_and_mins() -> Dict:
    """
    Get distributions of our maxes and mins for our permutations.

    Returns:
        Dict: Dummy pydoit dict.
    """
    def create_task(in_permutations: List[PathLike], out_maxes: PathLike, out_mins: PathLike, name: str) -> Dict:
        """
        Get distributions of our maxes and mins for our permutations.

        Args:
            in_permutations (List[PathLike]): List of paths to permutation images.
            out_maxes (PathLike): Table of max values from those images.
            out_mins (PathLike): Table of min values from those images.
            name (str): What to name the task.

        Returns:
            Dict: Dummy pydoit dict.
        """
        sources = [*in_permutations]
        targets = dict(out_maxes=out_maxes, out_mins=out_mins)
        kwargs = {**targets, "in_permutations": in_permutations}

        return dict(
            name=name,
            actions=[(get_maxes_and_mins.main, [], kwargs)],
            file_dep=sources,
            targets=list(targets.values()),
        )

    variable = "amplitude"

    start_volume_dict = {"alpha": 4, "ssvep": 5}
    for analysis, start_volume in start_volume_dict.items():
        # Get non-baseline-corrected.
        yield create_task(
            in_permutations=[fname.correlations_whole_brain_permutations_ttest(start_volume=start_volume, variable=variable, analysis=analysis, permutation=permutation) for permutation in PERMUTATIONS],
            out_maxes=fname.maxes_mins_table(start_volume=start_volume, variable=variable, analysis=analysis, outfile="maxes"),
            out_mins=fname.maxes_mins_table(start_volume=start_volume, variable=variable, analysis=analysis, outfile="mins"),
            name=f"non-baselined, start_volume--{start_volume}, variable--{variable}, analysis--{analysis}",
        )
        # Get baseline corrected.
        yield create_task(
            in_permutations=[fname.compared_permutations_to_canonical(start_volume=start_volume, variable=variable, analysis=analysis, permutation=permutation) for permutation in PERMUTATIONS],
            out_maxes=fname.maxes_mins_table(start_volume=start_volume, variable=variable, analysis=analysis, outfile="maxes"),
            out_mins=fname.maxes_mins_table(start_volume=start_volume, variable=variable, analysis=analysis, outfile="mins"),
            name=f"baselined, start_volume--{start_volume}, variable--{variable}, analysis--{analysis}",
        )
def task_plot_maxes_and_mins_distributions() -> Dict:
    """
    Plot the distributions of our maxes and mins.
    """
    def create_task(in_distribution: PathLike, out_plot: PathLike, title: str, xlim: Tuple[float, float], ylim: Tuple[float, float], name: str) -> Dict:
        """
        Make task to plot the distributions of our maxes and mins.
        """
        sources = dict(in_distribution=in_distribution)
        targets = dict(out_plot=out_plot)
        kwargs = {**targets, **sources, "title": title, "xlim": xlim, "ylim": ylim}

        return dict(
            name=name,
            actions=[(plot_distribution.main, [], kwargs)],
            file_dep=list(sources.values()),
            targets=list(targets.values()),
        )

    xlims = (-1, 101)
    variable = "amplitude"
    ylims = {"maxes": (0, .1), "mins": (-.1, 0)}
    start_volumes = {"alpha": 4, "ssvep": 5}

    for outfile, ylim in ylims.items():
        for analysis, start_volume in start_volumes.items():
            yield create_task(
                in_distribution=fname.maxes_mins_table(start_volume=start_volume, variable=variable, analysis=analysis, outfile=outfile),
                out_plot=fname.maxes_mins_plot(start_volume=start_volume, variable=variable, analysis=analysis, outfile=outfile),
                title=f"Distribution of {analysis} {outfile}",
                xlim=xlims,
                ylim=ylim,
                name=f"analysis--{analysis}, outfile--{outfile}",
            )
def task_threshold_results() -> Dict:
    """
    Read the results of our average averages. Rank them. Output the percentile the actual correlation is in.

    Returns:
        Dict: Dummy pydoit dict.
    """
    def create_task(in_permutation_results: List[PathLike], in_actual_result: PathLike, out_table: PathLike, name: str) -> Dict:
        """
        Create a task to test an actual result against permutation results.

        Args:
            in_permutation_results (List[PathLike]): List of paths to the permutation results.
            in_actual_result (PathLike): Path to our actual result.
            out_table (PathLike): Where to write a table to help us check our work.
            out_result (PathLike): Where to write the results of this task. Should include percentile and rank.
            name (str): What to name the subtask.

        Returns:
            Dict: Dummy pydoit dict. Ignore.
        """
        sources = [*in_permutation_results, in_actual_result]
        targets = dict(out_table=out_table)
        kwargs = {**targets, "in_permutation_results": in_permutation_results, "in_actual_result": in_actual_result}

        return dict(
            name=name,
            actions=[(test_permutations.main, [], kwargs)],
            file_dep=sources,
            targets=list(targets.values()),
        )

    for mask in "calcarine occipital".split():
        start_volume = 4
        variable = "amplitude"
        analysis = "alpha"
        yield create_task(
            in_permutation_results=[fname.average_averages_correlation_permutation(permutation=permutation, variable=variable, start_volume=start_volume, mask=mask, analysis=analysis, outfile="results") for permutation in PERMUTATIONS],
            in_actual_result=fname.average_averages_correlation(variable=variable, start_volume=start_volume, mask=mask, analysis=analysis, outfile="results"),
            out_table=fname.threshold_outfile(variable=variable, start_volume=start_volume, mask=mask, analysis=analysis, outfile="table"),
            name=f"startvolume--{start_volume}, mask--{mask}, variable--{variable}, analysis--{analysis}",
        )

        start_volume = 5
        analysis = "ssvep"
        yield create_task(
            in_permutation_results=[fname.average_averages_correlation_permutation(permutation=permutation, variable=variable, start_volume=start_volume, mask=mask, analysis=analysis, outfile="results") for permutation in PERMUTATIONS],
            in_actual_result=fname.average_averages_correlation(variable=variable, start_volume=start_volume, mask=mask, analysis=analysis, outfile="results"),
            out_table=fname.threshold_outfile(variable=variable, start_volume=start_volume, mask=mask, analysis=analysis, outfile="table"),
            name=f"startvolume--{start_volume}, mask--{mask}, variable--{variable}, analysis--{analysis}",
        )
def task_cat_permutations() -> Dict:
    """
    Concatenate our permutation correlations together.
    """
    def create_task(images: List[PathLike], out_prefix: PathLike, name: str) -> dict:
        """
        Use 3dTcat to concatenate permutation correlations together.

        Args:
            images (List[PathLike]): List of paths to our permutation correlations.
            out_prefix (PathLike): Where to write our concatenated permutations to.
            name (str): What to name the task.

        Returns:
            dict: Dummy dict for pydoit.
        """
        sources = dict(
            images=images,
        )

        targets = dict(
            out_prefix=out_prefix,
        )

        kwargs = dict(
            images=[f"{path}[0]" for path in sources["images"]],
            out_prefix=out_prefix
        )

        return dict(
            name=name,
            actions=[(cat_images.main, [], kwargs)],
            file_dep=list(sources.values())[0],
            targets=list(targets.values()),
        )

    start_volume = 4
    variable = "amplitude"
    analysis = "alpha"
    yield create_task(
        images=[fname.correlations_whole_brain_permutations_ttest(start_volume=start_volume, variable=variable, analysis=analysis, permutation=permutation) for permutation in PERMUTATIONS],
        out_prefix=fname.catenated_permutations(start_volume=start_volume, variable=variable, analysis=analysis),
        name=f"startvolume--{start_volume}, variable--{variable}, analysis--{analysis}",
    )

    start_volume = 5
    analysis = "ssvep"
    yield create_task(
        images=[fname.correlations_whole_brain_permutations_ttest(start_volume=start_volume, variable=variable, analysis=analysis, permutation=permutation) for permutation in PERMUTATIONS],
        out_prefix=fname.catenated_permutations(start_volume=start_volume, variable=variable, analysis=analysis),
        name=f"startvolume--{start_volume}, variable--{variable}, analysis--{analysis}",
    )
def task_calc_percentiles() -> Dict:
    """
    Calculate percentiles for each voxel based on the permutations we have.
    """
    def create_task(catenated_permutations: PathLike, actual_correlation: PathLike, out_percentiles: PathLike, name: str) -> Dict:
        """
        Create task to calculate a percentile for each voxel based on the permutations we have.

        Args:
            catenated_permutations (PathLike): Path to an image of catenated permutation correlations.
            actual_correlation (PathLike): Path to an image of actual correlations.
            out_percentiles (PathLike): Where to write our percentiles image to.
            name (str): What to name the task.

        Returns:
            Dict: Dummy pydoit dict.
        """
        sources = dict(catenated_permutations=catenated_permutations, actual_correlation=actual_correlation)
        targets = dict(out_percentiles=out_percentiles)
        kwargs = {**targets, **sources}

        return dict(
            name=name,
            actions=[(calc_percentiles.main, [], kwargs)],
            file_dep=list(sources.values()),
            targets=list(targets.values()),
        )

    start_volume = 4
    variable = "amplitude"
    analysis = "alpha"
    yield create_task(
        catenated_permutations=fname.catenated_permutations(start_volume=start_volume, variable=variable, analysis=analysis),
        actual_correlation=fname.correlations_whole_brain_alpha_ttest(start_volume=start_volume, data="values"),
        out_percentiles=fname.whole_brain_percentiles(start_volume=start_volume, variable=variable, analysis=analysis),
        name=f"startvolume--{start_volume}, variable--{variable}, analysis--{analysis}",
    )

    start_volume = 5
    analysis = "ssvep"
    yield create_task(
        catenated_permutations=fname.catenated_permutations(start_volume=start_volume, variable=variable, analysis=analysis),
        actual_correlation=fname.correlations_improved_whole_brain_ttest(start_volume=start_volume, variable="amplitudes"),
        out_percentiles=fname.whole_brain_percentiles(start_volume=start_volume, variable=variable, analysis=analysis),
        name=f"startvolume--{start_volume}, variable--{variable}, analysis--{analysis}",
    )


# Helper functions.
def _print_paths(paths: Iterable, name: str = None) -> None:
    """
    Prints an iterable and whether its paths exist.

    Useful for debugging.
    """
    if name:
        print(name)
    for path in paths:
        path = Path(path)
        exists = path.exists()
        exists_str = "exists" if exists else "doesn't exist"
        print(f"{path}  :  {exists_str}")
    print()
def get_prefix(filename: str) -> str:
    """
    Returns the prefix of an AFNI file. (Everything before the final "+".)
    """
    return "".join(filename.split("+")[:-1])
def get_matlab_prefix(filename: str) -> str:
    """
    Returns the prefix of MatLab files, which sometimes have multiple "." characters in their filenames.
    """
    return ".".join(filename.split(".")[:-2])
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
def encode_dict_in_bash(dictionary: dict) -> str:
    """
    Make dict into a series of key value pairs readable on command line.
    """
    pairs = []
    for key, value in dictionary.items():
        pairs.append(f"{key}={value}")

    return " ".join(pairs)
