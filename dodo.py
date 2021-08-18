#!/usr/bin/env python3
"""
Do-it script to execute the entire pipeline using the doit tool:
http://pydoit.org

All the filenames are defined in config.py
"""
# Import external modules and libraries.
from pathlib import Path
from typing import Dict, Iterable

# Import internal modules and libraries.
from config import fname, SUBJECTS, n_jobs, COMPONENTS_TO_REMOVE

# Import tasks
import create_bids_root
import bidsify_subject
import afniproc
import matlab
import trim_func_images
import align
import resample
import smooth
import scale
import deconvolve
import ttest_2_groups
import write_json
import correlate_eeg_fmri
import ttest
import average_freqtags
import combine_masks
import apply_mask
import clusterize


# Configuration for the pydoit tool.
DOIT_CONFIG = dict(
    # While running scripts, output everything the script is printing to the screen.
    verbosity=2,

    # Tell doit to use all processors I've said this machine has.
    num_process=n_jobs,

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
            name=subject,
            actions=[(bidsify_subject.main, [sources, targets])],
            file_dep=list(sources.values()),
            targets=list(targets.values()),
        )


# fMRI-only tasks.
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
            name=subject,
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
            name=subject,
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
            name=subject,
            actions=[(trim_func_images.main, [], kwargs)],
            file_dep=sources,
            targets=targets,
        )
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
            name=subject,
            actions=[(resample.main, [], kwargs)],
            file_dep=sources,
            targets=targets,
        )
def task_ttest_deconvolutions() -> Dict:
    """
    Here we t-test our 3dDeconvolve results. Onward!
    """
    for subbrick in range(1, 9):
        sources = [fname.afniproc_deconvolved_resampled(subject=subject) for subject in SUBJECTS]
        targets = [fname.afniproc_ttest_result(subbrick=subbrick)]

        kwargs = dict(
            images=[fname.afniproc_deconvolved_resampled(subject=subject) + f"[{subbrick}]" for subject in SUBJECTS],
            prefix=get_prefix(fname.afniproc_ttest_result(subbrick=subbrick)),
        )
        yield dict(
            name=f"subbrick-{subbrick}",
            actions=[(ttest.main, [], kwargs)],
            file_dep=sources,
            targets=targets,
        )


# EEG-only tasks.
def task_prepare_to_convert_eeg() -> Dict:
    """
    Write a JSON file that will be read later by a MatLab script I wrote to convert EEG files.
    """
    sources = [fname.brainvision_eeg(subject=subject) for subject in SUBJECTS]
    targets = [fname.converteeg_json]

    data = []
    for subject in SUBJECTS:
        data.append(dict(
            brainvision_dir=fname.brainvision_dir,
            brainvision_name=Path(fname.brainvision_eeg(subject=subject)).name,
            converted_path=fname.converted_eeg(subject=subject),
            setname=f"sub-{subject}",
        ))

    kwargs = dict(
        out_path=fname.converteeg_json,
        data=data,
    )

    return dict(
        actions=[(write_json.main, [], kwargs)],
        file_dep=sources,
        targets=targets,
    )
def task_convert_eeg() -> Dict:
    """
    Convert BrainVision EEG files into EEGLAB EEG files, which are easier to edit.

    Because MatLab is ~quirky~, we can't do multithreading on this one.
    We must convert all subjects in one glorious blaze!
    Sigh.
    """
    sources = [fname.converteeg_json]
    targets = [fname.converted_eeg(subject=subject) for subject in SUBJECTS]

    kwargs = dict(
        path_to_script="convert_eegs.m",
    )

    return dict(
        actions=[(matlab.main, [], kwargs)],
        file_dep=sources,
        targets=targets,
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

    kwargs = dict(
        path_to_script="preprocess_eegs.m",
    )

    return dict(
        actions=[(matlab.main, [], kwargs)],
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

    Because MatLab is ~quirky~, we can't do multithreading on this one.
    We must do all subjects in one glorious blaze!
    Sigh.
    """
    sources = [fname.segmenteeg_json]
    targets = [fname.segmented_eeg(subject=subject) for subject in SUBJECTS]

    kwargs = dict(
        path_to_script="segment_eegs.m",
    )

    return dict(
        actions=[(matlab.main, [], kwargs)],
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
            time_delta=-8,
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
    Trim our preprocessed EEG files to the time at which the fMRI turned on.

    Because MatLab is ~quirky~, we can't do multithreading on this one.
    We must do all subjects in one glorious blaze!
    Sigh.
    """
    sources = [fname.trimeeg_json]
    targets = [fname.trimmed_eeg(subject=subject) for subject in SUBJECTS]

    kwargs = dict(
        path_to_script="trim_eegs.m",
    )

    return dict(
        actions=[(matlab.main, [], kwargs)],
        file_dep=sources,
        targets=targets,
    )
def task_prepare_to_moving_moving_window_eeg() -> Dict:
    """
    Write a JSON file that will be read later by a MatLab script we wrote to run
    a moving moving window analysis. Welcome to the cutting edge!
    """
    sources = [fname.trimmed_eeg(subject=subject) for subject in SUBJECTS]
    targets = [fname.movingmovingwindoweeg_json]

    data = []
    for subject in SUBJECTS:
        out_path = Path(fname.moving_moving_windowed_eeg(subject=subject))
        data.append(dict(
            in_filename=Path(fname.trimmed_eeg(subject=subject)).name,
            in_dir=fname.trimeeg_dir,
            out_stem=str(out_path.parent / out_path.name.split(".")[0]),
            out_tsv_name=fname.out_tsv_name(subject=subject),
        ))

    kwargs = dict(
        out_path=fname.movingmovingwindoweeg_json,
        data=data,
    )

    return dict(
        actions=[(write_json.main, [], kwargs)],
        file_dep=sources,
        targets=targets,
    )
def task_moving_moving_window_eeg() -> Dict:
    """
    Run a moving moving window analysis on our trimmed, preprocessed EEG files.

    Because MatLab is ~quirky~, we can't do multithreading on this one.
    We must do all subjects in one glorious blaze!
    Sigh.
    """
    sources = [fname.movingmovingwindoweeg_json]
    targets = [fname.moving_moving_windowed_eeg(subject=subject) for subject in SUBJECTS]

    kwargs = dict(
        path_to_script="moving_moving_window_eegs.m",
    )

    return dict(
        actions=[(matlab.main, [], kwargs)],
        file_dep=sources,
        targets=targets,
    )
def task_prepare_to_freqtag_eeg() -> Dict:
    """
    Write a JSON file that will be read later by a MatLab script we wrote to run a frequency tagging analysis.
    """
    sources = [fname.segmented_eeg(subject=subject) for subject in SUBJECTS]
    targets = [fname.freqtageeg_json]

    data = []
    for subject in SUBJECTS:
        data.append(dict(
            in_eeg_name=Path(fname.segmented_eeg(subject=subject)).name,
            in_eeg_dir=fname.segmenteeg_dir,
            out_fft_path=fname.out_fft_path(subject=subject),
            out_hilbert_path=fname.out_hilbert_path(subject=subject),
        ))

    kwargs = dict(
        out_path=fname.freqtageeg_json,
        data=data,
    )

    return dict(
        actions=[(write_json.main, [], kwargs)],
        file_dep=sources,
        targets=targets,
    )
def task_freqtag_eeg() -> Dict:
    """
    Run a frequency tagging analysis!

    Because MatLab is ~quirky~, we can't do multithreading on this one.
    We must do all subjects in one glorious blaze!
    """
    sources = [fname.freqtageeg_json]
    targets = [fname.out_fft_path(subject=subject) for subject in SUBJECTS]
    targets += [fname.out_hilbert_path(subject=subject) for subject in SUBJECTS]

    kwargs = dict(
        path_to_script="freqtag_pipeline.m",
    )

    return dict(
        actions=[(matlab.main, [], kwargs)],
        file_dep=sources,
        targets=targets,
    )
def task_mean_mean_fft() -> Dict:
    """
    Average the mean FFT calculated for each subject. But this time, across ALL subjects!
    """
    sources = [fname.out_fft_path(subject=subject) for subject in SUBJECTS]
    targets = [fname.mean_mean_fft]

    kwargs = dict(
        in_tsvs=[fname.out_fft_path(subject=subject) for subject in SUBJECTS],
        out_tsv=fname.mean_mean_fft,
    )

    return dict(
        actions=[(average_freqtags.main, [], kwargs)],
        file_dep=sources,
        targets=targets,
    )
def task_mean_mean_hilbert() -> Dict:
    """
    Average the mean hilbert calculated for each subject. But this time, across ALL subjects!
    """
    sources = [fname.out_hilbert_path(subject=subject) for subject in SUBJECTS]
    targets = [fname.mean_mean_hilbert]

    kwargs = dict(
        in_tsvs=[fname.out_hilbert_path(subject=subject) for subject in SUBJECTS],
        out_tsv=fname.mean_mean_hilbert,
    )

    return dict(
        actions=[(average_freqtags.main, [], kwargs)],
        file_dep=sources,
        targets=targets,
    )


# fMRI/EEG correlation tasks.
def task_trim_func_images_again() -> Dict:
    """
    Trim func images one last time so we can correlate the BOLD response with EEG signal.
    """
    for subject in SUBJECTS:

        sources = [fname.trimmed_func(subject=subject)]

        for volumes_to_remove in range(1, 4):

            targets = [
                fname.final_func(subject=subject, start_volume=volumes_to_remove),
            ]

            kwargs = dict(
                new_start_volume=volumes_to_remove,
                func_path=fname.trimmed_func(subject=subject),
                out_prefix=get_prefix(fname.final_func(subject=subject, start_volume=volumes_to_remove)),
            )

            yield dict(
                name=f"sub-{subject}_startvolume-{volumes_to_remove}",
                actions=[(trim_func_images.main2, [], kwargs)],
                file_dep=sources,
                targets=targets,
            )
def task_correlate_eeg_fmri() -> Dict:
    """
    This is it! Huzzah! Correlate our EEG and fMRI data across the whole brain.
    """
    for subject in SUBJECTS:
        for start_volume in range(1, 4):

            sources = [
                fname.final_func(subject=subject, start_volume=start_volume),
                fname.out_tsv_name(subject=subject),
            ]
            targets = [
                fname.correlation_image(subject=subject, start_volume=start_volume),
            ]

            kwargs = dict(
                in_image_path=fname.final_func(subject=subject, start_volume=start_volume),
                in_eeg_path=fname.out_tsv_name(subject=subject),
                out_image_path=fname.correlation_image(subject=subject, start_volume=start_volume),
            )

            yield dict(
                name=f"sub-{subject}_startvolume-{start_volume}",
                actions=[(correlate_eeg_fmri.main, [], kwargs)],
                file_dep=sources,
                targets=targets,
            )
def task_ttest_eeg_fmri_correlations() -> Dict:
    """
    ttest the correlations we calculated.
    """
    for start_volume in range(1, 4):
        sources = [fname.correlation_image(subject=subject, start_volume=start_volume) for subject in SUBJECTS]
        targets = [fname.correlations_ttest(start_volume=start_volume)]

        kwargs = dict(
            images=[fname.correlation_image(subject=subject, start_volume=start_volume) for subject in SUBJECTS],
            prefix=get_prefix(fname.correlations_ttest(start_volume=start_volume)),
        )

        yield dict(
            name=f"startvolume-{start_volume}",
            actions=[(ttest.main, [], kwargs)],
            file_dep=sources,
            targets=targets,
        )
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
            name=mask["name"],
            actions=[(resample.main, [], kwargs)],
            file_dep=sources,
            targets=targets,
        )
def task_apply_masks_to_correlations() -> Dict:
    """
    Apply masks to our correlation results.
    """
    masks = [
        dict(name="calcarine", in_path=fname.resampled_mask(mask="calcarine")),
        dict(name="occipital", in_path=fname.resampled_mask(mask="occipital")),
    ]

    for start_volume in range(1, 4):
        for mask in masks:
            sources = [fname.correlations_ttest(start_volume=start_volume)]
            sources += [mask["in_path"]]

            targets = [fname.correlations_ttest_masked(start_volume=start_volume, mask=mask["name"])]

            kwargs = dict(
                in_image=fname.correlations_ttest(start_volume=start_volume),
                in_mask=mask["in_path"],
                out_prefix=get_prefix(fname.correlations_ttest_masked(start_volume=start_volume, mask=mask["name"])),
            )

            yield dict(
                name=f"startvolume-{start_volume}_mask-{mask['name']}",
                actions=[(apply_mask.main, [], kwargs)],
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
                name=f"sub {subject}, mask {mask['name']}",
                actions=[(apply_mask.main, [], kwargs)],
                file_dep=sources,
                targets=targets,
            )
def task_clusterize_irfs() -> Dict:
    """
    FOR each subject, FOR each masked IRF, GET all clusters.
    """
    for subject in SUBJECTS:
        masked_irfs = [
            dict(name="calcarine", in_path=fname.masked_irf(subject=subject, mask="calcarine")),
            dict(name="occipital", in_path=fname.masked_irf(subject=subject, mask="occipital")),
        ]

        for irf in masked_irfs:
            sources = [fname.masked_irf(subject=subject, mask=irf["name"])]
            targets = [fname.clusters(subject=subject, mask=irf["name"])]

            kwargs = dict(
                in_image=fname.masked_irf(subject=subject, mask=irf["name"]),
                out_prefix=get_prefix(fname.clusters(subject=subject, mask=irf["name"])),
                out_summary=fname.clusters_summary(subject=subject, mask=irf["name"]),
                subbrick=3,
            )

            yield dict(
                name=f"subject - {subject}, IRF - {irf['name']}",
                actions=[(clusterize.main, [], kwargs)],
                file_dep=sources,
                targets=targets,
            )


# Tasks to test afniproc vs fmriprep.
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
            fwhm="4.0",
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
def task_deconvolve_fmriprep() -> Dict:
    """
    Deconvolve our scaled fMRIPrep images.
    """
    for subject in SUBJECTS:
        sources = [
            fname.fmriprep_scaled(subject=subject),
        ]
        targets = [
            fname.fmriprep_deconvolved(subject=subject),
            fname.fmriprep_irf(subject=subject),
        ]

        kwargs = dict(
            regressors_tsv=fname.regressors_tsv(subject=subject),
            regressors_dir=fname.fmriprep_regressors_dir(subject=subject),
            deconvolved_prefix=get_prefix(fname.fmriprep_deconvolved(subject=subject)),
            IRF_prefix=get_prefix(fname.fmriprep_irf(subject=subject)),
            func_path=fname.fmriprep_scaled(subject=subject),
            events_tsv=fname.bids_events(subject=subject),
            events_dir=fname.fmriprep_deconvolve_dir(subject=subject)
        )

        yield dict(
            name=subject,
            actions=[(deconvolve.main, [], kwargs)],
            file_dep=sources,
            targets=targets,
        )
def task_align_fmriprep_irfs() -> Dict:
    """
    Align our fMRIPrep IRFs to the space of the Kastner cortex masks so we may compare them with afni_proc's IRFs.

    Note that the aligned IRFs appear in the fMRIPrep IRF dir, not the alignment template dir.
    I'm aware that this is quirky, but it's convenient.

    Also, we don't separate each subject into a sub-task for a very good reason. It's more efficient to calculate the alignment
    for one subject, then apply the transformation to all the others.
    """
    sources = [fname.fmriprep_irf(subject=subject) for subject in SUBJECTS]
    sources += [
        fname.atlas_template,
        fname.fmriprep_template,
    ]

    targets = [fname.fmriprep_aligned_irf(subject=subject) for subject in SUBJECTS]

    kwargs = dict(
        from_images=[fname.fmriprep_irf(subject=subject) for subject in SUBJECTS],
        from_template=fname.fmriprep_template,
        to_template=fname.atlas_template,
        to_dir=fname.fmriprep_alignment_dir,
    )

    return dict(
        actions=[(align.main, [], kwargs)],
        file_dep=sources,
        targets=targets,
    )
def task_resample_fmriprep_irfs() -> Dict:
    """
    Resample our IRFs to the space of the Kastner cortex masks so we may compare them with our afniproc IRFs.
    """
    for subject in SUBJECTS:
        sources = [
            fname.resampled_template,
            fname.fmriprep_aligned_irf(subject=subject)
        ]
        targets = [
            fname.fmriprep_resampled_irf(subject=subject)
        ]

        kwargs = dict(
            from_image=fname.fmriprep_aligned_irf(subject=subject),
            to_image=fname.resampled_template,
            to_prefix=get_prefix(fname.fmriprep_resampled_irf(subject=subject)),
        )

        yield dict(
            name=subject,
            actions=[(resample.main, [], kwargs)],
            file_dep=sources,
            targets=targets,
        )
def task_ttest_fmriprep_vs_afniproc() -> Dict:
    """
    ttest IRFs of fMRIPrep vs afniproc so we may find the TRUE pipeline of kings.
    """
    sources = [fname.fmriprep_resampled_irf(subject=subject) for subject in SUBJECTS]
    sources += [fname.afniproc_resampled_irf(subject=subject) for subject in SUBJECTS]

    for subbrick in range(1, 9):
        targets = [fname.ttest_result(subbrick=subbrick)]

        kwargs = dict(
            list_a=[fname.fmriprep_resampled_irf(subject=subject) + f"[{subbrick}]" for subject in SUBJECTS],
            list_b=[fname.afniproc_resampled_irf(subject=subject) + f"[{subbrick}]" for subject in SUBJECTS],
            prefix=get_prefix(fname.ttest_result(subbrick=subbrick)),
        )

        yield dict(
            name=f"subbrick {subbrick}",
            actions=[(ttest_2_groups.main, [], kwargs)],
            file_dep=sources,
            targets=targets,
        )


# Helper functions.
def _print_paths(paths: Iterable, name: str=None) -> None:
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
