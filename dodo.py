#!/usr/bin/env python3
"""
Do-it script to execute the entire pipeline using the doit tool:
http://pydoit.org

All the filenames are defined in config.py
"""
# Import external modules and libraries.
from pathlib import Path
from typing import Dict, Iterable
import textwrap

# Import internal modules and libraries.
from config import FREQUENCIES, fname, SUBJECTS, n_jobs, COMPONENTS_TO_REMOVE, START_VOLUMES

# Import tasks
import create_bids_root
import bidsify_subject
import afniproc
import trim_func_images
import align
import resample
import write_json
import correlate_eeg_fmri
import ttest
import combine_masks
import apply_mask
import clusterize
import average_voxels


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
            name=f"subject--{subject}",
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
    for subbrick in range(1, 9):
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

    path_to_script="convert_eegs.m"
    action = f"""
        python3 matlab.py
        --run_script_at {path_to_script}
    """.split()

    return dict(
        actions=[action],
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

    path_to_script="preprocess_eegs.m"
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

    path_to_script="segment_eegs.m"
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

    path_to_script="trim_eegs.m"
    action = f"""
        python3 matlab.py
        --run_script_at {path_to_script}
    """.split()

    return dict(
        actions=[action],
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
        for frequency in FREQUENCIES:
            out_path = Path(fname.moving_moving_windowed_eeg(subject=subject, frequency=frequency))
            data.append(dict(
                in_filename=Path(fname.trimmed_eeg(subject=subject)).name,
                in_dir=fname.trimeeg_dir,
                out_stem=str(out_path.parent / out_path.name.split(".")[0]),
                out_tsv_name=fname.out_tsv_name(subject=subject, frequency=frequency),
                frequency=frequency,
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
    for frequency in FREQUENCIES:

        sources = [fname.movingmovingwindoweeg_json]
        targets = [fname.moving_moving_windowed_eeg(subject=subject, frequency=frequency) for subject in SUBJECTS]

        path_to_script="moving_moving_window_eegs.m"
        action = f"""
            python3 matlab.py
            --run_script_at {path_to_script}
        """.split()

        yield dict(
            name=f"frequency--{frequency}",
            actions=[action],
            file_dep=sources,
            targets=targets,
        )
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
            dat=fname.bids_dat(subject=subject)
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

                %% Run functions.
                durations = get_durations('{sources["dat"]}')
                trialpow = [];
                winmat3d = [];
                phasestabmat = [];
                trialSNR = [];
                for i = 1:numel(dataset(1,1,:))
                    duration = durations(i)
                    frequency = 1./(duration/50)
                    [minitrialpow,miniwinmat3d,miniphasestabmat,minitrialSNR] = freqtag_slidewin(dataset(:,:,i), 0, stimulus_start:stimulus_end, stimulus_start:stimulus_end, frequency, 600, sampling_rate, 'TEMP');
                    trialpow = [trialpow, minitrialpow];
                    winmat3d = [winmat3d, miniwinmat3d];
                    phasestabmat = [phasestabmat, miniphasestabmat];
                    trialSNR = [trialSNR, minitrialSNR];
                end
                meanwinmat = mean(winmat3d, 3);

                %% Save output variables.
                save('{targets["trialpow"]}', 'trialpow');
                save('{targets["winmat3d"]}', 'winmat3d');
                save('{targets["phasestabmat"]}', 'phasestabmat');
                save('{targets["trialSNR"]}', 'trialSNR');
                save('{targets["meanwinmat"]}', 'meanwinmat');
            end

            function [durations] = get_durations(dat_path)
                datmat = importdata(dat_path)
                durations = datmat(1:end, 5);
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
def task_trim_func_images_again() -> Dict:
    """
    Trim func images one last time so we can correlate the BOLD response with EEG signal.
    """
    for subject in SUBJECTS:

        sources = [fname.trimmed_func(subject=subject)]

        for volumes_to_remove in START_VOLUMES:

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
def task_correlate_eeg_fmri() -> Dict:
    """
    This is it! Huzzah! Correlate our EEG and fMRI data across the whole brain.
    """
    for frequency in FREQUENCIES:
        for subject in SUBJECTS:
            for start_volume in START_VOLUMES:

                sources = [
                    fname.final_func(subject=subject, start_volume=start_volume),
                    fname.out_tsv_name(subject=subject, frequency=frequency),
                ]
                targets = [
                    fname.correlation_image(subject=subject, start_volume=start_volume, frequency=frequency),
                ]

                kwargs = dict(
                    in_image_path=fname.final_func(subject=subject, start_volume=start_volume),
                    in_eeg_path=fname.out_tsv_name(subject=subject, frequency=frequency),
                    out_image_path=fname.correlation_image(subject=subject, start_volume=start_volume, frequency=frequency),
                )

                yield dict(
                    name=f"sub--{subject}, startvolume--{start_volume}, frequency--{frequency}",
                    actions=[(correlate_eeg_fmri.main, [], kwargs)],
                    file_dep=sources,
                    targets=targets,
                )
def task_ttest_eeg_fmri_correlations() -> Dict:
    """
    ttest the correlations we calculated.
    """
    for frequency in FREQUENCIES:
        for start_volume in START_VOLUMES:
            sources = [fname.correlation_image(subject=subject, start_volume=start_volume, frequency=frequency) for subject in SUBJECTS]
            targets = [fname.correlations_ttest(start_volume=start_volume, frequency=frequency)]

            kwargs = dict(
                images=[fname.correlation_image(subject=subject, start_volume=start_volume, frequency=frequency) for subject in SUBJECTS],
                prefix=get_prefix(fname.correlations_ttest(start_volume=start_volume, frequency=frequency)),
            )

            yield dict(
                name=f"startvolume--{start_volume}, frequency--{frequency}",
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
            name=f"mask--{mask['name']}",
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
    
    for frequency in FREQUENCIES:
        for start_volume in START_VOLUMES:
            for mask in masks:
                sources = [fname.correlations_ttest(start_volume=start_volume, frequency=frequency)]
                sources += [mask["in_path"]]

                targets = [fname.correlations_ttest_masked(start_volume=start_volume, frequency=frequency, mask=mask["name"])]

                kwargs = dict(
                    in_image=fname.correlations_ttest(start_volume=start_volume, frequency=frequency),
                    in_mask=mask["in_path"],
                    out_prefix=get_prefix(fname.correlations_ttest_masked(start_volume=start_volume, frequency=frequency, mask=mask["name"])),
                )

                yield dict(
                    name=f"startvolume--{start_volume}, mask--{mask['name']}, frequency--{frequency}",
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
                --get_negative_cluster {True if region == "calcarine" else False,}
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
def task_correlate_microregions() -> Dict:
    """
    Correlate the time series of each microregion with its image's Oz data.
    """
    regions = "calcarine occipital".split()
    for frequency in FREQUENCIES:
        for subject in SUBJECTS:
            for region in regions:
                for start_volume in START_VOLUMES:
                    average = fname.microregion_average(subject=subject, mask=region, start_volume=start_volume)
                    moving_moving_window_data = fname.out_tsv_name(subject=subject, frequency=frequency)
                    correlations = fname.microregions_correlation_results(subject=subject, mask=region, start_volume=start_volume, frequency=frequency)
                    table = fname.microregions_and_amplitudes(subject=subject, mask=region, start_volume=start_volume, frequency=frequency)
                    scatter = fname.microregions_correlation_scatter_plot(subject=subject, mask=region, start_volume=start_volume, frequency=frequency)

                    sources = [average, moving_moving_window_data]
                    targets = [correlations, scatter, table]

                    action = f"""
                        python3 correlate_regions.py
                        --load_amplitudes_from {moving_moving_window_data}
                        --load_microROI_from {average}
                        --save_scatter_to {scatter}
                        --save_table_to {table}
                        --save_spearman_to {correlations}
                    """.split()

                    yield dict(
                        name=f"subject--{subject}, mask--{region}, frequency--{frequency}, start_volume--{start_volume}",
                        actions=[action],
                        file_dep=sources,
                        targets=targets,
                    )
def task_correlate_across_subjects() -> Dict:
    """
    How big is the correlation across ALL subjects for each microROI?
    """
    regions = "calcarine occipital".split()

    for frequency in FREQUENCIES:
        for region in regions:
            for start_volume in START_VOLUMES:

                # Get sources.
                tables_to_be_catenated = [fname.microregions_and_amplitudes(subject=subject, mask=region, start_volume=start_volume, frequency=frequency) for subject in SUBJECTS]
                sources = tables_to_be_catenated

                # Get targets.
                correlations = fname.correlation_across_subjects(mask=region, frequency=frequency, start_volume=start_volume)
                table = fname.correlation_across_subjects_table(mask=region, frequency=frequency, start_volume=start_volume)
                scatter = fname.correlation_across_subjects_scatter(mask=region, frequency=frequency, start_volume=start_volume)
                targets = [correlations, table, scatter]

                # Get action.
                action = f"""
                    python3 correlate_all_microregions.py
                    --load_tables_from {" ".join(tables_to_be_catenated)}
                    --save_scatter_to {scatter}
                    --save_table_to {table}
                    --save_spearman_to {correlations}
                """.split()

                yield dict(
                    name=f"mask--{region}, frequency--{frequency}, start_volume--{start_volume}",
                    actions=[action],
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
