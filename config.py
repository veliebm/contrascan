#!/usr/bin/env python3
"""
===========
Config file
===========

Configuration parameters for the study.
"""
import os
from fnames import FileNames

###############################
# Set number of cores to use and stuff.

raw_data_dir = "./data"
n_jobs = 1

# For BLAS to use the right amount of cores
os.environ["OMP_NUM_THREADS"] = str(n_jobs)


###############################################################################
# These are all the relevant parameters for the analysis.

# All subjects for whom our analysis can actually work.
SUBJECTS = "104 106 107 108 109 110 111 112 113 115 116 117 120 121 122 123 124 125".split()

# Which frequencies to use in our frequencies analyses.
FREQUENCIES = "12 24".split()

# Which volumes to start our EEG/fMRI correlation from.
START_VOLUMES = range(1, 5)

# Which ICA components to remove for each subject.
COMPONENTS_TO_REMOVE = {
    "102": [13],
    "104": [3, 4, 5, 6],
    "106": [1, 2, 4, 5, 6],
    "107": [2, 4, 5],
    "108": [2, 4, 10],
    "109": [2, 6, 15],
    "110": [2, 3, 9],
    "111": [2],
    "112": [4],
    "113": [2, 6, 11, 20, 24],
    "115": [3, 5, 6, 13, 14],
    "116": [5, 6, 7, 12],
    "117": [2, 3, 4, 9],
    "120": [1, 4, 6, 10],
    "121": [1, 3, 6, 7, 12, 13],
    "122": [5, 6, 7, 14, 20, 23],
    "123": [3, 16, 27],
    "124": [3, 10, 22],
    "125": [2, 5, 6, 10],
    "126": [1, 4, 9, 28]
}

###############################################################################
# Templates for filenames part of our main analysis.
#
# This part of the config file uses the FileNames class. It provides a small
# wrapper around string.format() to keep track of a list of filenames.
# See fnames.py for details on how this class works.
fname = FileNames()

# Raw data directories.
fname.add("raw_data_dir", raw_data_dir)
fname.add("raw_subject_dir", "{raw_data_dir}/subjects-complete/sub-{subject}")
fname.add("fmriprep_subject_dir", "{raw_data_dir}/fmriprep/sub-{subject}/fmriprep")
fname.add("brainvision_dir", "{raw_data_dir}/brainvision/EXPORT")

# Processed data directory.
fname.add("processed_data_dir", "./processed")

# Script directory.
fname.add("scripts_dir", ".")

# task_check:
fname.add("system_check", "./system_check.txt")

# task_create_bids_root:
fname.add("bids_dir", "{processed_data_dir}/bids")
fname.add("bids_description", "{bids_dir}/dataset_description.json")

# task_bidsify_subject:
fname.add("raw_eeg", "{raw_subject_dir}/contrascan_{subject}.eeg")
fname.add("raw_vmrk", "{raw_subject_dir}/contrascan_{subject}.vmrk")
fname.add("raw_vhdr", "{raw_subject_dir}/contrascan_{subject}.vhdr")
fname.add("raw_dat", "{raw_subject_dir}/SFcontrascan_{subject}.dat")
fname.add("raw_func", "{raw_subject_dir}/Keil_{subject}_EPI_2s_Gain.nii")
fname.add("raw_anat", "{raw_subject_dir}/Keil_{subject}_sT1W_3D_FFE_SAG.nii")
fname.add("bids_dat", "{bids_dir}/sourcedata/sub-{subject}_task-contrascan.dat")
fname.add("bids_subject_dir", "{bids_dir}/sub-{subject}")
fname.add("bids_anat", "{bids_subject_dir}/anat/sub-{subject}_T1w.nii")
fname.add("bids_eeg", "{bids_subject_dir}/eeg/sub-{subject}_task-contrascan_eeg.eeg")
fname.add("bids_vmrk", "{bids_subject_dir}/eeg/sub-{subject}_task-contrascan_eeg.vmrk")
fname.add("bids_vhdr", "{bids_subject_dir}/eeg/sub-{subject}_task-contrascan_eeg.vhdr")
fname.add("bids_func_json", "{bids_subject_dir}/func/sub-{subject}_task-contrascan_bold.json")
fname.add("bids_func", "{bids_subject_dir}/func/sub-{subject}_task-contrascan_bold.nii")
fname.add("bids_events", "{bids_subject_dir}/func/sub-{subject}_task-contrascan_events.tsv")

# task_afniproc:
fname.add("afniproc_dir", "{processed_data_dir}/afniproc")
fname.add("afniproc_subject_dir", "{afniproc_dir}/sub-{subject}")
fname.add("afniproc_log", "{afniproc_subject_dir}/output.proc.{subject}")
fname.add("afniproc_command", "{afniproc_subject_dir}/proc.{subject}")
fname.add("afniproc_deconvolved", "{afniproc_subject_dir}/{subject}.results/stats.{subject}+tlrc.HEAD")
fname.add("afniproc_irf", "{afniproc_subject_dir}/{subject}.results/iresp_stim.{subject}+tlrc.HEAD")
fname.add("afniproc_anat", "{afniproc_subject_dir}/{subject}.results/anat_final.{subject}+tlrc.HEAD")
fname.add("afniproc_preprocessed_func", "{afniproc_subject_dir}/{subject}.results/all_runs.{subject}+tlrc.HEAD")

# task_resample_deconvolutions:
fname.add("afniproc_deconvolved_resampled", "{afniproc_subject_dir}/{subject}.results/stats.{subject}_resampled+tlrc.HEAD")

# task_ttest_deconvolutions:
fname.add("afniproc_ttest_result", "{processed_data_dir}/afniproc_ttests/subbrick-{subbrick}+tlrc.HEAD")

# task_align_func_images:
fname.add("afniproc_template", "{raw_data_dir}/misc/MNI152_T1_2009c+tlrc.HEAD")
fname.add("alignment_dir", "{processed_data_dir}/func_alignment")
fname.add("aligned_func", "{afniproc_subject_dir}/{subject}.results/all_runs.{subject}_aligned+tlrc.HEAD")

# task_resample_template
fname.add("resampled_template", "{processed_data_dir}/kastner_resample/MNI152_T1_1mm_resampled+tlrc.HEAD")

# task_align_afniproc_irfs:
fname.add("atlas_template", "{raw_data_dir}/misc/kastner_cortex_masks/MNI152_T1_1mm.nii.gz")
fname.add("afniproc_alignment_dir", "{processed_data_dir}/IRF_alignment")
fname.add("afniproc_aligned_irf", "{afniproc_subject_dir}/{subject}.results/iresp_stim.{subject}_aligned+tlrc.HEAD")

# task_resample_afniproc_irfs:
fname.add("afniproc_resampled_irf", "{processed_data_dir}/IRF_resample/sub-{subject}_IRF_resampled+tlrc.HEAD")

# task_resample_func_images:
fname.add("resampled_func", "{processed_data_dir}/func_resample/sub-{subject}_func_resampled+tlrc.HEAD")

# task_trim_func_images:
fname.add("trimmed_dir", "{processed_data_dir}/funcs_trimmed_to_first_stim")
fname.add("trimmed_func", "{trimmed_dir}/sub-{subject}_func_trimmed+tlrc.HEAD")
fname.add("afniproc_onsets", "{afniproc_subject_dir}/onsets.tsv")
fname.add("afniproc_func", "{afniproc_subject_dir}/{subject}.results/all_runs.{subject}+tlrc.HEAD")

# task_trim_func_images_again:
fname.add("final_funcs_dir", "{processed_data_dir}/func_final_trim")
fname.add("final_func", "{final_funcs_dir}/sub-{subject}_startvolume-{start_volume}_func+tlrc.HEAD")

# task_prepare_to_convert_eeg:
fname.add("converteeg_dir", "{processed_data_dir}/eeg_convert_from_brainvision")
fname.add("converteeg_json", "{converteeg_dir}/parameters.json")
fname.add("brainvision_eeg", "{brainvision_dir}/contrascan_{subject}_Pulse Artifact Correction.vhdr")

# task_convert_eeg:
fname.add("converted_eeg", "{converteeg_dir}/sub-{subject}_eeg.set")

# task_prepare_to_preprocess_eeg:
fname.add("preprocesseeg_dir", "{processed_data_dir}/eeg_preprocess")
fname.add("preprocesseeg_json", "{preprocesseeg_dir}/parameters.json")

# task_preprocess_eeg:
fname.add("preprocessed_eeg", "{preprocesseeg_dir}/sub-{subject}_eeg.set")

# task_prepare_to_segment_eeg:
fname.add("segmenteeg_dir", "{processed_data_dir}/eeg_segment")
fname.add("segmenteeg_json", "{segmenteeg_dir}/parameters.json")

# task_segment_eeg:
fname.add("segmented_eeg", "{segmenteeg_dir}/sub-{subject}_eeg.set")

# task_prepare_to_trim_eeg:
fname.add("trimeeg_dir", "{processed_data_dir}/eeg_trimmed_to_trimmed_func_start")
fname.add("trimeeg_json", "{trimeeg_dir}/parameters.json")

# task_segment_eeg:
fname.add("trimmed_eeg", "{trimeeg_dir}/sub-{subject}_eeg.set")

# task_freqtag_calculate_parameters:
fname.add("freqtag_parameters_dir", "{processed_data_dir}/freqtag_parameters")
fname.add("freqtag_parameters_script", "{scripts_dir}/TEMP_get_params.m")
fname.add("freqtag_faxis", "{freqtag_parameters_dir}/faxis.mat")
fname.add("freqtag_faxisall", "{freqtag_parameters_dir}/faxisall.mat")
fname.add("freqtag_stimulus_start", "{freqtag_parameters_dir}/stimulus_start.mat")
fname.add("freqtag_stimulus_end", "{freqtag_parameters_dir}/stimulus_end.mat")
fname.add("freqtag_epoch_duration", "{freqtag_parameters_dir}/epoch_duration.mat")
fname.add("freqtag_sampling_rate", "{freqtag_parameters_dir}/sampling_rate.mat")

# task_freqtag_fft:
fname.add("freqtag_fft_script", "{scripts_dir}/TEMP_sub{subject}_freqtag_fft.m")
fname.add("freqtag_fft_dir", "{processed_data_dir}/freqtag_fft")
fname.add("freqtag_fft_phase", "{freqtag_fft_dir}/sub-{subject}_phase.mat")
fname.add("freqtag_fft_pow", "{freqtag_fft_dir}/sub-{subject}_pow.mat")
fname.add("freqtag_fft_freqs", "{freqtag_fft_dir}/sub-{subject}_freqs.mat")

# task_freqtag_3d_fft:
fname.add("freqtag_3d_fft_script", "{scripts_dir}/TEMP_sub{subject}_freqtag_3d_fft.m")
fname.add("freqtag_3d_fft_dir", "{processed_data_dir}/freqtag_3d_fft")
fname.add("freqtag_3d_fft_spec", "{freqtag_3d_fft_dir}/sub-{subject}_spec.mat")

# task_freqtag_sliding_window:
fname.add("freqtag_sliding_window_script", "{scripts_dir}/TEMP_sub{subject}_{frequency}Hz_sliding_window.m")
fname.add("freqtag_sliding_window_dir", "{processed_data_dir}/freqtag_sliding_window_dir")
fname.add("freqtag_sliding_window_trialpow", "{freqtag_sliding_window_dir}/sub-{subject}_frequency-{frequency}_trialpow.mat")
fname.add("freqtag_sliding_window_winmat3d", "{freqtag_sliding_window_dir}/sub-{subject}_frequency-{frequency}_winmat3d.mat")
fname.add("freqtag_sliding_window_phasestabmat", "{freqtag_sliding_window_dir}/sub-{subject}_frequency-{frequency}_phasestabmat.mat")
fname.add("freqtag_sliding_window_trialSNR", "{freqtag_sliding_window_dir}/sub-{subject}_frequency-{frequency}_trialSNR.mat")
fname.add("freqtag_sliding_window_outfile", "{freqtag_sliding_window_dir}/sub-{subject}_frequency-{frequency}_outfile.slidwin.mat")
fname.add("freqtag_sliding_window_meanwinmat", "{freqtag_sliding_window_dir}/sub-{subject}_frequency-{frequency}_meanwinmat.mat")

# task_freqtag_fft_sliding_window:
fname.add("freqtag_fft_sliding_window_script", "{scripts_dir}/TEMP_sub{subject}_{frequency}Hz_FFT_sliding_window.m")
fname.add("freqtag_fft_sliding_window_dir", "{processed_data_dir}/freqtag_FFT_sliding_window_dir")
fname.add("freqtag_fft_sliding_window_pow", "{freqtag_fft_sliding_window_dir}/sub-{subject}_frequency-{frequency}_pow.mat")
fname.add("freqtag_fft_sliding_window_phase", "{freqtag_fft_sliding_window_dir}/sub-{subject}_frequency-{frequency}_phase.mat")
fname.add("freqtag_fft_sliding_window_freqs", "{freqtag_fft_sliding_window_dir}/sub-{subject}_frequency-{frequency}_freqs.mat")

# task_freqtag_better_sliding_window:
fname.add("freqtag_better_sliding_window_script", "{scripts_dir}/TEMP_sub{subject}_better_sliding_window.m")
fname.add("freqtag_better_sliding_window_dir", "{processed_data_dir}/freqtag_better_sliding_window")
fname.add("freqtag_better_sliding_window_trialpow", "{freqtag_better_sliding_window_dir}/sub-{subject}_trialpow.mat")
fname.add("freqtag_better_sliding_window_winmat3d", "{freqtag_better_sliding_window_dir}/sub-{subject}_winmat3d.mat")
fname.add("freqtag_better_sliding_window_phasestabmat", "{freqtag_better_sliding_window_dir}/sub-{subject}_phasestabmat.mat")
fname.add("freqtag_better_sliding_window_trialSNR", "{freqtag_better_sliding_window_dir}/sub-{subject}_trialSNR.mat")
fname.add("freqtag_better_sliding_window_outfile", "{freqtag_better_sliding_window_dir}/sub-{subject}_outfile.slidwin.mat")
fname.add("freqtag_better_sliding_window_meanwinmat", "{freqtag_better_sliding_window_dir}/sub-{subject}_meanwinmat.mat")

# task_freqtag_better_fft_sliding_window:
fname.add("freqtag_better_fft_sliding_window_script", "{scripts_dir}/TEMP_sub{subject}_better_FFT_sliding_window.m")
fname.add("freqtag_better_fft_sliding_window_dir", "{processed_data_dir}/freqtag_better_FFT_sliding_window")
fname.add("freqtag_better_fft_sliding_window_pow", "{freqtag_better_fft_sliding_window_dir}/sub-{subject}_pow.mat")
fname.add("freqtag_better_fft_sliding_window_phase", "{freqtag_better_fft_sliding_window_dir}/sub-{subject}_phase.mat")
fname.add("freqtag_better_fft_sliding_window_freqs", "{freqtag_better_fft_sliding_window_dir}/sub-{subject}_freqs.mat")


# task_freqtag_hilbert:
fname.add("freqtag_hilbert_script", "{scripts_dir}/TEMP_sub{subject}_{frequency}Hz_hilbert.m")
fname.add("freqtag_hilbert_dir", "{processed_data_dir}/freqtag_hilbert")
fname.add("freqtag_hilbert_powermat", "{freqtag_hilbert_dir}/sub-{subject}_frequency-{frequency}_powermat.mat")
fname.add("freqtag_hilbert_phasemat", "{freqtag_hilbert_dir}/sub-{subject}_frequency-{frequency}_phasemat.mat")
fname.add("freqtag_hilbert_complexmat", "{freqtag_hilbert_dir}/sub-{subject}_frequency-{frequency}_complexmat.mat")

# task_mean_mean_fft:
fname.add("mean_mean_fft", "{processed_data_dir}/mean_mean_fft/frequency-{frequency}_mean_mean_FFT.tsv")

# task_mean_mean_hilbert:
fname.add("mean_mean_hilbert", "{processed_data_dir}/mean_mean_hilbert/frequency-{frequency}_mean_mean_hilbert.tsv")

# task_prepare_to_moving_moving_window_eeg:
fname.add("movingmovingwindoweeg_dir", "{processed_data_dir}/eeg_moving_moving_window")
fname.add("movingmovingwindoweeg_json", "{movingmovingwindoweeg_dir}/parameters.json")

# task_moving_moving_window_eeg:
fname.add("moving_moving_windowed_eeg", "{movingmovingwindoweeg_dir}/sub-{subject}_frequency-{frequency}_moving_moving_window_average.amp.at")
fname.add("amplitudes", "{movingmovingwindoweeg_dir}/sub-{subject}_frequency-{frequency}_moving_moving_window_average.tsv")

# task_correlate_eeg_fmri:
fname.add("correlation_dir", "{processed_data_dir}/correlation_whole_brain")
fname.add("correlation_image", "{correlation_dir}/sub-{subject}_frequency-{frequency}_startvolume-{start_volume}_correlation.nii")

# task_ttest_eeg_fmri_correlations:
fname.add("correlations_ttest_dir", "{processed_data_dir}/correlation_whole_brain_ttest")
fname.add("correlations_ttest", "{correlations_ttest_dir}/frequency-{frequency}_startvolume-{start_volume}_correlations_ttest+tlrc.HEAD")

# task_get_occipital_mask:
fname.add("kastner_mask_dir", "{raw_data_dir}/misc/kastner_cortex_masks/subj_vol_all")
fname.add("kastner_mask", "{kastner_mask_dir}/perc_VTPM_vol_roi{roi_number}_{hemisphere}h.nii.gz")
fname.add("masks_dir", "{processed_data_dir}/masks")
fname.add("occipital_pole_mask", "{masks_dir}/occipital_pole+tlrc.HEAD")

# task_get_calcarine_mask:
fname.add("calcarine_mask", "{masks_dir}/calcarine+tlrc.HEAD")

# task_resample_masks:
fname.add("resampled_mask", "{processed_data_dir}/masks_resampled/mask-{mask}_resampled+tlrc.HEAD")

# task_apply_masks_to_correlations:
fname.add("correlations_ttest_masked", "{correlations_ttest_dir}/startvolume-{start_volume}_frequency-{frequency}_masked-{mask}_correlations_ttest+tlrc.HEAD")

# task_apply_masks_to_irfs:
fname.add("masked_irf", "{processed_data_dir}/masked_irfs/sub-{subject}_masked-{mask}_irf+tlrc.HEAD")

# task_clusterize_irfs:
fname.add("clusters_dir", "{processed_data_dir}/clusters")
fname.add("clusters", "{clusters_dir}/sub-{subject}_source-{mask}_clusters+tlrc.HEAD")
fname.add("clusters_summary", "{clusters_dir}/sub-{subject}_source-{mask}_summary.1D")

# task_make_micromasks:
fname.add("micromask", "{processed_data_dir}/micromasks/sub-{subject}_source-{mask}_micromask+tlrc.HEAD")

# task_apply_micromasks_to_trimmed_trimmed_funcs:
fname.add("micromasked_func", "{processed_data_dir}/micromasked_funcs/sub-{subject}_source-{mask}_startvolume-{start_volume}_bold_micromasked+tlrc.HEAD")

# task_average_microregion_voxels:
fname.add("microregion_average", "{processed_data_dir}/microregion_averages/sub-{subject}_source-{mask}_startvolume-{start_volume}_average.txt")

# task_correlate_microregions:
fname.add("microregions_correlation_dir", "{processed_data_dir}/microregion_correlations")
fname.add("microregions_correlation_results", "{microregions_correlation_dir}/sub-{subject}_frequency-{frequency}_source-{mask}_startvolume-{start_volume}_microregion_correlations.txt")
fname.add("microregions_and_amplitudes", "{microregions_correlation_dir}/sub-{subject}_frequency-{frequency}_source-{mask}_startvolume-{start_volume}_microregion+amplitudes.csv")
fname.add("microregions_correlation_scatter_plot", "{microregions_correlation_dir}/sub-{subject}_frequency-{frequency}_source-{mask}_startvolume-{start_volume}_scatter.png")

# task_correlate_across_subjects:
fname.add("correlation_across_subjects_dir", "{processed_data_dir}/correlation_across_subjects")
fname.add("correlation_across_subjects", "{correlation_across_subjects_dir}/frequency-{frequency}_source-{mask}_startvolume-{start_volume}_correlation.txt")
fname.add("correlation_across_subjects_scatter", "{correlation_across_subjects_dir}/frequency-{frequency}_source-{mask}_startvolume-{start_volume}_scatter.png")
fname.add("correlation_across_subjects_table", "{correlation_across_subjects_dir}/frequency-{frequency}_source-{mask}_startvolume-{start_volume}_microregions+amplitudes.csv")
