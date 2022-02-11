#!/usr/bin/env python3
"""
===========
Config file
===========

Configuration parameters for the study.
"""
from fnames import FileNames

###############################
# Set raw data dir.

raw_data_dir = "./data"


###############################################################################
# These are all the relevant parameters for the analysis.

# All subjects for whom our analysis can actually work.
SUBJECTS = "104 106 107 108 109 110 111 112 113 115 116 117 120 121 122 123 124 125".split()

# Which frequencies to use in our frequencies analyses.
FREQUENCIES = "12 24".split()

# Number of permutations to use for thresholding.
PERMUTATIONS = range(1, 201)

# Which volumes to start our EEG/fMRI correlation from.
START_VOLUMES = range(1, 5)
EXPANDED_START_VOLUMES = range(1, 9)

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

# task_ttest_IRFs:
fname.add("ttested_IRF", "{processed_data_dir}/IRF_ttests/subbrick-{subbrick}+tlrc.HEAD")

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

# task_get_IRF_mean:
fname.add("IRF_mean", "{processed_data_dir}/IRF_mean/IRF_mean+tlrc.HEAD")
fname.add("deconvolve_mean", "{processed_data_dir}/IRF_mean/deconvolve_mean+tlrc.HEAD")

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

# task_get_trial_func_amplitudes:
fname.add("func_trial_amplitudes_dir", "{processed_data_dir}/func_trial_amplitudes")
fname.add("func_trial_amplitudes", "{func_trial_amplitudes_dir}/sub-{subject}_trials.nii")

# task_convert_eeg:
fname.add("converteeg_dir", "{processed_data_dir}/eeg_convert_from_brainvision")
fname.add("converteeg_script", "{scripts_dir}/TEMP_convert_eeg_sub{subject}.m")
fname.add("converteeg_json", "{converteeg_dir}/parameters.json")
fname.add("brainvision_eeg", "{brainvision_dir}/contrascan_{subject}_Pulse Artifact Correction.vhdr")
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

# task_trim_eeg:
fname.add("trimmed_eeg", "{trimeeg_dir}/sub-{subject}_eeg.set")

# task_eeg_get_flicker_frequencies:
fname.add("eeg_flicker_frequencies_dir", "{processed_data_dir}/eeg_flicker_frequencies")
fname.add("eeg_flicker_frequencies_script", "{scripts_dir}/TEMP_sub{subject}_eeg_get_flicker_frequencies.m")
fname.add("eeg_flicker_frequencies", "{eeg_flicker_frequencies_dir}/sub-{subject}_variable-{variable}.mat")

# task_eeg_get_alphas:
fname.add("eeg_alpha_dir", "{processed_data_dir}/eeg_alphas")
fname.add("eeg_alpha_script", "{scripts_dir}/TEMP_sub{subject}_eeg_get_alphas.m")
fname.add("eeg_alpha", "{eeg_alpha_dir}/sub-{subject}_data-{data}_alphas.mat")
fname.add("average_power", "{eeg_alpha_dir}/sub-{subject}_average_power.mat")

# task_eeg_get_trial_by_trial_alpha:
fname.add("eeg_trial_alpha_dir", "{processed_data_dir}/eeg_trial_alphas")
fname.add("eeg_trial_alpha_script", "{scripts_dir}/TEMP_sub{subject}_eeg_get_trial_alphas.m")
fname.add("eeg_trial_alpha", "{eeg_trial_alpha_dir}/sub-{subject}_data-{data}_trial_alpha.mat")

# task_correlate_alpha_and_snr:
fname.add("eeg_correlation_alpha_snr_dir", "{processed_data_dir}/correlation_alpha_and_SNRs")
fname.add("eeg_correlation_alpha_snr_table", "{eeg_correlation_alpha_snr_dir}/sub-{subject}_table.csv")
fname.add("eeg_correlation_alpha_snr_results", "{eeg_correlation_alpha_snr_dir}/sub-{subject}_results.txt")
fname.add("eeg_correlation_alpha_snr_scatter", "{eeg_correlation_alpha_snr_dir}/sub-{subject}_scatter.png")

# task_correlate_alpha_and_snr_across_subjects
fname.add("eeg_correlation_across_subjects_alpha_snr_table", "{eeg_correlation_alpha_snr_dir}/across_subjects_table.csv")
fname.add("eeg_correlation_across_subjects_alpha_snr_results", "{eeg_correlation_alpha_snr_dir}/across_subjects_results.txt")
fname.add("eeg_correlation_across_subjects_alpha_snr_scatter", "{eeg_correlation_alpha_snr_dir}/across_subjects_scatter.png")

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
fname.add("freqtag_better_sliding_window_channel", "{freqtag_better_sliding_window_dir}/sub-{subject}_variable-{variable}_channel-{channel}.mat")

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

# task_eeg_sliding_sliding_window:
fname.add("eeg_sliding_sliding_window_script", "{scripts_dir}/TEMP_sub{subject}_{frequency}Hz_eeg_sliding_sliding_window.m")
fname.add("eeg_sliding_sliding_window_dir", "{processed_data_dir}/eeg_sliding_sliding_window")
fname.add("eeg_sliding_sliding_window_amplitudes", "{eeg_sliding_sliding_window_dir}/sub-{subject}_frequency-{frequency}_amplitudes.fftamp.mat")
fname.add("eeg_sliding_sliding_window_SNR", "{eeg_sliding_sliding_window_dir}/sub-{subject}_frequency-{frequency}_SNR.mat")
fname.add("eeg_sliding_sliding_window_oz_SNR", "{eeg_sliding_sliding_window_dir}/sub-{subject}_frequency-{frequency}_oz_SNR.mat")
fname.add("eeg_sliding_sliding_window_oz_amplitudes", "{eeg_sliding_sliding_window_dir}/sub-{subject}_frequency-{frequency}_oz_amplitudes.mat")

# task_improve_sliding_sliding_window:
fname.add("eeg_sliding_sliding_window_improved_dir", "{processed_data_dir}/eeg_sliding_sliding_window_improved")
fname.add("eeg_sliding_sliding_window_improved", "{eeg_sliding_sliding_window_improved_dir}/sub-{subject}_variable-{variable}_improved.mat")

# task_correlate_whole_brain:
fname.add("correlation_dir", "{processed_data_dir}/correlation_whole_brain")
fname.add("correlation_image", "{correlation_dir}/sub-{subject}_frequency-{frequency}_startvolume-{start_volume}_correlation.nii")
fname.add("correlation_whole_brain_SNR_image", "{correlation_dir}/sub-{subject}_frequency-{frequency}_startvolume-{start_volume}_SNR_correlation.nii")
fname.add("correlation_whole_brain_alpha", "{correlation_dir}/sub-{subject}_startvolume-{start_volume}_data-{data}_alphas.nii")
fname.add("correlation_whole_brain_improved", "{correlation_dir}/sub-{subject}_startvolume-{start_volume}_variable-{variable}_improved.nii")
fname.add("correlation_whole_brain_trials", "{correlation_dir}/sub-{subject}_variable-{variable}_trials.nii")
fname.add("correlation_whole_brain_canonical", "{correlation_dir}/sub-{subject}_canonical.nii")
fname.add("correlation_whole_brain_permutation", "{correlation_dir}/permutation_testing/sub-{subject}/{analysis}/sub-{subject}_startvolume-{start_volume}_variable-{variable}_{analysis}_{permutation}.nii.gz")

# task_fisher_transform_whole_brain:
fname.add("correlation_fisher_dir", "{processed_data_dir}/correlation_whole_brain_fisher")
fname.add("correlation_fisher", "{correlation_fisher_dir}/subject-{subject}_startvolume-{start_volume}_variable-{variable}_{analysis}_fisher.nii")

# task_ttest_whole_brain_fishers:
fname.add("correlations_ttest_fisher_dir", "{processed_data_dir}/correlation_whole_brain_fisher_ttest")
fname.add("correlations_ttest_fisher", "{correlations_ttest_fisher_dir}/startvolume-{start_volume}_variable-{variable}_{analysis}_fisher_ttest.nii")

# task_cohens_d_whole_brain:
fname.add("correlations_cohens_fisher_dir", "{processed_data_dir}/correlation_whole_brain_fisher_cohens")
fname.add("correlations_cohens_fisher", "{correlations_cohens_fisher_dir}/startvolume-{start_volume}_variable-{variable}_{analysis}_fisher_cohens.nii")

# task_ttest_whole_brain_correlations:
fname.add("correlations_ttest_dir", "{processed_data_dir}/correlation_whole_brain_ttest")
fname.add("correlations_ttest", "{correlations_ttest_dir}/frequency-{frequency}_startvolume-{start_volume}_correlations_ttest+tlrc.HEAD")
fname.add("correlations_SNR_ttest", "{correlations_ttest_dir}/frequency-{frequency}_startvolume-{start_volume}_correlations_ttest_SNR+tlrc.HEAD")
fname.add("correlations_whole_brain_alpha_ttest", "{correlations_ttest_dir}/startvolume-{start_volume}_data-{data}_correlations_alpha_ttest+tlrc.HEAD")
fname.add("correlations_improved_whole_brain_ttest", "{correlations_ttest_dir}/startvolume-{start_volume}_variable-{variable}_improved_correlations_ttest+tlrc.HEAD")
fname.add("correlations_whole_brain_trials_ttest", "{correlations_ttest_dir}/variable-{variable}_trials_correlations_ttest+tlrc.HEAD")
fname.add("correlations_whole_brain_canonical_ttest", "{correlations_ttest_dir}/canonical_correlations_ttest+tlrc.HEAD")
fname.add("correlations_whole_brain_permutations_ttest", "{correlations_ttest_dir}/permutation_testing/{analysis}/startvolume-{start_volume}_variable-{variable}_{analysis}_{permutation}+tlrc.HEAD")

# task_get_occipital_mask:
fname.add("kastner_mask_dir", "{raw_data_dir}/misc/kastner_cortex_masks/subj_vol_all")
fname.add("kastner_mask", "{kastner_mask_dir}/perc_VTPM_vol_roi{roi_number}_{hemisphere}h.nii.gz")
fname.add("masks_dir", "{processed_data_dir}/masks")
fname.add("occipital_pole_mask", "{masks_dir}/occipital_pole+tlrc.HEAD")

# task_get_calcarine_mask:
fname.add("calcarine_mask", "{masks_dir}/calcarine+tlrc.HEAD")

# task_resample_masks:
fname.add("resampled_mask", "{processed_data_dir}/masks_resampled/mask-{mask}_resampled+tlrc.HEAD")

# task_apply_masks_to_irfs:
fname.add("masked_irf", "{processed_data_dir}/masked_irfs/sub-{subject}_masked-{mask}_irf+tlrc.HEAD")

# task_get_average_IRF_in_ROTs:
fname.add("average_irf_in_roi_dir", "{processed_data_dir}/IRF_averages_in_rois")
fname.add("average_irf_in_roi", "{average_irf_in_roi_dir}/sub-{subject}_masked-{mask}.txt")

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

# task_correlate_eeg_with_average_microregion_timeseries:
fname.add("microregions_correlation_dir", "{processed_data_dir}/microregion_correlations")
fname.add("microregions_correlation_results", "{microregions_correlation_dir}/sub-{subject}_frequency-{frequency}_source-{mask}_startvolume-{start_volume}_microregion_correlations.txt")
fname.add("microregions_and_amplitudes", "{microregions_correlation_dir}/sub-{subject}_frequency-{frequency}_source-{mask}_startvolume-{start_volume}_microregion+amplitudes.csv")
fname.add("microregions_correlation_scatter_plot", "{microregions_correlation_dir}/sub-{subject}_frequency-{frequency}_source-{mask}_startvolume-{start_volume}_scatter.png")
fname.add("microregions_correlation_SNR_results", "{microregions_correlation_dir}/sub-{subject}_frequency-{frequency}_source-{mask}_startvolume-{start_volume}_microregion_SNR_correlations.txt")
fname.add("microregions_and_SNR_amplitudes", "{microregions_correlation_dir}/sub-{subject}_frequency-{frequency}_source-{mask}_startvolume-{start_volume}_microregion+SNRamplitudes.csv")
fname.add("microregions_correlation_SNR_scatter_plot", "{microregions_correlation_dir}/sub-{subject}_frequency-{frequency}_source-{mask}_startvolume-{start_volume}_SNR_scatter.png")
fname.add("microregions_correlation_alpha_results", "{microregions_correlation_dir}/sub-{subject}_source-{mask}_startvolume-{start_volume}_data-{data}_microregion_alpha_correlations.txt")
fname.add("microregions_and_alpha_amplitudes", "{microregions_correlation_dir}/sub-{subject}_source-{mask}_startvolume-{start_volume}_data-{data}_microregion+alpha_table.csv")
fname.add("microregions_correlation_alpha_scatter_plot", "{microregions_correlation_dir}/sub-{subject}_source-{mask}_startvolume-{start_volume}_data-{data}_alpha_scatter.png")
fname.add("microregions_correlation_better", "{microregions_correlation_dir}/sub-{subject}_source-{mask}_startvolume-{start_volume}_variable-{variable}_better_correlations.txt")
fname.add("microregions_and_amplitudes_better", "{microregions_correlation_dir}/sub-{subject}_source-{mask}_startvolume-{start_volume}_variable-{variable}_better_table.csv")
fname.add("microregions_correlation_scatter_plot_better", "{microregions_correlation_dir}/sub-{subject}_source-{mask}_startvolume-{start_volume}_variable-{variable}_better_scatter.png")

# task_correlate_eeg_with_average_microregion_timeseries_across_subjects:
fname.add("correlation_across_subjects_dir", "{processed_data_dir}/correlation_across_subjects")
fname.add("correlation_across_subjects", "{correlation_across_subjects_dir}/frequency-{frequency}_source-{mask}_startvolume-{start_volume}_correlation.txt")
fname.add("correlation_across_subjects_scatter", "{correlation_across_subjects_dir}/frequency-{frequency}_source-{mask}_startvolume-{start_volume}_scatter.png")
fname.add("correlation_across_subjects_table", "{correlation_across_subjects_dir}/frequency-{frequency}_source-{mask}_startvolume-{start_volume}_microregions+amplitudes.csv")
fname.add("correlation_across_subjects_SNR", "{correlation_across_subjects_dir}/frequency-{frequency}_source-{mask}_startvolume-{start_volume}_correlation_SNR.txt")
fname.add("correlation_across_subjects_SNR_scatter", "{correlation_across_subjects_dir}/frequency-{frequency}_source-{mask}_startvolume-{start_volume}_scatter_SNR.png")
fname.add("correlation_across_subjects_SNR_table", "{correlation_across_subjects_dir}/frequency-{frequency}_source-{mask}_startvolume-{start_volume}_microregions+amplitudes_SNR.csv")
fname.add("correlation_across_subjects_alpha", "{correlation_across_subjects_dir}/source-{mask}_startvolume-{start_volume}_data-{data}_correlation_alpha.txt")
fname.add("correlation_across_subjects_alpha_scatter", "{correlation_across_subjects_dir}/source-{mask}_startvolume-{start_volume}_data-{data}_scatter_alpha.png")
fname.add("correlation_across_subjects_alpha_table", "{correlation_across_subjects_dir}/source-{mask}_startvolume-{start_volume}_data-{data}_microregions+alpha_table.csv")
fname.add("correlation_across_subjects_better", "{correlation_across_subjects_dir}/source-{mask}_startvolume-{start_volume}_variable-{variable}_correlation_better.txt")
fname.add("correlation_across_subjects_scatter_better", "{correlation_across_subjects_dir}/source-{mask}_startvolume-{start_volume}_variable-{variable}_scatter_better.png")
fname.add("correlation_across_subjects_table_better", "{correlation_across_subjects_dir}/source-{mask}_startvolume-{start_volume}_variable-{variable}_table_better.csv")

# task_mask_and_average_correlations:
fname.add("correlation_averages_dir", "{processed_data_dir}/correlation_averages")
fname.add("correlation_alpha_average", "{correlation_averages_dir}/sub-{subject}_startvolume-{start_volume}_source-{mask}_data-{data}_correlation_average_alpha.txt")
fname.add("correlation_improved_sliding_sliding_window_average", "{correlation_averages_dir}/sub-{subject}_startvolume-{start_volume}_source-{mask}_variable-{variable}_correlation_average_improved_slidslidwin.txt")
fname.add("correlation_permutation_average", "{correlation_averages_dir}/permutations/sub-{subject}/{analysis}/sub-{subject}_startvolume-{start_volume}_source-{mask}_variable-{variable}_analysis-{analysis}_correlation_average_{permutation}.txt")

# task_ttest_averages:
fname.add("ttest_averages_dir", "{processed_data_dir}/correlation_averages_ttests")
fname.add("ttest_averages", "{ttest_averages_dir}/variable-{variable}_startvolume-{start_volume}_mask-{mask}_analysis-{analysis}_outfile-{outfile}.csv")

# task_average_averages:
fname.add("average_averages_dir", "{processed_data_dir}/correlation_average_averages")
fname.add("average_averages_correlation", "{average_averages_dir}/variable-{variable}_startvolume-{start_volume}_mask-{mask}_analysis-{analysis}_{outfile}.csv")
fname.add("average_averages_correlation_permutation", "{average_averages_dir}/permutations/{analysis}/variable-{variable}_startvolume-{start_volume}_mask-{mask}_analysis-{analysis}_{outfile}_{permutation}.csv")

# task_get_canonical_bold:
fname.add("canonical_dir", "{processed_data_dir}/canonical_bold_response")
fname.add("canonical", "{canonical_dir}/sub-{subject}_canonical.mat")
fname.add("canonical_script", "{scripts_dir}/TEMP_sub{subject}_canonical.m")

# task_trim_canonical_bold:
fname.add("canonical_trimmed_dir", "{processed_data_dir}/canonical_bold_response_trimmed")
fname.add("canonical_trimmed", "{canonical_trimmed_dir}/sub-{subject}_canonical.mat")

# task_trim_canonical_bold_again:
fname.add("canonical_lagged_dir", "{processed_data_dir}/canonical_bold_response_lagged")
fname.add("canonical_lagged", "{canonical_lagged_dir}/sub-{subject}_startvolume-{start_volume}_canonical.mat")

# task_subtract_canonical_bold:
fname.add("compared_to_canonical_dir", "{processed_data_dir}/compare_to_canonical")
fname.add("compared_to_canonical", "{compared_to_canonical_dir}/startvolume-{start_volume}_variable-{variable}_{analysis}_baselinecorrected.nii")
fname.add("compared_permutations_to_canonical", "{compared_to_canonical_dir}/permutations/{analysis}/startvolume-{start_volume}_variable-{variable}_{analysis}_baselinecorrected_{permutation}.nii")

# task_calculate_variance:
fname.add("variance_whole_brain_dir", "{processed_data_dir}/variance_whole_brain")
fname.add("variance_whole_brain", "{variance_whole_brain_dir}/startvolume-{start_volume}_variable-{variable}_baselined-false_{analysis}.nii.gz")

# task_subtract_canonical_variance:
fname.add("variance_whole_brain_baselined", "{variance_whole_brain_dir}/startvolume-{start_volume}_variable-{variable}_baselined-true_{analysis}.nii.gz")

# task_scramble_data:
fname.add("scrambled_series_dir", "{processed_data_dir}/scrambled_series")
fname.add("scrambled_series", "{scrambled_series_dir}/sub-{subject}/{analysis}/sub-{subject}_startvolume-{start_volume}_variable-{variable}_{analysis}_scrambled_{permutation}.mat")

# task_get_maxes_and_mins:
fname.add("maxes_mins_dir", "{processed_data_dir}/maxes_and_mins_distributions")
fname.add("maxes_mins_table", "{maxes_mins_dir}/startvolume-{start_volume}_variable-{variable}_baselined-{baselined}_{analysis}_{outfile}_table.csv")

# plot_maxes_and_mins_distributions:
fname.add("maxes_mins_plot", "{maxes_mins_dir}/startvolume-{start_volume}_variable-{variable}_{analysis}_{outfile}_plot.png")

# task_threshold_results:
fname.add("thresholding_results_dir", "{processed_data_dir}/thresholding_results")
fname.add("threshold_outfile", "{thresholding_results_dir}/variable-{variable}_startvolume-{start_volume}_mask-{mask}_analysis-{analysis}_{outfile}.csv")

# task_cat_permutations:
fname.add("catenated_permutations_dir", "{processed_data_dir}/catenated_permutations")
fname.add("catenated_permutations", "{catenated_permutations_dir}/variable-{variable}_startvolume-{start_volume}_analysis-{analysis}_catenated.nii.gz")

# task_calc_percentiles:
fname.add("whole_brain_percentiles_dir", "{processed_data_dir}/whole_brain_percentiles")
fname.add("whole_brain_percentiles", "{whole_brain_percentiles_dir}/variable-{variable}_startvolume-{start_volume}_analysis-{analysis}_percentiles.nii.gz")
