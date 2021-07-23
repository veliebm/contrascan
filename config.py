#!/usr/bin/env python3
"""
===========
Config file
===========

Configuration parameters for the study.
"""
import os
import getpass
from socket import getfqdn
from fnames import FileNames

###############################################################################
# Determine which user is running the scripts on which machine and set the path
# where the data is stored and how many CPU cores to use.

user = getpass.getuser()  # Username of the user running the scripts
host = getfqdn()  # Hostname of the machine running the scripts

# You want to add your machine to this list
if user == "csea":
    # CSEA desktop
    raw_data_dir = "./data"
    n_jobs = 4  # The CSEA desktop has 12 logical processors, but I don't want to use all of them at once. Gotta save some for other jobs.
else:
    # Defaults
    raw_data_dir = "./data"
    n_jobs = 1

# For BLAS to use the right amount of cores
os.environ["OMP_NUM_THREADS"] = str(n_jobs)


###############################################################################
# These are all the relevant parameters for the analysis.

# All subjects for whom our analysis can actually work.
SUBJECTS = "104 106 107 108 109 110 111 112 113 115 116 117 120 121 122 123 124 125".split()


###############################################################################
# Templates for filenames
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

# task_align_afniproc_irfs:
fname.add("atlas_template", "{raw_data_dir}/misc/kastner_cortex_masks/MNI152_T1_1mm.nii.gz")
fname.add("afniproc_template", "{raw_data_dir}/misc/MNI152_T1_2009c+tlrc.HEAD")
fname.add("afniproc_alignment_dir", "{processed_data_dir}/afniproc_alignment")
fname.add("afniproc_aligned_irf", "{afniproc_subject_dir}/{subject}.results/iresp_stim.{subject}_aligned+tlrc.HEAD")

# task_resample_afniproc_irfs:
fname.add("afniproc_resampled_irf", "{processed_data_dir}/afniproc_resample/sub-{subject}_IRF_resampled+tlrc.HEAD")

# task_smooth_fmriprep:
fname.add("fmriprep_func", "{fmriprep_subject_dir}/sub-{subject}/func/sub-{subject}_task-gabor_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz")
fname.add("fmriprep_smoothed", "{processed_data_dir}/fmriprep_smoothed/sub-{subject}_func_smoothed+tlrc.HEAD")

# task_scale_fmriprep:
fname.add("fmriprep_scaled", "{processed_data_dir}/fmriprep_scaled/sub-{subject}_func_scaled+tlrc.HEAD")

# task_deconvolve_fmriprep:
fname.add("fmriprep_deconvolve_dir", "{processed_data_dir}/fmriprep_deconvolved/sub-{subject}")
fname.add("fmriprep_deconvolved", "{fmriprep_deconvolve_dir}/sub-{subject}_deconvolved+tlrc.HEAD")
fname.add("fmriprep_irf", "{fmriprep_deconvolve_dir}/sub-{subject}_IRF+tlrc.HEAD")
fname.add("fmriprep_regressors_dir", "{fmriprep_deconvolve_dir}/regressors")
fname.add("regressors_tsv", "{fmriprep_subject_dir}/sub-{subject}/func/sub-{subject}_task-gabor_desc-confounds_timeseries.tsv")

# task_align_fmriprep_irfs:
fname.add("fmriprep_template", "{raw_data_dir}/misc/tpl-MNI152NLin2009cAsym_res-01_desc-brain_T1w.nii.gz")
fname.add("fmriprep_alignment_dir", "{processed_data_dir}/fmriprep_alignment")
fname.add("fmriprep_aligned_irf", "{fmriprep_deconvolve_dir}/sub-{subject}_IRF_aligned+tlrc.HEAD")

# task_resample_fmriprep_irfs:
fname.add("fmriprep_resampled_irf", "{processed_data_dir}/fmriprep_resample/sub-{subject}_IRF_resampled+tlrc.HEAD")

# task_ttest_fmriprep_vs_afniproc:
fname.add("ttest_result", "{processed_data_dir}/fmriprep_vs_afniproc_ttests/subbrick-{subbrick}_seta-fmriprep_setb-afniproc_ttest+tlrc.HEAD")

# task_trim_func_images:
fname.add("trimmed_dir", "{processed_data_dir}/trimmedfuncs")
fname.add("trimmed_func", "{trimmed_dir}/sub-{subject}_func_trimmed+tlrc.HEAD")
fname.add("afniproc_onsets", "{afniproc_subject_dir}/onsets.tsv")
fname.add("afniproc_func", "{afniproc_subject_dir}/{subject}.results/all_runs.{subject}+tlrc.HEAD")

# task_convert_eeg
fname.add("converteeg_dir", "{processed_data_dir}/converteeg")
fname.add("converted_eeg", "{converteeg_dir}/sub-{subject}_eeg.set")
fname.add("brainvision_eeg", "{brainvision_dir}/contrascan_{subject}_Pulse Artifact Correction.vhdr")
