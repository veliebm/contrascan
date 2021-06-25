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
    n_jobs = 12  # The CSEA desktop has 12 logical processors
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

# Task: Create BIDS root.
fname.add("bids_dir", "{processed_data_dir}/bids")
fname.add("bids_description", "{bids_dir}/dataset_description.json")

# Task: Bidsify subjects.
fname.add("raw_eeg", "{raw_subject_dir}/contrascan_{subject}.{eeg_type}")
fname.add("raw_dat", "{raw_subject_dir}/SFcontrascan_{subject}.dat")
fname.add("raw_func", "{raw_subject_dir}/Keil_{subject}_EPI_2s_Gain.nii")
fname.add("raw_anat", "{raw_subject_dir}/Keil_{subject}_sT1W_3D_FFE_SAG.nii")
fname.add("bids_dat", "{bids_dir}/sourcedata/sub-{subject}_task-contrascan.dat")
fname.add("bids_subject_dir", "{bids_dir}/sub-{subject}")
fname.add("bids_anat", "{bids_subject_dir}/anat/sub-{subject}_T1w.nii")
fname.add("bids_eeg", "{bids_subject_dir}/eeg/sub-{subject}_task-contrascan_eeg.{eeg_type}")
fname.add("bids_func_json", "{bids_subject_dir}/func/sub-{subject}_task-contrascan_bold.json")
fname.add("bids_func", "{bids_subject_dir}/func/sub-{subject}_task-contrascan_bold.nii")
fname.add("bids_events", "{bids_subject_dir}/func/sub-{subject}_task-contrascan_events.tsv")

# The data files that are used and produced by the analysis steps
fname.add("output", "{processed_data_dir}/output-{subject}.txt")
fname.add("grand_average", "{processed_data_dir}/grand_average.txt")

# File produced by check_system.py
fname.add("system_check", "./system_check.txt")
