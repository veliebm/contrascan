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

# Processed data directories.
fname.add("processed_data_dir", "./processed")

# Task: Create BIDS root.
fname.add("bids_root", "{processed_data_dir}/bids")
fname.add("bids_description", "{bids_root}/dataset_description.json")

# The data files that are used and produced by the analysis steps
fname.add("output", "{processed_data_dir}/output-{subject}.txt")
fname.add("grand_average", "{processed_data_dir}/grand_average.txt")

# File produced by check_system.py
fname.add("system_check", "./system_check.txt")
