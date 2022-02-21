%% random field thresholding for contrascan whole brain in MNI space

% number of voxels in search: here we're using the whole entire brain.
% I calculated that using AFNI. Recorded the process in my notes.
numvox = 135137;

% search volume in mm^3
voxel_size = 2.5^3;
search_vol = numvox * voxel_size;

% FWHM smoothing kernel, from preprocessing
kernel = 4;

% degrees of freedom = number of subjects
num_subjects = 18;
df = num_subjects;

% run stat_threshold.
[peak_threshold_tscore, extent_threshold] = stat_threshold(search_vol,numvox,kernel,df);

% peak_threshold_correlation = 7.393187804421838 -> 7.3932
% extent_threshold = 125.3679030879341 -> 125.3679
