%% random field thresholding for contrascan

% number of voxels in search: here we're using just the visual cortex via Kastner's atlas
% I calculated that using AFNI. Recorded the process in my personal notes.
numvox = 3815;

% search volume in mm^3
voxel_size = 2.5;
search_vol = numvox * (voxel_size^3);

% FWHM smoothing kernel, from preprocessing
kernel = 4;

% degrees of freedom (TRs correlated over x Subj avg to get data map)
alpha_num_trs = 366;
ssvep_num_trs = 365;
num_subjects = 18;
df_alpha = alpha_num_trs * num_subjects; 
df_ssvep = ssvep_num_trs * num_subjects;

% run stat_threshold. Last 3 parameters Maeve said to leave at default she specified.
[alpha_peak_threshold_tscore, alpha_extent_threshold] = stat_threshold(search_vol,numvox,kernel,df_alpha,0.05,0.001,0.05);
[ssvep_peak_threshold_tscore, ssvep_extent_threshold] = stat_threshold(search_vol,numvox,kernel,df_ssvep,0.05,0.001,0.05);

% convert t values into correlation coefficients
alpha_peak_threshold_correlation = sqrt(alpha_peak_threshold^2 / (alpha_peak_threshold^2 + df_alpha));
ssvep_peak_threshold_correlation = sqrt(ssvep_peak_threshold^2 / (ssvep_peak_threshold^2 + df_ssvep));

% alpha_peak_threshold_correlation = r of 0.051728349056573 -> 0.0517 -> .052
% ssvep_peak_threshold_correlation = r of 0.051798971471031 -> 0.0518 -> .052
% alpha_extent_threshold = 64.414123679024570 vox -> 64 vox
% ssvep_extent_threshold = 64.414123679024570 vox -> 64 vox
