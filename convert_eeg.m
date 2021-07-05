% Converts our brainvision files into regular EEGLAB files.
%
% Created 6/28/21 by Benjamin Velie.
% ---------------------------------------------------------

function convert_eeg(in_dir, in_filename, out_path, setname)
    % Convert a subject from BrainVision format to EEGLAB format.
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    EEG = eeg_checkset( EEG );
    EEG = pop_loadbv(in_dir, in_filename);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname',setname,'savenew',out_path,'gui','off');
end
