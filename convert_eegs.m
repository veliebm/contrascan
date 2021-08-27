% Converts a list of brainvision files into regular EEGLAB files.
%
% Created 7/23/21 by Benjamin Velie.
% ---------------------------------------------------------

parameters_file = 'processed/converteeg/parameters.json';
do_all(parameters_file)
delete_lock_file(mfilename('fullpath'))


function do_all(parameters_file)
    % Convert many subjects from BrainVision format to EEGLAB format. Reads a JSON file to get parameters.

    all_parameters = read_json(parameters_file);

    for i = 1:numel(all_parameters)
        parameters = all_parameters(i)
        do_one(parameters.brainvision_dir, parameters.brainvision_name, parameters.converted_path, parameters.setname)
    end
end

function do_one(in_dir, in_filename, out_path, setname)
    % Convert a subject from BrainVision format to EEGLAB format.
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    EEG = eeg_checkset( EEG );
    EEG = pop_loadbv(in_dir, in_filename);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname',setname,'savenew',out_path,'gui','off');
end

function [data] = read_json(in_path)
    % Read a JSON file.
    fname = in_path; 
    fid = fopen(fname); 
    raw = fread(fid,inf); 
    str = char(raw'); 
    fclose(fid); 
    data = jsondecode(str);
end
