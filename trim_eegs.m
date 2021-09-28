% Truncate an EEG dataset to start exactly 8s before the fMRI turned on.
%
% Created by Benjamin Velie on 7/27/2021.
% veliebm@ufl.edu
% --------------------------------------------------------------------------

parameters_file = 'processed/eeg_trimmed_to_trimmed_func_start/parameters.json';
do_all(parameters_file)
delete_lock_file(mfilename('fullpath'))

function do_all(parameters_file)
    % Trim all subjects whose metadata we've stored in a JSON file.
    all_parameters = read_json(parameters_file);

    for i = 1:numel(all_parameters)
        parameters = all_parameters(i)
        trim(parameters.subject, parameters.in_filename, parameters.in_dir, parameters.out_dir);
    end
end

function EEG = trim(subject, in_filename, in_dir, out_dir)
    eeglab;
    EEG = load_dataset(in_filename, in_dir);
    
    EEG = pop_rmdat(EEG, {'S  2'}, [0 999999], 0); % Trim so EEG starts at first stimulus.
    EEG = eeg_checkset( EEG );

    EEG = pop_rmdat(EEG, {'S  2'}, [-999999 10], 0); % Trim so EEG ends 10s after the final stimulus.
    EEG = eeg_checkset( EEG );
    
    EEG.setname=sprintf('sub-%s_truncated to fmri start time', subject);
    EEG = eeg_checkset( EEG );
    
    EEG = save_dataset(EEG, EEG.filename, out_dir);
end

function data = load_dataset(file_name, directory)
    % Load a dataset.
    data = pop_loadset('filename',file_name, 'filepath',directory);
    data = eeg_checkset( data );
end
function data = save_dataset(data, file_name, directory)
    % Save a dataset.
    data = pop_saveset(data, 'filename',file_name, 'filepath',directory);
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
