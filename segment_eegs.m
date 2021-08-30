% Segment a preprocessed EEG subject.
%
% Created by Ben Velie on 7/26/21.
% veliebm@ufl.edu
% --------------------------------------------------------------------------

parameters_file = 'processed/segmenteeg/parameters.json';
do_all(parameters_file)
delete_lock_file(mfilename('fullpath'))

function do_all(parameters_file)
    % Preprocesses all subjects whose metadata we've stored in a JSON file.
    all_parameters = read_json(parameters_file);

    for i = 1:numel(all_parameters)
        parameters = all_parameters(i)
        do_one(parameters.in_name, parameters.in_directory, parameters.out_directory);
    end
end
function do_one(file_name, in_directory, out_directory)
    eeglab;
    EEG = load_dataset(file_name, in_directory);
    
    EEG = extract_epochs(EEG);
        
    EEG = save_dataset(EEG, EEG.filename, out_directory);
end

function data = extract_epochs(data)
    % Split our dataset into epochs anchored at "S  2" stimuli with the following range: (-800ms, +4500ms)
    data = pop_epoch( data, {  'S  2'  }, [-0.8           4.5], 'epochinfo', 'yes');
    data = eeg_checkset( data );
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
