% Preprocesses our EEG files using ICA.
%
% Created by Ben Velie on 4/6/2021.
% veliebm@ufl.edu
% --------------------------------------------------------------------------

parameters_file = 'processed/preprocesseeg/parameters.json';
preprocess_all(parameters_file)

function preprocess_all(parameters_file)
    % Preprocesses all subjects whose metadata we've stored in a JSON file.
    all_parameters = read_json(parameters_file);

    for i = 1:numel(all_parameters)
        parameters = all_parameters(i)
        filter_and_ICA_and_rereference_and_save_dataset(parameters.subject, parameters.in_file, parameters.in_directory, parameters.out_directory, parameters.components_to_remove);
    end
end


function filter_and_ICA_and_rereference_and_save_dataset(subject, in_file, in_directory, out_directory, components_to_remove)
    eeglab;
    EEG = load_dataset(in_file, in_directory);
    
    %plot_dataset(EEG);
    %savefig(sprintf('sub-%s_before ICA.fig', subject))
    
    EEG = fir_filter(EEG);
    EEG = auto_ICA(EEG, components_to_remove);
    
    %plot_dataset(EEG);
    %savefig(sprintf('sub-%s_after ICA.fig', subject))
    
    EEG = rereference(EEG);
    
    %plot_dataset(EEG);
    %savefig(sprintf('sub-%s_after rereferencing.fig', subject))
    
    EEG = save_dataset(EEG, EEG.filename, out_directory);
end

function values = extract_json(path)
    % Extract values from a json file. Returns a structure.
    fid = fopen(path); 
    raw = fread(fid,inf); 
    string = char(raw'); 
    fclose(fid); 
    values = jsondecode(string);
end
function data = rereference(data)
    % Rereference the dataset.
    data = pop_reref( data, []);
    data = eeg_checkset( data );
end
function data = auto_ICA(data, components_to_remove)
    % Run a SOBI ICA analysis on the dataset.
    data = pop_runica(data, 'icatype', 'sobi');
    data = eeg_checkset(data);
    data = pop_subcomp(data, components_to_remove, 0);
    data = eeg_checkset(data);
end
function data = fir_filter(data)
    % FIR filter a dataset. The low frequency cutoff is 1 Hz, and the high frequency cutoff is 36 Hz.
    data = pop_eegfiltnew(data, 'locutoff',1,'hicutoff',36,'plotfreqz',0);
    data = eeg_checkset( data );
end
function plot_dataset(data)
    % Plot a dataset.
    pop_eegplot(data, 1, 1, 1);
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
