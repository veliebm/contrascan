% Creates a MOVING moving window average for each subject.
%
% Driver code created by Benjamin Velie on 7/27/2021.
% veliebm@ufl.edu
%
% stead2singtrialsCont and regressionMAT created by Andreas Keil in May 2021.
% akeil@ufl.edu
% --------------------------------------------------------------------------

parameters_file = 'processed/eeg_moving_moving_window/parameters.json';
do_all(parameters_file)
delete_lock_file(mfilename('fullpath'))

%% Structural functions.
function do_all(parameters_file)
    % Moving moving window all subjects whose metadata we've stored in a JSON file.
    all_parameters = read_json(parameters_file);

    for i = 1:numel(all_parameters)
        parameters = all_parameters(i)
        do_one(parameters.in_filename, parameters.in_dir, parameters.out_stem, parameters.out_tsv_name, str2num(parameters.frequency));
    end
end
function do_one(in_filename, in_dir, out_stem, out_tsv_name, frequency)
    % Run stead2singtrialsCont on a subject.
    eeglab;
    stead2singtrialsCont(in_filename, in_dir, 0, 1:1000, 1:1000, frequency, 600, 500, 1, out_stem);
    convert_to_tsv(strcat(out_stem, '.amp.at'), out_tsv_name)
end

%% Input/output functions.
function [data] = read_json(in_path)
    % Read a JSON file.
    fname = in_path; 
    fid = fopen(fname); 
    raw = fread(fid,inf); 
    str = char(raw'); 
    fclose(fid); 
    data = jsondecode(str);
end
function convert_to_tsv(in_filename, out_filename)
    % Convert an emegs file to tsv.
    amplitude = ReadAvgFile(in_filename);
    dlmwrite(out_filename, amplitude, '\t');
end
