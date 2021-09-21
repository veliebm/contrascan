function [filename, directory] = split_path(path)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Returns the filename and directory components of a path.
    %
    % Inputs
    % ------
    % path : character vector
    %   Path to a file.
    %
    % Returns
    % -------
    % filename : character vector
    %   Name of the file.
    % directory : character vector
    %   Path to the directory containing the file.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [directory, stem, suffix] = fileparts(path);
    filename = strcat(stem, suffix);
end
