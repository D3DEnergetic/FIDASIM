function data = get_data_adf21(file_path)
%GET_DATA_ADF21 Read a single ADF21 file and return a unified data structure.
%
%   DATA = GET_DATA_ADF21(FILE_PATH) reads one ADAS ADF21 beam-stopping
%   dataset and returns all relevant information in a single structure.
%
%   Input
%   -----
%   file_path : char or string
%       Full or relative path to an ADF21 data file.
%
%   Output
%   ------
%   data : struct with fields
%       .file
%           .path   full file path
%           .name   file name only
%
%       .beam
%           Output of PARSE_ADF21_BEAM_SPECIES (beam symbol, mass, etc.)
%
%       .svt
%           Temperature scan structure from READ_ADF21
%
%       .sven
%           Energy–density scan structure from READ_ADF21
%
%   Notes
%   -----
%   - This function operates on a *single* ADF21 file.
%   - Looping over multiple files is intentionally left to the caller.
%   - Requires: read_adf21.m, parse_adf21_beam_species.m

    % Validate input:
    if ~(ischar(file_path) || (isstring(file_path) && isscalar(file_path)))
        error('file_path must be a char or string scalar');
    end

    file_path = char(file_path);

    % File info:
    [folder, name, ext] = fileparts(file_path);
    data.file.path = file_path;
    data.file.name = [name ext];
    data.file.dir  = folder;

    % Beam info:
    data.beam = parse_adf21_beam_species(file_path);

    % Read ADF21 data:
    [data.svt, data.sven] = read_adf21(file_path);
end