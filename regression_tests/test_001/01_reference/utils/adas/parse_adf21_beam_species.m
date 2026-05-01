function beam = parse_adf21_beam_species(filename)
%PARSE_ADF21_BEAM_SPECIES Extract beam species and mass from ADF21 filename
%
% Supported forms:
%   h_h1.dat
%   he_ne10.dat
%   bms97#h_h1.dat
%   bms10#he_ne10.dat
%
% Output:
%   beam.symbol   (char)
%   beam.mass_amu (double)

    % strip path and extension
    [~, name, ~] = fileparts(filename);

    % if dataset prefix exists, remove it
    parts = split(name, '#');
    if numel(parts) == 2
        core = parts{2};
    else
        core = parts{1};
    end

    % split beam and target
    species = split(core, '_');
    if numel(species) < 2
        error('Filename "%s" does not match <beam>_<target>.dat pattern', filename);
    end

    beam_str = lower(species{1});

    % map beam species to atomic mass (amu)
    switch beam_str
        case 'h'
            beam.symbol   = 'H';
            beam.mass_amu = 1.00784;
        case 'd'
            beam.symbol   = 'D';
            beam.mass_amu = 2.01410;
        case 't'
            beam.symbol   = 'T';
            beam.mass_amu = 3.01605;
        case 'he'
            beam.symbol   = 'He';
            beam.mass_amu = 4.00260;
        default
            error('Unknown beam species "%s" in filename "%s"', beam_str, filename);
    end
end