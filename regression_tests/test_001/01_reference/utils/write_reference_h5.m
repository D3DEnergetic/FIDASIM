function write_reference_h5(out_file, adf21_file, ne_vec, eb_vec, beam_mass, t_vec, Ls)
%WRITE_REFERENCE_H5: Write reference.h5 with metadata as attributes
%
% Numeric datasets:
%   /ne, /eb, /ab, /T, /mfp
%
% File-level attributes:
%   description, creation_date, adf21_file
%
% Expected shapes:
%   ne_vec : (nden,1)
%   en_vec : (neb,1)
%   t_vec  : (nt,1)
%   Ls     : (nden, neb, nt)  corresponds to (ne, eb, T)

    arguments
        out_file   (1,1) string
        adf21_file (1,1) string
        ne_vec     (:,1) double
        eb_vec     (:,1) double
        beam_mass  (1,1) double
        t_vec      (:,1) double
        Ls         (:,:,:) double
    end

    if isfile(out_file)
        delete(out_file); % avoid h5create collisions
    end

    % Density grid:
    h5create(out_file, "/ne", size(ne_vec));
    h5write (out_file, "/ne", ne_vec);
    h5writeatt(out_file, "/ne", "description", "electron density vector");
    h5writeatt(out_file, "/ne", "units", "[cm^-3]");

    % Beam energy grid:
    h5create(out_file, "/eb", size(eb_vec));
    h5write (out_file, "/eb", eb_vec);
    h5writeatt(out_file, "/eb", "description", "beam energy vector in eV/amu");
    h5writeatt(out_file, "/eb", "units", "[eV/amu]");

    % Beam mass:
    h5create(out_file, "/ab", 1);
    h5write (out_file, "/ab", beam_mass);
    h5writeatt(out_file, "/ab", "description", "beam mass in amu");
    h5writeatt(out_file, "/ab", "units", "amu");

    % Plasma temperature grid:
    h5create(out_file, "/T", size(t_vec));
    h5write (out_file, "/T", t_vec);
    h5writeatt(out_file, "/T", "description", "plasma temperature in eV (Te = Ti)");
    h5writeatt(out_file, "/T", "units", "[eV]");

    % Mean free path array:
    h5create(out_file, "/mfp", size(Ls));
    h5write (out_file, "/mfp", Ls);
    h5writeatt(out_file, "/mfp", "description", "mean free path");
    h5writeatt(out_file, "/mfp", "units", "[cm]");
    h5writeatt(out_file, "/mfp", "dimension_order", "(ne, eb, T) = (nden, neb, nt)");

    % Root-level attributes:
    h5writeatt(out_file, "/", "description", ...
        "Contains mean free paths derived from an ADF21 dataset.");
    h5writeatt(out_file, "/", "creation_date", ...
        char(datetime("now","Format","yyyy-MM-dd HH:mm:ss")));
    h5writeatt(out_file, "/", "adf21_file", char(adf21_file));
end