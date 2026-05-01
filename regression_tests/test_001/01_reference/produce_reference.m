% PRODUCE_REFERENCE  Generate ADF21-based stopping mean-free-path reference table.
%
%   This script reads an ADAS ADF21 beam stopping dataset, constructs a 3D grid
%   of plasma conditions (electron density, beam energy per amu, and temperature),
%   computes the corresponding stopping mean free path λ_s for each condition,
%   saves the resulting reference table to an HDF5 file (reference.h5), and
%   produces diagnostic plots of λ_s vs temperature.
%
%   The ADF21 dataset is treated as the authoritative reference. Mean free paths
%   are derived from the ADF21 effective stopping coefficient by:
%
%       ks  = effective stopping coefficient interpolated from ADF21
%       vb  = sqrt(2 * e * Eb / M)   (converted to cm/s)
%       λ_s = vb / (ne * ks)
%
%   where:
%       Eb  is the beam energy per amu [eV/amu]
%       ne  is the electron density [cm^-3]
%       M   is the beam particle mass [kg]
%       e   is the elementary charge [C]
%
% INPUTS
%   input.nml
%     Fortran namelist defining the grid extents and resolution:
%       dene      : array of densities used only to set [nemin, nemax] [cm^-3]
%       nden      : number of density points
%       tmin,tmax : temperature bounds [eV]
%       nt        : number of temperature points
%       emin,emax : beam-energy-per-amu bounds [eV/amu]
%       neb       : number of energy points
%       beam_mass : beam mass [amu] (stored to output metadata)
%
%   ADF21 file (under input_dir, e.g. ./adas/)
%     Example: bms93#h_h1.dat
%     Read via get_data_adf21(), interpolated via interp_adf21().
%
% OUTPUTS
%   reference.h5  (when save_data = 1)
%     Numeric datasets:
%       /ne   : ne_vec     electron density vector [cm^-3]
%       /eb   : eb_vec     beam energy per amu vector [eV/amu]
%       /ab   : beam_mass  beam mass [amu]
%       /T    : t_vec      temperature vector [eV]  (Te = Ti assumption)
%       /mfp  : Ls         mean free path [cm], size (nden, neb, nt)
%                           with indexing Ls(i,j,k) ↔ (ne_vec(i), eb_vec(j), t_vec(k))
%     Root attributes (metadata):
%       description, creation_date, adf21_file
%     Dataset attributes include units/description and /mfp dimension order.
%
%   Figures (when save_figure = 1)
%     One figure per density value, plotting λ_s vs T for multiple beam energies.
%     Saved under ./figures/ using save_figure_to_file().
%
% PHYSICAL ASSUMPTIONS / NOTES
%   - Te = Ti (temperature equality) to match the assumptions used to produce the
%     ADF21-derived stopping coefficient tables.
%   - Energy input is in eV/amu; ensure consistency when changing beam isotope.
%
% DEPENDENCIES
%   ./utils/ (added to path)
%     - get_data_adf21()
%     - read_fortran_namelist()
%     - interp_adf21()
%     - write_reference_h5()
%     - save_figure_to_file()
%
% USAGE
%   1) Define ADfF21 file via file_name_adf21
%   2) Populate input.nml.
%   3) Run:
%        produce_reference

clear all
close all
clc

disp("Running script: ")
disp(mfilename + ".m")

% flags:
save_figure = 1;
save_data = 1;

% Define directories:
output_dir = "./";
input_dir = "./adas/";

% Include toolbox:
addpath(genpath("./utils/"))

%% DEFINE physical constants:
m_p = 1.6726e-27; % [kg]
e_c = 1.6020e-19; % [C]

%% DEFINE ADF21 file to use:
file_name_adf21 = "bms93#h_h1.dat"; % 1998
file_id = 'bms93_h_h1';

%% EXTRACT data:
% data.svt: Stopping coefficient over temperature
% data.sven: Stopping coefficient array over energy and density (neb, ndens) 
% data.beam: information on beam species and mass

file_path = fullfile(input_dir,file_name_adf21);
data = get_data_adf21(file_path);

%% READ input namelist:
input = read_fortran_namelist("input.nml");

nemin = min(input.dene); % [cm^-3]
nemax = max(input.dene);
n_dene = input.nden;

tmin = input.tmin; % [eV]
tmax = input.tmax;
n_temp = input.nt;
 
emin = input.emin; % [eV/amu]
emax = input.emax;
n_ener = input.neb;

%% VALIDATE inputs:

% Validate against ADF21 data limits:
nemin0 = min(data.sven.dens);
message = 'nemin needs to be >= %.2e [cm^-3] ';
assert(nemin >= nemin0, message, nemin0)

nemax0 = max(data.sven.dens);
message = 'nemax needs to be <= %.2e [cm^-3] ';
assert(nemax <= nemax0, message, nemax0)

tmin0 = min(data.svt.temp);
message = 'tmin needs to be >= %.2e [eV]';
assert(tmin >= tmin0, message, tmin0)

tmax0 = max(data.svt.temp);
message = 'tmax needs to be <= %.2e [eV]';
assert(tmax <= tmax0, message, tmax0)

emin0 = min(data.sven.eb);
message = 'emin needs to be >= %.2e [eV/amu]';
assert(emin >= emin0, message, emin0)

emax0 = max(data.sven.eb);
message = 'emax needs to be <= %.2e [eV/amu]';
assert(emax <= emax0, message, emax0)

% TODO: needs to validate against atomic_tables.h5 data too

%% CREATE grids:

ne_vec = input.dene;
t_vec  = linspace(tmin,tmax,n_temp);
eb_vec = linspace(emin,emax,n_ener);

%% COMPUTE: mean free path table (3D)

Ls = zeros(n_dene,n_ener,n_temp);
M = m_p*data.beam.mass_amu;
for dd = 1:n_dene
    for ee = 1:n_ener
        for tt = 1:n_temp

            % Conditions to interpolate at:
            EB = eb_vec(ee);
            NE = ne_vec(dd);
            TEI = t_vec(tt);

            % Get interpolated stopping rate coeff:
            ks = interp_adf21(data,EB,NE,TEI);

            % Neutral velocity:
            vb = (1e2)*sqrt(2*e_c*EB/M); % [cm/s]

            % Mean free path:
            Ls(dd,ee,tt) = vb/(NE*ks); % [cm]
        end
    end
end

%% SAVE data to HDF5:

if save_data
    disp(" ")
    disp("Saving data to HDF5 file 'reference.h5' ...")
    write_reference_h5("reference.h5", file_name_adf21, ...
        ne_vec, eb_vec, input.beam_mass, t_vec, Ls)
    disp("Data saving complete!")
end

%% PLOT: mean free path vrs temperature
font_size.ax = 14;
font_size.title = 16;
font_size.label = 16;
font_size.legend = 14;

line_color = {"k","r","g","bl","m","c"};

% Loop over density vales:
for dd = 1:numel(ne_vec)
    figure('color','w')
    hold on
    box on
    set(gca,'FontSize',font_size.ax)
    ymax = 0;
    ee_list = 1:1:n_ener;

    % Loop over energy:
    for ss = 1:numel(ee_list)
        ee = ee_list(ss);
        yval = (squeeze(Ls(dd,ee,:)))';
        ymax = max(max(yval),ymax);

        % Plot vs temperature:
        hls(ss) = plot(t_vec*1e-3,yval,line_color{ee},'linewidth',3);
        leg_str{ss} = "$n_e$: " + sprintf('%.1e',ne_vec(dd)) + " [cm$^{-3}$], " + ...
                      "$E_b$: " + sprintf('%.1e',eb_vec(ee)*1e-3) + " [keV]";
    end
    grid on
    
    hleg = legend(hls,leg_str);
    set(hleg,'Interpreter','latex','FontSize',font_size.legend)
    set(hleg,'Location','eastoutside')
    set(hleg,'Location','southoutside')
    
    title("Mean free path $\lambda_s$ [cm]",'Interpreter','latex',FontSize=font_size.title)
    ylabel("$\lambda_s$ [cm]", 'Interpreter','latex','FontSize',font_size.label)
    ylim([0,1.2]*ymax)
    xlabel("T [keV]", 'Interpreter','latex','FontSize',font_size.label)

    if save_figure
        disp(" ")
        disp("Saving figure ...")
        figure_name = file_id + "_MFP_" + dd;
        save_figure_to_file("./figures/",figure_name)
        disp("Saving complete!")
    end

end