function [svt, sven] = read_adf21(file_name)
%READ_ADF21 Read an ADAS ADF21 effective beam-stopping rate coefficient file.
%
% The beam stopping rate coefficient is represented by "k" in the
% following:
% dNb/dx = -ne*k*Nb
%
%   [SVT, SVEN] = READ_ADF21(FILE_NAME) parses an ADAS ADF21 file and returns
%   two scan structures:
%
%   SVEN : 2D scan over beam energy Eb and electron density ne at a reference
%          temperature (tref).
%          Fields:
%            .z           target plasma ion charge state (integer)
%            .spec        3-char plasma target species label
%            .svref       reference S*v coefficient [cm^3/s]
%            .tref        reference temperature [eV]
%            .eb          beam energy grid [eV]          (1 x neb)
%            .dens        electron density grid [cm^-3]  (1 x ndens)
%            .rate_coeff  S*v(Eb, ne) [cm^3/s]           (neb x ndens)
%
%   SVT  : 1D scan over temperature at reference Eb and ne.
%          Fields:
%            .ebref       reference beam energy [eV]
%            .denref      reference density [cm^-3]
%            .temp        temperature grid [eV]          (1 x ntemp)
%            .rate_coeff  S*v(T) [cm^3/s]                (1 x ntemp)
%
% Notes:
% - "SV” is commonly called the stopping value in ADAS beam-stopping tables.
% - Many ADF files use Fortran 'D' exponents; this reader converts D->E.
%
% Input
%   file_name (char/string): path to ADF21 file
%
% Output
%   svt, sven (struct)
%

file_name = char(file_name);

fid = fopen(file_name, 'r');
if fid < 0
    error('read_adf21:FileOpenFailed', 'Could not open file: %s', file_name);
end

% Automatic cleanup of fid file:
c = onCleanup(@() fclose(fid));

% --- helpers ---
    function line = nextLine()
        line = fgetl(fid);
        if ~ischar(line)
            error('read_adf21:UnexpectedEOF', 'Unexpected end-of-file in %s', file_name);
        end
        % Convert Fortran D exponent to E exponent for sscanf.
        line = strrep(line, 'D+', 'E+');
        line = strrep(line, 'd+', 'E+');
        line = strrep(line, 'D-', 'E-');
        line = strrep(line, 'd-', 'E-');
    end

    function vals = readNValues(n)
        % Read floats from as many lines as needed until we have n values.
        vals = zeros(1, n);
        got = 0;
        while got < n
            line = nextLine();
            v = sscanf(line, '%f').';
            if ~isempty(v)
                take = min(numel(v), n - got);
                vals(got+1:got+take) = v(1:take);
                got = got + take;
            end
        end
    end

% ===========================
% 2D energy and density scan
% ===========================

line = nextLine();

% Keep your fixed-column parse, but guard length.
if numel(line) < 22
    error('read_adf21:ParseError', 'Header line too short in %s', file_name);
end

sven.z = sscanf(line, '%d', 1);

% SCREF / SVREF (either is acceptable)
if contains(line,'SCREF=')
    sven.svref = read_keyword(line,'SCREF','%f');
elseif contains(line,'SVREF=')
    sven.svref = read_keyword(line,'SVREF','%f');
else
    error('ADF21 header does not contain SCREF or SVREF');
end

% SPEC is optional
sven.spec = read_keyword(line,'SPEC','%s','');

% Discard separator/comment line
nextLine();

% Read Eb and dens vectors:
line = nextLine();
neb   = str2double(strtrim(line(1:5)));
ndens = str2double(strtrim(line(6:10)));
sven.tref = read_keyword(line,'TREF','%f');
nextLine(); % Discard
sven.eb   = readNValues(neb);
sven.dens = readNValues(ndens);
nextLine(); % Discard

% Read rate coefficients: ndens columns, each column listed as neb values
sven.rate_coeff = zeros(neb, ndens);
for k = 1:ndens
    sven.rate_coeff(:, k) = readNValues(neb).';
end

nextLine(); % discard

% ===========================
% 1D temperature scan
% ===========================

line = nextLine();
ntemp     = str2double(strtrim(line(1:5)));
svt.ebref  = read_keyword(line,'EREF','%f');
svt.denref = read_keyword(line,'NREF','%f');
nextLine(); % discard
svt.temp = readNValues(ntemp);
nextLine(); % discard
svt.rate_coeff = readNValues(ntemp);

end
