function k = interp_adf21(data, eb, ne, T)
%INTERP_ADF21 Interpolate effective stopping rate coefficient k(E,ne,T) [cm^3/s]
%
% Uses common ADAS factorization:
%   k(E,ne,T) = (k_T(T) / k_ref) * k_EN(E,ne)
%
% Inputs:
%   data:     : outputs from read_adf21, structure containing substructures:
%               svt, sven, beam, file 
%   eb        : beam energy [eV]
%   ne        : electron density [cm^-3]
%   T         : electron temperature [eV]
%
% Output:
%   k         : effective stopping rate coefficient [cm^3/s]

k_T     = interp1(data.svt.temp, data.svt.rate_coeff, T, 'linear', 'extrap');
k_ref   = data.sven.svref;
k_EN    = interp2(data.sven.eb, data.sven.dens, data.sven.rate_coeff.', eb, ne, 'linear');

k = (k_T / k_ref) * k_EN;
end
