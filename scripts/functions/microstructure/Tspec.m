function SK=Tspec(spc,a,b,c,d,q)
%TSPEC Computes theoretical spectrum.
%
%   INPUTS:
%   spc (string): type of spectrum, 'K' (Kraichnan) or 'B' (Batchelor).
%   a (double or 1*2 double array): Xi parameter if double, [Xi,kB]
%   parameters if double array; with Xi=dissipation of temperature variance
%   and kB=Batchelor wavenumber.
%   b (double or 1*n double array): kB parameter is a is double, K coordinate (component of the wavenumber vector along one direction) if a is a double
%   array.
%   c (1*n double array, not needed if a is a double array): K coordinate
%   d (double, not needed if a is a double array): molecular thermal diffusivity k 
%
%
%   OUTPUTS:
%   SK (1*n double array): theoretical spectral densities. 
%
% T. Doda based on S. Piccolroaz, 06.02.2025
%%

if nargin>=4
    Xi=a;
    KB=b;
    K=c;
    k=d;
elseif nargin==3 % In that case, the thermal diffusivity k is missing! (T. Doda, 06.02.25)
    error('Not enough arguments!') % added by T. Doda
    % Xi=a(1);
    % KB=a(2);
    % K=b;
end

% See Sanchez et al. (2011) DOI: 10.1175/JPO-D-11-025.1
if spc=='K'  % Kraichnan
    phi=KB/sqrt(2*q);
    y=K/phi;
    f=y.*exp(-sqrt(3)*y); %Kraichnan
    SK=Xi/(2*k*KB)*sqrt(2*q)*f;
elseif spc=='B' % Batchelor
    phi=KB/sqrt(2*q);
    y=K/phi;
    f=y.*(exp(-y.^2/2)-y.*sqrt(pi()/2).*(1-erf(y/sqrt(2)))); %Batchelor
    SK=Xi/(2*k*KB)*sqrt(2*q)*f;
end
end
