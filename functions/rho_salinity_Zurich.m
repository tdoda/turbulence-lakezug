function [rhoTS,Sal,depth]=rho_salinity_Zurich(Temp,Cond,Press)
%
%******************************INFO ***********************************
%
% [rhoTS,Sal,depth]=rho_salinity_Geneva(Temp,Cond)
% 
% I N P U T S
% Temp: Temperature [deg C]
% Cond: Conductivity [uS/cm] (standard for freshwaters)
% Press: Pressure [dbar] (already corrected for atmosferic pressure)
% 
% O T P U T S
% rhoTS: density [kg m^-3] as a function of Temperature and Salinity
% Salinity: Salinity [�]
% depth: depth [m]
%
% Modified by Bieito Fernandez
% Oscar Sep�lveda Steiner (oscar.sepulvedasteiner@epfl.ch) APHYS - EPFL
%
%***********************************************************************

% constants
BetaS=0.807e-3; % Haline contraction coefficient [�^-1]

% Convert conductivity (C) to C_20�C
fT=1.8626 - 0.0052908*Temp + 0.00093057.*Temp.^2;
Cond20=fT.*Cond; % C_20�C / Cond MUST be in uS/cm

% calculate salinity
Sal=0.7999e-3*Cond20; % [�]

% caluclate rho0
%rho0=1000-7e-3*(Temp-4).^2; %[kg m^-3]
rho0 = 999.84298 + 1e-3*(65.4891*Temp - 8.56272*Temp.^2 + 0.059385*Temp.^3);% doesn't work for temperature shiger than 25�C

% calculate rho as a function of Temperature and Salinity
rhoTS=rho0.*(1+BetaS*Sal); %[kg m^-3] 

% calculate depth
mrho = cumsum(rhoTS)./[1:length(rhoTS)]';
depth=10000*Press./(mrho*9.81);
