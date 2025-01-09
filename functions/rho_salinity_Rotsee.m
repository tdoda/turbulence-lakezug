function [rhoTS,Sal,depth]=rho_salinity_Rotsee(Temp,Cond,Press)

%
%******************************INFO ***********************************
%
% [rhoTS,Sal,depth]=rho_salinity_Garda(Temp,Cond)
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
%***********************************************************************

%Conversion from conductivity to salinity (specific to Rotsee):
%
%cond25 = fT_function25(temp).*cond
%sal25 = C25_to_S *cond25
%
%with: 
%C25_to_S = 0.79
%fT_function25(t)=1.8692 -0.055181.*Temp + 0.0011096.*Temp.^2-1.1762E-005.*Temp.^3

% constants
BetaS=0.807e-3; % Haline contraction coefficient [�^-1]

% Convert conductivity (C) to C_20�C
fT=1.8692 -0.055181.*Temp + 0.0011096.*Temp.^2-1.1762E-005.*Temp.^3;
Cond25=fT.*Cond; % C_20�C / Cond MUST be in uS/cm

% calculate salinity
Sal=0.79*Cond25/1000; % [g/kg]

% caluclate rho0
%rho0=1000-7e-3*(Temp-4).^2; %[kg m^-3]
rho0 = 999.84298 + 1e-3*(65.4891*Temp - 8.56272*Temp.^2 + 0.059385*Temp.^3); %doesn't work for temperature shiger than 25�C

% calculate rho as a function of Temperature and Salinity
rhoTS=rho0.*(1+BetaS*Sal); %[kg m^-3] 

% calculate depth
mrho = cumsum(rhoTS)./[1:length(rhoTS)]';
depth=10000*Press./(mrho*9.81);
