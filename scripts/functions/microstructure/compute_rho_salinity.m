function [rhoTS,Cmatch,Cond20,Sal,depth] = compute_rho_salinity(lakename,Temp,Cond,Press,JAC_correction)
%COMPUTE_RHO_SALINITY Computes salinity, density and depth based on
%temperature, conductivity and pressure for a given lake.
%
%   INPUTS:
%   lakename (str): name of the lake (options: "Geneva", "default")
%   Temp (1*n double array): Temperature [deg C]
%   Cond (1*n double array): Conductivity [uS/cm] (standard for freshwater)
%   Press (1*n double array): Pressure [dbar] (already corrected for atmospheric pressure)
%   JAC_Correction (boolean) [optional, default = false]: run the ODAS
%   function salinity_JAC to match the JAC_C signal to the JAC_T signal and
%   reduce errors.
%
%   OUTPUTS:
%   rhoTS (1*n double array): density [kg m^-3] as a function of Temperature and Salinity
%   Cmatch (1*n double array): corrected conductivity [uS/cm] to match
%   temperature data response
%   cond20 (1*n double array): conductivity at 20°C [uS/cm] 
%   Sal (1*n double array): Salinity [g kg^-1]
%   depth (1*n double array): depth [m] (not corrected for atmospheric
%   pressure)
%
% T.Doda, based on S. Piccolroaz and B. Fernandez Castro, 10.12.2024
%% Parameters (*** add values used to compute cond20***)
if nargin<5
    JAC_correction=false;
end

Cond(Cond<0)=0; % Replace negative conductivity values by zero (T. Doda)
Cond(~isfinite(Cond))=0; % Replace non finite conductivity values by zero (T. Doda)
Cond20=NaN(size(Cond));

if strcmp(lakename,'default')
    [Sal,Cmatch]=salinity_JAC(Press,Temp,Cond);
    Cmatch=Cmatch*1000; % [uS/cm]
    rhoTS=1000+sigma_p(Temp,Sal,Press); % [kg/m3]
    mrho = cumsum(rhoTS)./(1:length(rhoTS))';
    depth=10000*Press./(mrho*9.81);
    return

elseif strcmp(lakename,'Geneva')
    BetaS=0.807e-3; % Haline contraction coefficient [kg/g]
    gS=0.874e-3; % Coefficient converting conductivity at 20°C to salinity [g cm kg^-1 uS^-1]
    g=9.81; % gravitational acceleration, [m s^-2]
elseif strcmp(lakename,'Zug')
    % *********** NEED TO CHECK AND UPDATE THE ESTIMATES **************
    gS=0.79e-3; %  Coefficient converting conductivity at 20°C to salinity [g cm kg^-1 uS^-1]  (Boisgontier, 2016)
    
    rho0 = 999.84298 + 1e-3*(65.4891*Temp - 8.56272*Temp.^2 + 0.059385*Temp.^3); % [kg/m3], doesn't work for temperature higher than 25°C
    if JAC_correction
        [~,Cmatch] = salinity_JAC(Press,Temp,Cond);
        Cmatch=Cmatch*1000; % [uS/cm]
    else
        Cmatch=Cond;
    end
    % Convert conductivity (C) to C_20°C
    [Cond20,~,~]=conductivity_20(Temp,Cmatch); % [uS/cm]
    Sal=gS*Cond20; % [g/kg]
    rhoTS=rho0.*(1+gS*1e-3*Cond20); % [kg/m3]
    mrho = cumsum(rhoTS)./(1:length(rhoTS))';
    depth=10000*Press./(mrho*9.81);
    return
else
    error('Wrong lake name!')
end

%% Computation (*** could be updated based on Python scripts ***)
if JAC_correction
    [~,Cmatch] = salinity_JAC(Press,Temp,Cond);
    Cmatch=Cmatch*1000; % [uS/cm]
else
    Cmatch=Cond;
end

% Convert conductivity (C) to C_20°C
[Cond20,~,~]=conductivity_20(Temp,Cmatch); % [uS/cm]

% calculate salinity
Sal=gS*Cond20; % [g/kg]

% calculate rho0
%rho0=1000-7e-3*(Temp-4).^2; %[kg m^-3]
rho0 = 999.84298 + 1e-3*(65.4891*Temp - 8.56272*Temp.^2 + 0.059385*Temp.^3); %doesn't work for temperature higher than 25°C

% calculate rho as a function of Temperature and Salinity
rhoTS=rho0.*(1+BetaS*Sal); %[kg m^-3]

% calculate depth
mrho = cumsum(rhoTS)./(1:length(rhoTS))';
depth=10000*Press./(mrho*g);


end

function [Cond20,alpha,dCdT]=conductivity_20(Temp,Cond)
% CAUTION: could vary with ionic composition and should ideally be updated
% rom one lake to another (T. Doda)
    %
    %******************************INFO ***********************************
    %
    % [Cond20,alpha,dCdT]=conductivity_20(Temp,Cond)
    % 
    % I N P U T S
    % Temp: Temperature [deg C]
    % Cond: Conductivity [uS/cm] (standard for freshwaters)
    % 
    % O T P U T S
    % Cond20: Conductivity [uS/cm] at 20°C
    % alpha: temperature coefficient of variation [1/deg C]
    % dCdT: variation of C with T [uS/cm/deg C]
    
    A=1.684; B=-0.04645; C=0.000602;
    fT= A + B*Temp + C.*Temp.^2;
    Cond20=fT.*Cond; % C_20°C / Cond MUST be in uS/cm
    
    % comparing the 2nd order equation above with Cond20=Cond/(1+alpha(T-20)):
    Tmean=mean(Temp); Cmean=mean(Cond);
    alpha = -(A+B*Tmean+C*Tmean^2-1)/((Tmean-20)*(A+Tmean.*(B+C*Tmean)));
    dCdT = alpha*Cmean; % change of C with T
end