function [Cond20,alpha,dCdT]=conductivity_20(Temp,Cond)
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
Cond20=fT.*Cond; % C_20�C / Cond MUST be in uS/cm

% comparing the 2nd order equation above with Cond20=Cond/(1+alpha(T-20)):
Tmean=mean(Temp); Cmean=mean(Cond);
alpha = -(A+B*Tmean+C*Tmean^2-1)/((Tmean-20)*(A+Tmean.*(B+C*Tmean)));
dCdT = alpha*Cmean; % change of C with T