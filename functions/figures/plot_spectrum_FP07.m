function [fig1] = plot_spectrum_FP07(data_prof_bin,data_spectrum,param,ind_bin,sensor_name)
%PLOT_SPECTRUM_FP07 Plot FP07 temperature spectrum from a bin.
%   Detailed explanation goes here

%% Parameters
SPECTRUM=data_spectrum(ind_bin); % Spectral parameters
ikfit=SPECTRUM.ind_fit;

% Get parameters from the bin
Xi_T=data_prof_bin.(['Xi_',sensor_name])(ind_bin);
kB=data_prof_bin.(['kB_',sensor_name])(ind_bin);
eps_T=data_prof_bin.(['eps_',sensor_name])(ind_bin);
LR=data_prof_bin.(['LR_',sensor_name])(ind_bin);
MAD_T=data_prof_bin.(['MAD_',sensor_name])(ind_bin);
kP=data_prof_bin.(['kpeak_',sensor_name])(ind_bin);
MADc=data_prof_bin.MADc(ind_bin); % Threshold for the Mean Absolute Deviation between observed and empirical spectra (independent of the sensor used)

% Other parameters
fn=param.info.fAA;
W=data_prof_bin.WW(ind_bin);
Snfact=param.info.Snfact;

%% Plot

fig1=figure('Units', 'Centimeters', 'Position', [1,1,18,8]);

% 1. Raw PSD
s_raw=loglog(SPECTRUM.k0,SPECTRUM.PSD_raw, '-k','linewidth',0.5); hold on

% 2. Corrected PSD for sensor time response
s_corr=loglog(SPECTRUM.k0,SPECTRUM.PSD0,'-','color',[75 97 209]/255,'linewidth',0.75);

% 3. Corrected PSD over the wavenumber range used for fitting
loglog(SPECTRUM.k(ikfit),SPECTRUM.PSD(ikfit),'ok','markerfacecolor',[75 97 209]/255, 'markersize',5)

% 3. Full theoretical spectrum 
s_fit=loglog(SPECTRUM.k0,SPECTRUM.PSD_theo, 'color',[75 97 209]/255, 'linewidth',2);

% 4. Power fitted spectrum (including noise spectrum) used to compute the
% likelihood ratio (ratio between the likelihood from theoretical spectrum
% and likelihood power-law fitting):
powerfit=loglog(SPECTRUM.k(ikfit),SPECTRUM.Sm(ikfit)+SPECTRUM.Sn(ikfit),'-','color',[174  38  43]/255,'linewidth',2);

% 5. Noise spectrum corrected for time response (all wavenumbers):
s_noise=loglog(SPECTRUM.k0,SPECTRUM.Sn0, ':k', 'linewidth',0.7);

ylim([10^-9,1e3])

% Display values of estimated parameters:
text(1.2,1e-4,['\chi_{T} = ', num2str(Xi_T,'%1.3e'), '^{\circ}C^2 s^{-1}'],'horizontalalignment','left','fontsize',8)
text(1.2,10^-4.9,['k_B = ', num2str(kB,'%3.0f'),' cpm'],'horizontalalignment','left','fontsize',8)
text(1.2,10^-5.7,['\epsilon_T = ', num2str(eps_T,'%1.2e'),' m^2 s^{-3}'],'horizontalalignment','left','fontsize',8)
text(1.2,10^-6.6,['LR = ', num2str(LR,'%1.1f'),'(2)'],'horizontalalignment','left','fontsize',8) % Likelihood ratio as Log10
text(1.2,10^-7.45,['MAD = ', num2str(MAD_T,'%1.2f'),' (', num2str(2*MADc, '%1.2f'),')'],'horizontalalignment','left','fontsize',8)

% Indicate Batchelor wavenumber:
plot(kB,1.4E-9,'^','markeredgecolor',[75 97 209]/255,'markerfacecolor',[75 97 209]/255,'markersize',5)
text(kB,0.8E-9,'k_B','HorizontalAlignment','center','VerticalAlignment','top','color',[75 97 209]/255)

% Indicate upper wavenumber limit based on the anti-aliasing filter:
plot(fn/W,1.4E-9,'>','markeredgecolor',[0.8 0 0],'markerfacecolor',[0.8 0 0],'markersize',5)
text(fn/W,0.8E-9,'k_{AA}','HorizontalAlignment','center','VerticalAlignment','top','color',[0.8 0 0])

% Indicate wavenumber above which noise dominates (SNR<Snfact):
itmp=find(SPECTRUM.PSD0<Snfact*SPECTRUM.Sn0 ,1,'first');
if isempty(itmp)
    itmp=length(SPECTRUM.k0);
end
plot(SPECTRUM.k0(itmp),1.4E-9,'<','markerfacecolor','k','markeredgecolor','k','markersize',5)
text(SPECTRUM.k0(itmp),1.4E-9,'k_{SNR}','HorizontalAlignment','center','VerticalAlignment','bottom')

% Indicate wavenumber above which H<H_lim:
plot(SPECTRUM.k0(find(SPECTRUM.H<SPECTRUM.H_lim ,1,'first')),1.4E-9,'s','markeredgecolor',[0.75 0.75 0.75],'markerfacecolor',[0.75 0.75 0.75],'markersize',5)
text(SPECTRUM.k0(find(SPECTRUM.H<SPECTRUM.H_lim ,1,'first')),0.8E-9,'k_{H}','HorizontalAlignment','center','VerticalAlignment','top','color',[0.75 0.75 0.75])

% Indicate peak wavenumber:
plot(kP,1.4E-9,'o','markeredgecolor',[75 97 209]/255,'markerfacecolor',[75 97 209]/255,'markersize',5)
text(kP,0.8E-9,'k_{P}','HorizontalAlignment','center','VerticalAlignment','top','color',[75 97 209]/255)

xlabel('k [cpm]','fontsize',10)
ylabel('\Psi_T [(^{\circ}C m^{-1})^{-2} cpm^{-1}]','fontsize',10)
xlim([1,1000])
set(gca,'ytick',[1e-9 1e-7 1e-5 1e-3 1e-1 1e1 1e3]);
legend([s_raw,s_corr,s_fit,powerfit,s_noise],{'Raw spec.','Corrected spec.','Kraichnan spec. (fit)','Power law (fit)','Noise spec.'},...
    'location','North','orientation','horizontal','NumColumns',2)
grid on; 
set(gca,'GridAlpha',0.1);
set(gca, 'xminorgrid', 'on')
set(gca,'MinorGridLineStyle','-','MinorGridAlpha',0.1)
set(gca,'GridLineStyle','-')

title(sprintf('Bin %d/%d (%0.2f m), %s',ind_bin,length(data_spectrum),data_prof_bin.depth(ind_bin),sensor_name))

end