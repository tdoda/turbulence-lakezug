function [fig1] = plot_spectrum_FP07(data_prof_bin,data_prof_fast,ind_bin)
%PLOT_SPECTRUM_FP07 Plot FP07 temperature spectrum from a bin.
%   Detailed explanation goes here

fig1=figure('Units', 'Centimeters', 'Position', [1,1,18,8]);
s_raw=loglog(k0,PSD_raw, '-k','linewidth',0.5); hold on
s_corr=loglog(k0,PSD0,'-','color',[75 97 209]/255,'linewidth',0.75);
if cont
    loglog(k(ikfit),PSD(ikfit),'ok','markerfacecolor',[75 97 209]/255, 'markersize',5)
    s_fit=loglog(k0,Tspec(Tdis, Xi_T, kB,k0,D,q), 'color',[75 97 209]/255, 'linewidth',2);
    powerfit=loglog(k(ikfit),Sm(ikfit)+Sn(ikfit),'-','color',[174  38  43]/255,'linewidth',2);
end
s_noise=loglog(k0,Sn0, ':k', 'linewidth',0.7);
ylim([10^-9,1e3])

text(1.2,1e-4,['\chi_{T} = ', num2str(Xi_T,'%1.3e'), '^{\circ}C^2 s^{-1}'],'horizontalalignment','left','fontsize',8)
text(1.2,10^-4.9,['k_B = ', num2str(kB,'%3.0f'),' cpm'],'horizontalalignment','left','fontsize',8)
text(1.2,10^-5.7,['\epsilon_T = ', num2str(eps_T,'%1.2e'),' m^2 s^{-3}'],'horizontalalignment','left','fontsize',8)
text(1.2,10^-6.6,['LR = ', num2str(LR,'%1.1f'),'(2)'],'horizontalalignment','left','fontsize',8)
text(1.2,10^-7.45,['MAD = ', num2str(MAD_T,'%1.2f'),' (', num2str(2*MADc, '%1.2f'),')'],'horizontalalignment','left','fontsize',8)

plot(kB,1.4E-9,'^','markeredgecolor','k','markerfacecolor','k','markersize',5)
plot(fn/W,1.4E-9,'^','markeredgecolor',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5],'markersize',5)
itmp=find(PSD0<Snfact*Sn0 ,1,'first');
if isempty(itmp)
    itmp=length(k0);
end
plot(k0(itmp),1.4E-9,'^','markerfacecolor',[0.75 0.75 0.75],'markeredgecolor',[0.75 0.75 0.75],'markersize',5)
plot(k0(find(H<H_lim ,1,'first')),1.4E-9,'^','markeredgecolor','k','markerfacecolor',[1 1 1],'markersize',5)
plot(kP,1.4E-9,'^','markeredgecolor',[75 97 209]/255,'markerfacecolor',[75 97 209]/255,'markersize',5)
xlabel('k [cpm]','fontsize',10)
ylabel('\Psi_T [(^{\circ}C m^{-1})^{-2} cpm^{-1}]','fontsize',10)
xlim([1,1000])
set(gca,'ytick',[1e-9 1e-7 1e-5 1e-3 1e-1 1e1 1e3]);
legend([s_raw,s_corr,s_fit,powerfit,s_noise],{'Raw spec.','Corrected spec.','Kraichnan spec. (fit)','Power law (fit)','Noise spec.'},...
    'location','North','orientation','horizontal','NumColumns',2)
grid on; ax.GridAlpha=0.1;
set(gca, 'xminorgrid', 'on')
set(gca,'MinorGridLineStyle','-','MinorGridAlpha',0.1)
set(gca,'GridLineStyle','-')

end