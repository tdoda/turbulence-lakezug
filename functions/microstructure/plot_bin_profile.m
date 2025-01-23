function [] = plot_bin_profile(SLOW,BIN,param)
%PLOT_BIN_PROFILE Make plots of the binned profiles.
%   Detailed explanation goes here

%%
info=param.info;

pmin = info.pmin - info.dpD/2;
pmax = BIN.pressure(find(isfinite(BIN.temperature),1,'last')) + info.dpD/2;

figure('Units','centimeters','Position',[1 1 18 18])
clf
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperSize', [29 20]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition', [0 0 29 20]);

% 1. Temperature and salinity profiles
ax1=subplot(2,4,1);
plot(SLOW.temperature, SLOW.pressure,'-k','linewidth',1)
xlabel('T (°C)')
ylabel('p (db)')
set(gca,'YDir','reverse')
ylim([pmin,pmax])
ax2=axes('Position',ax1.Position);
set(ax1,'box','off')
plot(SLOW.salinity/1000, SLOW.pressure,'r', 'parent' , ax2,'linewidth',1)
set(ax2,'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'XColor','r','YColor','k');
yticklabels([])
ylim([pmin,pmax])
xlabel('Salinity (g/kg)')
set(gca,'yticklabel',[])
set(gca,'YDir','reverse')

% 2. Speed, Chl-a & turbidity in bins
ax3=subplot(2,4,2);
p1=plot(BIN.speed,BIN.pressure,'.-k','linewidth',1,'markersize',4);
xlabel('W (db/s)')
ylim([pmin,pmax])
set(gca,'yticklabel',[])
set(gca,'YDir','reverse')
grid('on')
ax4=axes('Position',ax3.Position);
set(ax3,'box','off')
p2=plot(BIN.chlorophyll, BIN.pressure,'g', 'parent' , ax4,'linewidth',1);
hold on
p3=plot(BIN.turbidity, BIN.pressure,'r', 'parent' , ax4,'linewidth',1);
set(ax4,'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'XColor','g','YColor','k');
xlabel(sprintf('Chl [%s]-Turb [%s]',param.unit_Chl,param.unit_Turb),'color','k')
yticklabels([])
ylim([pmin,pmax])
set(gca,'yticklabel',[])
set(gca,'YDir','reverse')
legend([p1,p2,p3],{'W','Chl-a','Turb'})

% 3. Epsilon (T1, T2, S1, S2) in bins
subplot(2,4,3)
leg_entries={};
if param.config.T1
    plot(BIN.eps_T1,BIN.pressure,'.-','linewidth',1,'markersize',4)
    hold on
    plot(BIN.eps_T1(BIN.flag_T1==0),BIN.pressure(BIN.flag_T1==0),'ob','markersize',3)
    leg_entries(end+1:end+2)={'T01','T01_flagged'};
end

if param.config.T2
    plot(BIN.eps_T2,BIN.pressure,'.-','linewidth',1,'markersize',4)
    hold on
    plot(BIN.eps_T2(BIN.flag_T2==0),BIN.pressure(BIN.flag_T2==0),'or','markersize',3)
    leg_entries(end+1:end+2)={'T02','T02_flagged'};
end

if param.config.S1
    plot(BIN.eps_S1,BIN.pressure,'.-','linewidth',1,'markersize',4)
    hold on
    leg_entries{end+1}='sh1';
end

if param.config.S2
    plot(BIN.eps_S2,BIN.pressure,'.-','linewidth',1,'markersize',4)
    hold on
    leg_entries{end+1}='sh2';
end

if param.config.T1 && param.config.T2
    plot(0.5*(BIN.epsT1max+BIN.epsT2max),BIN.pressure,'-','linewidth',1,'markersize',4,'color',[0.5,0.5,0.5])
    leg_entries{end+1}='maxT';
end

leg=legend(leg_entries,'fontsize',7,'location','best');
leg.ItemTokenSize = [7,7];
xlabel('\epsilon (m^2 s^{-3})')
set(gca,'yticklabel',[])
set(gca,'YDir','reverse')
set(gca, 'xscale','log')
xlim([1e-12,1e-4])
ylim([pmin,pmax])
set(gca,'xtick',10.^(-12:2:-4))
grid('on')
title(datestr(SLOW.datenum,'yyyy-mm-dd HH:MM:SS'))

% 4. Xi (T1, T2, ST1, ST2) in bins
subplot(2,4,4)
leg_entries={};
if param.config.T1
    p1=plot(BIN.Xi_T1,BIN.pressure,'.-','linewidth',1,'markersize',4);
    hold on
    plot(BIN.Xi_ST1,BIN.pressure,'--','linewidth',1,'color',p1.Color)
    leg_entries(end+1:end+2)={'T01','ST01'};
end

if param.config.T2
    p2=plot(BIN.Xi_T2,BIN.pressure,'.-','linewidth',1,'markersize',4);
    hold on
    plot(BIN.Xi_ST2,BIN.pressure,'--','linewidth',1,'color',p2.Color)
    leg_entries(end+1:end+2)={'T02','ST02'};
end
leg=legend(leg_entries,'fontsize',7,'location','best');
leg.ItemTokenSize = [7,7];
xlabel('\chi (K^2 s^{-1})')
set(gca,'yticklabel',[])
set(gca,'YDir','reverse')
set(gca,'xscale','log')
xlim([1e-12,1e-3])
ylim([pmin,pmax])
set(gca,'xtick',10.^(-11:2:-3))
grid('on')

% 5. Kz (T1, T2, S1, S2) in bins
subplot(2,4,5)
leg_entries={};
if param.config.T1
    plot(BIN.KOsbornCox_T1,BIN.pressure,'.-','linewidth',1,'markersize',4)
    hold on
    leg_entries{end+1}='T01 (O&C)';
end

if param.config.T2
    plot(BIN.KOsbornCox_T2,BIN.pressure,'.-','linewidth',1,'markersize',4)
    hold on
    leg_entries{end+1}='T02 (O&C)';
end

if param.config.S1
    plot(BIN.KOsborn_S1,BIN.pressure,'.-','linewidth',1,'markersize',4)
    hold on
    leg_entries{end+1}='sh1 (O)';
end

if param.config.S2
    plot(BIN.KOsborn_S2,BIN.pressure,'.-','linewidth',1,'markersize',4)
    hold on
    leg_entries{end+1}='sh2 (O)';
end
xlabel('K (m^2 s^{-1})')
ylabel('p (db)')
set(gca,'YDir','reverse')
set(gca,'xscale','log')
xlim([1e-9,1e0])
%line([D,D],ylim(),'color',[0.5,0.5,0.5])
plot(BIN.Diff,BIN.pressure,'color',[0.5,0.5,0.5]) % Plot molecular diffusivity profile
leg_entries{end+1}='molecular';
leg=legend(leg_entries,'fontsize',7,'location','best');
leg.ItemTokenSize = [7,7];
ylim([pmin,pmax])
set(gca,'xtick',10.^(-7:2:0))
grid('on')

% 6. Thorpe and Ozmidov scales in bins
subplot(2,4,6)
% $$$     l1=plot(BIN.kB_T1,BIN.pressure,'.-','linewidth',1,'markersize',4);
% $$$     hold on
% $$$     l2=plot(BIN.kB_T2,BIN.pressure,'.-','linewidth',1,'markersize',4);
% $$$     plot(BIN.KBSH,BIN.pressure,'.-','linewidth',1,'markersize',4)
% $$$     plot(BIN.kU_T1,BIN.pressure,'--','linewidth',1,'color',l1.Color())
% $$$     plot(BIN.kU_T2,BIN.pressure,'--','linewidth',1,'color',l2.Color())
% $$$     xlabel('K_B (m)')
% $$$     set(gca,'YDir','reverse')
% $$$     %set(gca,'xscale','log')
% $$$     grid('on')
% $$$     ylim([pmin,pmax])

leg_LT={};

if param.config.T1
    plot(BIN.LTuT1,BIN.pressure,'.-','linewidth',1,'markersize',4)
    hold on
    leg_LT{end+1}='L_T^{T1}';
end

if param.config.T2
    plot(BIN.LTuT2,BIN.pressure,'.-','linewidth',1,'markersize',4)
    hold on
    leg_LT{end+1}='L_T^{T2}';
end

if param.config.S1
    plot(BIN.LO_S1,BIN.pressure,'.-','linewidth',1,'markersize',4)
    hold on
    leg_LT{end+1}='L_O^{sh1}';
end

if param.config.S2
    plot(BIN.LO_S2,BIN.pressure,'.-','linewidth',1,'markersize',4)
    hold on
    leg_LT{end+1}='L_O^{sh2}';
end

leg=legend(leg_LT,'fontsize',7,'location','best');
leg.ItemTokenSize = [7,7];
xlabel('L_T, L_O (m)')
set(gca,'YDir','reverse')
set(gca,'xscale','log')
grid('on')
ylim([pmin,pmax])


% 7. Performance estimates: Mean Absolute Deviation
subplot(2,4,7)
if param.config.T1
    plot(BIN.MAD_T1, BIN.pressure,'.-','linewidth',1,'markersize',4)
    hold on
end

if param.config.T2
    plot(BIN.MAD_T2, BIN.pressure,'.-','linewidth',1,'markersize',4)
    hold on
end

if param.config.S1
    plot(BIN.MAD_S1, BIN.pressure,'.-','linewidth',1,'markersize',4)
    hold on
end

if param.config.S2
    plot(BIN.MAD_S2, BIN.pressure,'.-','linewidth',1,'markersize',4)
    hold on
end
line(2*[min(BIN.MADc),min(BIN.MADc)],ylim, 'color','k')
line([min(BIN.MADc),min(BIN.MADc)],ylim, 'color','k')
xlabel('MAD')
xlim([0,2])
set(gca,'yticklabel',[])
set(gca,'YDir','reverse')
grid('on')
ylim([pmin,pmax])

% 8. Performance estimates: Likelihood Ratio
subplot(2,4,8)
if param.config.T1
    plot(BIN.LR_T1, BIN.pressure,'.-','linewidth',1,'markersize',4)
    hold on
end

if param.config.T2
    plot(BIN.LR_T2, BIN.pressure,'.-','linewidth',1,'markersize',4)
    hold on
end
xlabel('Likelihood ratio')
set(gca,'yticklabel',[])
set(gca,'YDir','reverse')
grid('on')
ylim([pmin,pmax])




%hold on
%semilogx(epsilon2,-pres)
%semilogx(epsilonN,-pres)


end