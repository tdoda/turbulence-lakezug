% Compare turbulence  profiles

close all
clear
clc

%% Parameters

savefig=false;
campaign_name='20250211';
datafolder=['..\..\data\VMP\',campaign_name,'\Level2\'];
folder_param_name='down1.0_K5.3_L_NAS0.0058_singlepole_nfft3_EPFL_Goodman';
filenames={'VMP003','VMP004','VMP007','VMP008','VMP009','VMP010','VMP012','VMP013','VMP014'}; % all files
profnames={'PSB01','PC01','PC02','PC03','PC04','PC05','PC06','PC07','PSB02'};

addpath("..\..\functions\figures\") % Add functions
addpath("..\..\functions\figures\export_fig-master")

%% Load the data
disp('Loading the data...')


for kf=1:length(filenames)
    datavar=load([datafolder,filenames{kf},'_',folder_param_name,'\results_',filenames{kf},'_down.mat']);
    datavar.filename=filenames{kf};
    DATA(kf)=datavar;

end

%% 1. CTD data
indprof=2:length(filenames);
col=winter(length(indprof));
col(end,:)=[0,0,0];
select_profile='longest';

fig_CTD=figure('Units','centimeters','Position',[1 1 18 10]);

dy=0; % To reduce height of plots
nrow=1; ncol=3;

ax1=subplot(nrow,ncol,1); % Temperature
hp=plot_prof_slow(DATA,filenames,filenames(indprof),col,select_profile,'temperature',1);
ylabel('Depth [m]')
xlabel('CTD-T [°C]')
set(ax1,'ydir','reverse')
hleg=legend(hp,profnames(indprof));
title(hleg,DATA(1).FAST{1}.date)
axpos=ax1.Position;
ax1.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax2=subplot(nrow,ncol,2); % Chl-a
plot_prof_slow(DATA,filenames,filenames(indprof),col,select_profile,'chlorophyll',50);
set(ax2,'ydir','reverse')
xlabel('Chl-a [ppb]')
axpos=ax2.Position;
ax2.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax3=subplot(nrow,ncol,3); % Turb
plot_prof_slow(DATA,filenames,filenames(indprof),col,select_profile,'turbidity',50);
set(ax3,'ydir','reverse')
xlabel('Turb [FTU]')
axpos=ax3.Position;
ax3.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

linkaxes([ax1,ax2,ax3],'y')

if savefig
    saveas(fig_CTD,['CTD_profiles_',campaign_name,'_',filenames{:},'.fig'])
    exportgraphics(fig_CTD,['CTD_profiles_',campaign_name,'_',filenames{:},'.png'],'Resolution',400)

end

%% 2. Fast sensors
indprof=[2,3,6,9];
col=winter(length(indprof));
col(end,:)=[0,0,0];
select_profile='longest';

fig_fast=figure('Units','centimeters','Position',[1 1 18 18]);

dy=0; % To reduce height of plots
nrow=2; ncol=4;

ax1=subplot(nrow,ncol,1); % T1
hp=plot_prof_fast(DATA,filenames,filenames(indprof),col,select_profile,'fast_T1',1);
ylabel('Depth [m]')
xlabel('T1 [°C]')
set(ax1,'ydir','reverse')
hleg=legend(hp,profnames(indprof));
title(hleg,DATA(1).FAST{1}.date)
axpos=ax1.Position;
ax1.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax2=subplot(nrow,ncol,ncol+1); % T2
plot_prof_fast(DATA,filenames,filenames(indprof),col,select_profile,'fast_T2',1);
set(ax2,'ydir','reverse')
xlabel('T2 [°C]')
ylabel('Depth [m]')
axpos=ax2.Position;
ax2.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax3=subplot(nrow,ncol,2); % gradT1
plot_prof_fast(DATA,filenames,filenames(indprof),col,select_profile,'grad_T1',1000,true);
% plot_prof_fast(DATA,filenames_all,filenames,col,select_profile,'grad_T1',1);
set(ax3,'ydir','reverse')
xlabel('|dT1/dz| [°C/m]')
set(gca,'xscale','log','xtick',10.^(-10:10))
axpos=ax3.Position;
ax3.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax4=subplot(nrow,ncol,ncol+2); % gradT2
plot_prof_fast(DATA,filenames,filenames(indprof),col,select_profile,'grad_T2',1000,true);
% plot_prof_fast(DATA,filenames_all,filenames,col,select_profile,'grad_T2',1);
set(ax4,'ydir','reverse')
xlabel('|dT2/dz| [°C/m]')
set(gca,'xscale','log','xtick',10.^(-10:10))
axpos=ax4.Position;
ax4.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax5=subplot(nrow,ncol,3); % S1
plot_prof_fast(DATA,filenames,filenames(indprof),col,select_profile,'fast_S1',1000,true);
% hp=plot_prof_fast(DATA,filenames_all,filenames,col,select_profile,'fast_S1',1);
xlabel('|S1| [s^{-1}]')
set(ax5,'ydir','reverse')
%xlim([-20 20])
set(gca,'xscale','log','xtick',10.^(-10:10))
axpos=ax5.Position;
ax5.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax6=subplot(nrow,ncol,ncol+3); % S2
plot_prof_fast(DATA,filenames,filenames(indprof),col,select_profile,'fast_S2',1000,true);
% plot_prof_fast(DATA,filenames_all,filenames,col,select_profile,'fast_S2',1);
set(ax6,'ydir','reverse')
xlabel('|S2| [s^{-1}]')
%xlim([-20 20])
set(gca,'xscale','log','xtick',10.^(-10:10))
axpos=ax6.Position;
ax6.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax7=subplot(nrow,ncol,4); % Ax
plot_prof_fast(DATA,filenames,filenames(indprof),col,select_profile,'A_x',1);
set(ax7,'ydir','reverse')
xlabel('Ax [m s^{-2}]')
axpos=ax7.Position;
ax7.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax8=subplot(nrow,ncol,ncol+4); % Ay
plot_prof_fast(DATA,filenames,filenames(indprof),col,select_profile,'A_y',1);
set(ax8,'ydir','reverse')
xlabel('Ay [m s^{-2}]')
axpos=ax8.Position;
ax8.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];


linkaxes([ax1,ax2,ax3,ax4,ax5,ax6],'y')

if savefig
    saveas(fig_fast,['fast_',campaign_name,'_',filenames{:},'.fig'])
    exportgraphics(fig_fast,['fast_',campaign_name,'_',filenames{:},'.png'],'Resolution',400)
end

%% 3a. Turbulence data
indprof=[2,3,6,9];
col=winter(length(indprof));
col(end,:)=[0,0,0];
select_profile='longest';

fig_turb=figure('Units','centimeters','Position',[1 1 18 18]);
dy=0; % To reduce height of plots
nrow=2; ncol=4;

flag_T1='flag_T1';
flag_T2='flag_T2';

ax1=subplot(nrow,ncol,1); % Eps-T1
hp=plot_prof_bin(DATA,filenames,filenames(indprof),col,select_profile,'eps_T1',1,flag_T1);
ylabel('Depth [m]')
xlabel('\epsilon_{T1} [m^2 s^{-3}]')
set(ax1,'ydir','reverse','xscale','log','xtick',10.^(-14:2:0),'XTickLabelRotation',45)
hleg=legend(hp,profnames(indprof));
title(hleg,DATA(1).FAST{1}.date)
axpos=ax1.Position;
ax1.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax5=subplot(nrow,ncol,ncol+1); % Eps-T2
plot_prof_bin(DATA,filenames,filenames(indprof),col,select_profile,'eps_T2',1,flag_T2);
xlabel('\epsilon_{T2} [m^2 s^{-3}]')
ylabel('Depth [m]')
set(ax5,'ydir','reverse','xscale','log','xtick',10.^(-14:2:0),'XTickLabelRotation',45)
xlim(ax1.XLim)
axpos=ax5.Position;
ax5.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

% ax9=subplot(nrow,ncol,2); % Eps-S1
% hp=plot_prof_ql(DATA,filenames_all,filenames,col,select_profile,'e',1,1);
% xlabel('\epsilon_{S1} [m^2 s^{-3}]')
% set(ax9,'ydir','reverse','xscale','log','xtick',10.^(-14:2:0),'XTickLabelRotation',45)
% axpos=ax9.Position;
% ax9.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];
% xlim(ax1.XLim)
% 
% ax10=subplot(nrow,ncol,ncol+2); % Eps-S2
% plot_prof_ql(DATA,filenames_all,filenames,col,select_profile,'e',2,1);
% xlabel('\epsilon_{S2} [m^2 s^{-3}]')
% set(ax10,'ydir','reverse','xscale','log','xtick',10.^(-14:2:0),'XTickLabelRotation',45)
% xlim(ax1.XLim)
% axpos=ax10.Position;
% ax10.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax2=subplot(nrow,ncol,2); % X-T1
plot_prof_bin(DATA,filenames,filenames(indprof),col,select_profile,'Xi_T1',1,flag_T1);
xlabel('\chi_{T1} [K^2 s^{-1}]')
set(ax2,'ydir','reverse','xscale','log','xtick',10.^(-14:2:0),'XTickLabelRotation',45)
axpos=ax2.Position;
ax2.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax6=subplot(nrow,ncol,ncol+2); % X-T2
plot_prof_bin(DATA,filenames,filenames(indprof),col,select_profile,'Xi_T2',1,flag_T2);
xlabel('\chi_{T2} [K^2 s^{-1}]')
set(ax6,'ydir','reverse','xscale','log','xtick',10.^(-14:2:0),'XTickLabelRotation',45)
axpos=ax6.Position;
xlim(ax2.XLim)
ax6.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax3=subplot(nrow,ncol,3); % Kz-T1
plot_prof_bin(DATA,filenames,filenames(indprof),col,select_profile,'KOsbornCox_T1',1,flag_T1);
xlabel('K_{z,T1} [m^2 s^{-1}]')
set(ax3,'ydir','reverse','xscale','log','xtick',10.^(-14:2:0),'XTickLabelRotation',45)
axpos=ax3.Position;
ax3.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax7=subplot(nrow,ncol,ncol+3); % Kz-T2
plot_prof_bin(DATA,filenames,filenames(indprof),col,select_profile,'KOsbornCox_T2',1,flag_T2);
xlabel('K_{z,T2} [m^2 s^{-1}]')
set(ax7,'ydir','reverse','xscale','log','xtick',10.^(-14:2:0),'XTickLabelRotation',45)
axpos=ax7.Position;
xlim(ax3.XLim)
ax7.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax4=subplot(nrow,ncol,4); % LT1
plot_prof_bin(DATA,filenames,filenames(indprof),col,select_profile,'LTuT1',1,''); % No flagging
xlabel('L_{T,T1} [m]')
set(ax4,'ydir','reverse','xscale','log','xtick',10.^(-14:2:0),'XTickLabelRotation',45)
axpos=ax4.Position;
ax4.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax8=subplot(nrow,ncol,ncol+4); % LT2
plot_prof_bin(DATA,filenames,filenames(indprof),col,select_profile,'LTuT2',1,''); % No flagging
xlabel('L_{T,T2} [m]')
set(ax8,'ydir','reverse','xscale','log','xtick',10.^(-14:2:0),'XTickLabelRotation',45)
axpos=ax8.Position;
xlim(ax4.XLim)
ax8.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8],'y')

if savefig
    saveas(fig_turb,['turbulence_',campaign_name,'_',filenames{:},'.fig'])
    exportgraphics(fig_turb,['turbulence_',campaign_name,'_',filenames{:},'.png'],'Resolution',400)
end

%% 3b. Turbulence T1 data with median
indprof=[2,3,6,9];
col=winter(length(indprof));
col(end,:)=[0,0,0];
select_profile='longest';

fig_turbmed=figure('Units','centimeters','Position',[1 1 18 9]);
dy=0.1; % To reduce height of plots
nrow=1; ncol=3;

ax1=subplot(nrow,ncol,1); % Eps-T1
[depthC_all,epsC_all,depthC_cells,logmed_epsC] = compute_logmedian(DATA,5,'eps_T1',[2,3,6]);
[depthS_all,epsS_all,depthS_cells,logmed_epsS] = compute_logmedian(DATA,5,'eps_T1',9);
plot(epsC_all,depthC_all,'o','markersize',1,'markeredgecolor',[1 0.75 0.6],'markerfacecolor',[1 0.75 0.6])
hold on
plot(10.^(logmed_epsC),depthC_cells,'-','color',[1 0.4 0.15])
plot(epsS_all,depthS_all,'o','markersize',1,'markeredgecolor',[0.7 0.7 0.7],'markerfacecolor',[0.7 0.7 0.7])
plot(10.^(logmed_epsS),depthS_cells,'k-')
ylabel('Depth [m]')
xlabel('\epsilon [m^2 s^{-3}]')
set(ax1,'ydir','reverse','xscale','log','xtick',10.^(-14:2:0),'XTickLabelRotation',45)
axpos=ax1.Position;
ax1.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];
xlim([10^(-11),10^(-4)])

ax2=subplot(nrow,ncol,2); % X-T1
[depthC_all,XiC_all,depthC_cells,logmed_XiC] = compute_logmedian(DATA,5,'Xi_T1',[2,3,6]);
[depthS_all,XiS_all,depthS_cells,logmed_XiS] = compute_logmedian(DATA,5,'Xi_T1',9);
plot(XiC_all,depthC_all,'o','markersize',1,'markeredgecolor',[1 0.75 0.6],'markerfacecolor',[1 0.75 0.6])
hold on
plot(10.^(logmed_XiC),depthC_cells,'-','color',[1 0.4 0.15])
plot(XiS_all,depthS_all,'o','markersize',1,'markeredgecolor',[0.7 0.7 0.7],'markerfacecolor',[0.7 0.7 0.7])
plot(10.^(logmed_XiS),depthS_cells,'k-')
ylabel('Depth [m]')
xlabel('\chi [K^2 s^{-1}]')
set(ax2,'ydir','reverse','xscale','log','xtick',10.^(-14:2:0),'XTickLabelRotation',45)
axpos=ax2.Position;
ax2.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];
xlim([10^(-13),10^(-6)])

ax3=subplot(nrow,ncol,3); % Kz-T1
[depthC_all,KzC_all,depthC_cells,logmed_KzC] = compute_logmedian(DATA,5,'KOsbornCox_T1',[2,3,6]);
[depthS_all,KzS_all,depthS_cells,logmed_KzS] = compute_logmedian(DATA,5,'KOsbornCox_T1',9);
plot(KzC_all,depthC_all,'o','markersize',1,'markeredgecolor',[1 0.75 0.6],'markerfacecolor',[1 0.75 0.6])
hold on
plot(10.^(logmed_KzC),depthC_cells,'-','color',[1 0.4 0.15])
plot(KzS_all,depthS_all,'o','markersize',1,'markeredgecolor',[0.7 0.7 0.7],'markerfacecolor',[0.7 0.7 0.7])
plot(10.^(logmed_KzS),depthS_cells,'k-')
ylabel('Depth [m]')
xlabel('K_{z} [m^2 s^{-1}]')
set(ax3,'ydir','reverse','xscale','log','xtick',10.^(-14:2:0),'XTickLabelRotation',45)
axpos=ax3.Position;
ax3.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];
%xlim([10^(-11),10^(-4)])

linkaxes([ax1,ax2,ax3],'y')

if savefig
    saveas(fig_turbmed,['turbulence_median_',campaign_name,'_',filenames{:},'.fig'])
    exportgraphics(fig_turbmed,['turbulence_median_',campaign_name,'_',filenames{:},'.png'],'Resolution',400)
    export_fig(fig_turbmed,['turbulence_median_',campaign_name,'_',filenames{:},'.eps'],'-transparent','-painters','-nocrop')
end

%% 4. Thorpe length scale
% indprof=[2,3,6,9];
indprof=[3,9];
col=winter(length(indprof));
col(end,:)=[0,0,0];
select_profile='longest';
ylimval=[];
Tlimval=[];
% Segment 1:
% ylimval=[63,75];
% Tlimval=[5.73,5.83];
% Zoom 1:
% ylimval=[67,72];
% Tlimval=[5.76,5.79];
% Segment 1:
% ylimval=[120,140];
% Tlimval=[5.36,5.41];


fig_Thorpe=figure('Units','centimeters','Position',[1 1 9 18]);
dy=0; % To reduce height of plots
nrow=2; ncol=2;

flag_T1='flag_T1';
flag_T2='flag_T2';

ax1=subplot(nrow,ncol,1); % T1
hp=plot_prof_fast(DATA,filenames,filenames(indprof),col,select_profile,'fast_T1',1);
ylabel('Depth [m]')
xlabel('T1 [°C]')
set(ax1,'ydir','reverse')
hleg=legend(hp,profnames(indprof));
title(hleg,DATA(1).FAST{1}.date)
axpos=ax1.Position;
ax1.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax2=subplot(nrow,ncol,ncol+1); % T2
plot_prof_fast(DATA,filenames,filenames(indprof),col,select_profile,'fast_T2',1);
set(ax2,'ydir','reverse')
xlabel('T2 [°C]')
ylabel('Depth [m]')
axpos=ax2.Position;
ax2.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax3=subplot(nrow,ncol,2); % LT1
plot_prof_bin(DATA,filenames,filenames(indprof),col,select_profile,'LTuT1',1,''); % No flagging
xlabel('L_{T,T1} [m]')
set(ax3,'ydir','reverse','xscale','log','xtick',10.^(-14:2:0),'XTickLabelRotation',45)
axpos=ax3.Position;
ax3.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax4=subplot(nrow,ncol,ncol+2); % LT2
plot_prof_bin(DATA,filenames,filenames(indprof),col,select_profile,'LTuT2',1,''); % No flagging
xlabel('L_{T,T2} [m]')
set(ax4,'ydir','reverse','xscale','log','xtick',10.^(-14:2:0),'XTickLabelRotation',45)
axpos=ax4.Position;
ax4.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];
xlim(ax3.XLim)


linkaxes([ax1,ax2,ax3,ax4],'y')
if ~isempty(ylimval)
    set(ax1,'Ylim',ylimval)
end

if ~isempty(Tlimval)
    set(ax1,'Xlim',Tlimval)
    set(ax2,'Xlim',Tlimval)
end


if savefig
    saveas(fig_Thorpe,['Thorpe_',campaign_name,'_',filenames{:},'.fig'])
    exportgraphics(fig_Thorpe,['Thorpe_',campaign_name,'_',filenames{:},'.png'],'Resolution',400)
end

%% Figures turbulence
% indprof=[2,3,6,9];
% col_day=winter(length(indprof));
% col_day(end,:)=[0,0,0];
% [fig1,fig2,fig3,fig4]=plot_comparison_turbulence(DATA(indprof),filenames(indprof),'longest',col_day);
% 
% if savefig
%     saveas(fig1,['CTD_profiles_',campaign_name,'_',filenames{:},'.fig'])
%     exportgraphics(fig1,['CTD_profiles_',campaign_name,'_',filenames{:},'.png'],'Resolution',400)
% 
%     saveas(fig2,['fast_profiles_',campaign_name,'_',filenames{:},'.fig'])
%     exportgraphics(fig2,['fast_profiles_',campaign_name,'_',filenames{:},'.png'],'Resolution',400)
% 
%     saveas(fig3,['turbulence_profiles_',campaign_name,'_',filenames{:},'.fig'])
%     exportgraphics(fig3,['turbulence_profiles_',campaign_name,'_',filenames{:},'.png'],'Resolution',400)
% end

