function [fig_CTD,fig_fast,fig_turb_all,fig_turbS] = plot_comparison_turbulence(DATA,filenames,select_profile,col,leg)

% select_profile is optional: ='all' (all profiles from the file, default) or 'longest' (longest
% profile of the file)
% col is optional
% leg is optional (filename by default)


%%
filenames_all={DATA(:).filename};

% Compute number of profiles
nprof=0;
profname={};
for kf=1:length(filenames) % Loop over files
    ind_file=find(strcmp(filenames_all,filenames{kf}),1);
    nprof=nprof+length(DATA(ind_file).BINNED);
    for kprof=1:length(DATA(ind_file).BINNED)
        profname{end+1}=[filenames{kf},'-P',num2str(kprof)];
    end
end

if nargin<3
    select_profile='all';
end

if strcmp(select_profile,'all')
    if nargin<4
        col=lines(nprof);
    elseif size(col,1)<nprof
        error('Not enough colors are specified.')
    end

    if nargin<5
        leg=profname;
    elseif length(leg)<nprof
        error('Not enough legend entries are specified.')
    end
else
    if nargin<4
        col=lines(length(filenames));
    end

    if nargin<5
        leg=filenames;
    end

end

%% 1. CTD data
fig_CTD=figure('Units','centimeters','Position',[1 1 14 18]);

dy=0; % To reduce height of plots
nrow=2; ncol=3;

ax1=subplot(nrow,ncol,1); % Eps-T1
hp=plot_prof_slow(DATA,filenames_all,filenames,col,select_profile,'temperature',1);
ylabel('Depth [m]')
xlabel('CTD-T [°C]')
set(ax1,'ydir','reverse')
hleg=legend(hp,leg);
title(hleg,DATA(1).FAST{1}.date)
axpos=ax1.Position;
ax1.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

%subplot(nrow,ncol,1) % Accurate S: didn't work for Lake Zug

ax2=subplot(nrow,ncol,2); % Chl-a
plot_prof_slow(DATA,filenames_all,filenames,col,select_profile,'chlorophyll',1);
set(ax2,'ydir','reverse')
xlabel('Chl-a [ppb]')
axpos=ax2.Position;
ax2.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax3=subplot(nrow,ncol,3); % Turb
plot_prof_slow(DATA,filenames_all,filenames,col,select_profile,'turbidity',1);
set(ax3,'ydir','reverse')
xlabel('Turb [FTU]')
axpos=ax3.Position;
ax3.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax4=subplot(nrow,ncol,4); % velocity
plot_prof_slow(DATA,filenames_all,filenames,col,select_profile,'velocity',1);
set(ax4,'ydir','reverse')
xlabel('W [m/s]')
axpos=ax4.Position;
ax4.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax5=subplot(nrow,ncol,5); % Inclination X
plot_prof_slow(DATA,filenames_all,filenames,col,select_profile,'Incl_x',1);
set(ax5,'ydir','reverse')
xlabel('Incl X [°]')
axpos=ax5.Position;
ax5.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];


ax6=subplot(nrow,ncol,6); % Inclination Y
plot_prof_slow(DATA,filenames_all,filenames,col,select_profile,'Incl_y',1);
set(ax6,'ydir','reverse')
xlabel('Incl Y [°]')
axpos=ax6.Position;
ax6.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];



linkaxes([ax1,ax2,ax3,ax4,ax5,ax6],'y')

%% 2. Fast temperature & shear
fig_fast=figure('Units','centimeters','Position',[1 1 18 18]);

dy=0; % To reduce height of plots
nrow=2; ncol=4;

ax1=subplot(nrow,ncol,1); % T1
hp=plot_prof_fast(DATA,filenames_all,filenames,col,select_profile,'fast_T1',1);
ylabel('Depth [m]')
xlabel('T1 [°C]')
set(ax1,'ydir','reverse')
hleg=legend(hp,leg);
title(hleg,DATA(1).FAST{1}.date)
axpos=ax1.Position;
ax1.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax2=subplot(nrow,ncol,5); % T2
plot_prof_fast(DATA,filenames_all,filenames,col,select_profile,'fast_T2',1);
set(ax2,'ydir','reverse')
xlabel('T2 [°C]')
axpos=ax2.Position;
ax2.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax3=subplot(nrow,ncol,2); % gradT1
plot_prof_fast(DATA,filenames_all,filenames,col,select_profile,'grad_T1',1000,true);
% plot_prof_fast(DATA,filenames_all,filenames,col,select_profile,'grad_T1',1);
set(ax3,'ydir','reverse')
xlabel('|dT1/dz| [°C/m]')
set(gca,'xscale','log','xtick',10.^(-10:10))
axpos=ax3.Position;
ax3.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax4=subplot(nrow,ncol,6); % gradT2
plot_prof_fast(DATA,filenames_all,filenames,col,select_profile,'grad_T2',1000,true);
% plot_prof_fast(DATA,filenames_all,filenames,col,select_profile,'grad_T2',1);
set(ax4,'ydir','reverse')
xlabel('|dT2/dz| [°C/m]')
set(gca,'xscale','log','xtick',10.^(-10:10))
axpos=ax4.Position;
ax4.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax5=subplot(nrow,ncol,3); % S1
plot_prof_fast(DATA,filenames_all,filenames,col,select_profile,'fast_S1',1000,true);
% hp=plot_prof_fast(DATA,filenames_all,filenames,col,select_profile,'fast_S1',1);
ylabel('Depth [m]')
xlabel('|S1| [s^{-1}]')
set(ax5,'ydir','reverse')
%xlim([-20 20])
set(gca,'xscale','log','xtick',10.^(-10:10))
axpos=ax5.Position;
ax5.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax6=subplot(nrow,ncol,7); % S2
plot_prof_fast(DATA,filenames_all,filenames,col,select_profile,'fast_S2',1000,true);
% plot_prof_fast(DATA,filenames_all,filenames,col,select_profile,'fast_S2',1);
set(ax6,'ydir','reverse')
xlabel('|S2| [s^{-1}]')
%xlim([-20 20])
set(gca,'xscale','log','xtick',10.^(-10:10))
axpos=ax6.Position;
ax6.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax7=subplot(nrow,ncol,4); % Ax
plot_prof_fast(DATA,filenames_all,filenames,col,select_profile,'A_x',1);
set(ax7,'ydir','reverse')
xlabel('Ax [m s^{-2}]')
axpos=ax7.Position;
ax7.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax8=subplot(nrow,ncol,8); % Ay
plot_prof_fast(DATA,filenames_all,filenames,col,select_profile,'A_y',1);
set(ax8,'ydir','reverse')
xlabel('Ay [m s^{-2}]')
axpos=ax8.Position;
ax8.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];


linkaxes([ax1,ax2,ax3,ax4,ax5,ax6],'y')

%% 3. Shear
% fig_fastS=figure('Units','centimeters','Position',[1 1 9 8]);
% 
% dy=0; % To reduce height of plots
% nrow=1; ncol=2;
% 
% ax1=subplot(nrow,ncol,1); % S1
% hp=plot_prof_fast(DATA,filenames_all,filenames,col,select_profile,'fast_S1',1000,true);
% % hp=plot_prof_fast(DATA,filenames_all,filenames,col,select_profile,'fast_S1',1);
% ylabel('Depth [m]')
% xlabel('|S1| [s^{-1}]')
% set(ax1,'ydir','reverse')
% %xlim([-20 20])
% hleg=legend(hp,leg);
% title(hleg,DATA(1).FAST{1}.date)
% set(gca,'xscale','log','xtick',10.^(-10:10))
% axpos=ax1.Position;
% ax1.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];
% 
% ax2=subplot(nrow,ncol,2); % S2
% plot_prof_fast(DATA,filenames_all,filenames,col,select_profile,'fast_S2',1000,true);
% % plot_prof_fast(DATA,filenames_all,filenames,col,select_profile,'fast_S2',1);
% set(ax2,'ydir','reverse')
% xlabel('|S2| [s^{-1}]')
% %xlim([-20 20])
% set(gca,'xscale','log','xtick',10.^(-10:10))
% axpos=ax2.Position;
% ax2.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];
% 
% linkaxes([ax1,ax2],'y')



%% 2. FP07 data + shear epsilon
fig_turb_all=figure('Units','centimeters','Position',[1 1 18 18]);
dy=0; % To reduce height of plots
nrow=2; ncol=5;

flag_T1='flag_T1';
flag_T2='flag_T2';
% flag_T1='';
% flag_T2='';

ax1=subplot(nrow,ncol,1); % Eps-T1
hp=plot_prof_bin(DATA,filenames_all,filenames,col,select_profile,'eps_T1',1,flag_T1);
ylabel('Depth [m]')
xlabel('\epsilon_{T1} [m^2 s^{-3}]')
set(ax1,'ydir','reverse','xscale','log','xtick',10.^(-14:2:0),'XTickLabelRotation',45)
hleg=legend(hp,leg);
title(hleg,DATA(1).FAST{1}.date)
axpos=ax1.Position;
ax1.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax5=subplot(nrow,ncol,ncol+1); % Eps-T2
plot_prof_bin(DATA,filenames_all,filenames,col,select_profile,'eps_T2',1,flag_T2);
xlabel('\epsilon_{T2} [m^2 s^{-3}]')
set(ax5,'ydir','reverse','xscale','log','xtick',10.^(-14:2:0),'XTickLabelRotation',45)
xlim(ax1.XLim)
axpos=ax5.Position;
ax5.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax9=subplot(nrow,ncol,2); % Eps-S1
hp=plot_prof_ql(DATA,filenames_all,filenames,col,select_profile,'e',1,1);
xlabel('\epsilon_{S1} [m^2 s^{-3}]')
set(ax9,'ydir','reverse','xscale','log','xtick',10.^(-14:2:0),'XTickLabelRotation',45)
axpos=ax9.Position;
ax9.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];
xlim(ax1.XLim)

ax10=subplot(nrow,ncol,ncol+2); % Eps-S2
plot_prof_ql(DATA,filenames_all,filenames,col,select_profile,'e',2,1);
xlabel('\epsilon_{S2} [m^2 s^{-3}]')
set(ax10,'ydir','reverse','xscale','log','xtick',10.^(-14:2:0),'XTickLabelRotation',45)
xlim(ax1.XLim)
axpos=ax10.Position;
ax10.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax2=subplot(nrow,ncol,3); % X-T1
plot_prof_bin(DATA,filenames_all,filenames,col,select_profile,'Xi_T1',1,flag_T1);
xlabel('\chi_{T1} [K^2 s^{-1}]')
set(ax2,'ydir','reverse','xscale','log','xtick',10.^(-14:2:0),'XTickLabelRotation',45)
axpos=ax2.Position;
ax2.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax6=subplot(nrow,ncol,ncol+3); % X-T2
plot_prof_bin(DATA,filenames_all,filenames,col,select_profile,'Xi_T2',1,flag_T2);
xlabel('\chi_{T2} [K^2 s^{-1}]')
set(ax6,'ydir','reverse','xscale','log','xtick',10.^(-14:2:0),'XTickLabelRotation',45)
axpos=ax6.Position;
xlim(ax2.XLim)
ax6.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax3=subplot(nrow,ncol,4); % Kz-T1
plot_prof_bin(DATA,filenames_all,filenames,col,select_profile,'KOsbornCox_T1',1,flag_T1);
xlabel('K_{z,T1} [m^2 s^{-1}]')
set(ax3,'ydir','reverse','xscale','log','xtick',10.^(-14:2:0),'XTickLabelRotation',45)
axpos=ax3.Position;
ax3.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax7=subplot(nrow,ncol,ncol+4); % Kz-T2
plot_prof_bin(DATA,filenames_all,filenames,col,select_profile,'KOsbornCox_T2',1,flag_T2);
xlabel('K_{z,T2} [m^2 s^{-1}]')
set(ax7,'ydir','reverse','xscale','log','xtick',10.^(-14:2:0),'XTickLabelRotation',45)
axpos=ax7.Position;
xlim(ax3.XLim)
ax7.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax4=subplot(nrow,ncol,5); % Kz-T1
plot_prof_bin(DATA,filenames_all,filenames,col,select_profile,'LTuT1',1,flag_T1);
xlabel('L_{T,T1} [m]]')
set(ax4,'ydir','reverse','xscale','log','xtick',10.^(-14:2:0),'XTickLabelRotation',45)
axpos=ax4.Position;
ax4.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax8=subplot(nrow,ncol,ncol+5); % Kz-T2
plot_prof_bin(DATA,filenames_all,filenames,col,select_profile,'LTuT2',1,flag_T2);
xlabel('L_{T,T2} [m]')
set(ax8,'ydir','reverse','xscale','log','xtick',10.^(-14:2:0),'XTickLabelRotation',45)
axpos=ax8.Position;
xlim(ax4.XLim)
ax8.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10],'y')
%% 3. Shear turbulence

fig_turbS=figure('Units','centimeters','Position',[1 1 9 9]);

dy=0.07; % To reduce height of plots

ax1=subplot(1,2,1); % Eps-S1
hp=plot_prof_ql(DATA,filenames_all,filenames,col,select_profile,'e',1,1);
ylabel('Depth [m]')
xlabel('\epsilon_{S1} [m^2 s^{-3}]')
set(ax1,'ydir','reverse','xscale','log','xtick',10.^(-14:2:0),'XTickLabelRotation',45)
hleg=legend(hp,leg);
title(hleg,DATA(1).FAST{1}.date)
axpos=ax1.Position;
ax1.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

ax2=subplot(1,2,2); % Eps-S2
plot_prof_ql(DATA,filenames_all,filenames,col,select_profile,'e',2,1);
xlabel('\epsilon_{S2} [m^2 s^{-3}]')
set(ax2,'ydir','reverse','xscale','log','xtick',10.^(-14:2:0),'XTickLabelRotation',45)
xlim(ax1.XLim)
axpos=ax2.Position;
ax2.Position=[axpos(1),axpos(2)+dy,axpos(3),axpos(4)-dy];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hp]=plot_prof_slow(DATA,filenames_all,filenames,col,select_profile,varname,smoothing_window,absolute)
% Plot all the profiles of a slow variable

hp=[];

if nargin<8
    absolute=false;
end


indprof=0;
for kf=1:length(filenames)
    ind_file=find(strcmp(filenames_all,filenames{kf}),1);
    prof_data=DATA(ind_file).SLOW;
    if strcmp(select_profile,'all')
        for kprof=1:length(prof_data)
            indprof=indprof+1;

            if absolute
                data_plot=movmean(abs(prof_data{kprof}.(varname)),smoothing_window);
            else
                data_plot=movmean(prof_data{kprof}.(varname),smoothing_window);
            end
            hp(end+1)=plot(data_plot,prof_data{kprof}.depth,'Color',col(indprof,:),'LineWidth',1);
            hold on
        end
    else % Plot the longest profile only
        depth_range=NaN(1,length(prof_data)); % [m]
        for kprof=1:length(prof_data)
            depth_range(kprof)=max(prof_data{kprof}.depth)-min(prof_data{kprof}.depth);
        end
        [~,indmax]=max(depth_range);
        if absolute
            data_plot=movmean(abs(prof_data{indmax}.(varname)),smoothing_window);
        else
            data_plot=movmean(prof_data{indmax}.(varname),smoothing_window);
        end
        hp(end+1)=plot(data_plot,prof_data{indmax}.depth,'Color',col(kf,:),'LineWidth',1);
        hold on
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hp]=plot_prof_fast(DATA,filenames_all,filenames,col,select_profile,varname,smoothing_window,absolute)
% Plot all the profiles of a fast variable

hp=[];

if nargin<8
    absolute=false;
end

indprof=0;
for kf=1:length(filenames)
    ind_file=find(strcmp(filenames_all,filenames{kf}),1);
    prof_data=DATA(ind_file).FAST;
    if strcmp(select_profile,'all')
        for kprof=1:length(prof_data)
            indprof=indprof+1;
            if absolute
                %scatter(abs(prof_data{kprof}.(varname)),prof_data{kprof}.depth,'o','SizeData',3,'MarkerFaceColor',col(indprof,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
                %hold on
                hp(end+1)=plot(movmean(abs(prof_data{kprof}.(varname)),smoothing_window),prof_data{kprof}.depth,'Color',col(indprof,:),'LineWidth',1);
                hold on
            else
                hp(end+1)=plot(movmean(prof_data{kprof}.(varname),smoothing_window),prof_data{kprof}.depth,'Color',col(indprof,:),'LineWidth',1);
                hold on
            end
        end
    else  % Plot the longest profile only
        depth_range=NaN(1,length(prof_data)); % [m]
        for kprof=1:length(prof_data)
            depth_range(kprof)=max(prof_data{kprof}.depth)-min(prof_data{kprof}.depth);
        end
        [~,indmax]=max(depth_range);
        if absolute
            %scatter(abs(prof_data{indmax}.(varname)),prof_data{indmax}.depth,'o','SizeData',3,'MarkerFaceColor',col(kf,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.1);
            %hold on
            hp(end+1)=plot(movmean(abs(prof_data{indmax}.(varname)),smoothing_window),prof_data{indmax}.depth,'Color',col(kf,:),'LineWidth',1);
            hold on
        else
            hp(end+1)=plot(movmean(prof_data{indmax}.(varname),smoothing_window),prof_data{indmax}.depth,'Color',col(kf,:),'LineWidth',1);
            hold on
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hp]=plot_prof_bin(DATA,filenames_all,filenames,col,select_profile,varname,smoothing_window,flagname,absolute)
% Plot all the profiles of a binned variable

hp=[];

if nargin<8
    flagname='';
end

if nargin<9
    absolute=false;
end

indprof=0;
for kf=1:length(filenames)
    ind_file=find(strcmp(filenames_all,filenames{kf}),1);
    prof_data=DATA(ind_file).BINNED;
    if strcmp(select_profile,'all')
        for kprof=1:length(prof_data)
            % datavar=prof_data{kprof}.(varname);
            if ~isempty(flagname)
                %datavar(prof_data{kprof}.(flagname)==0)=NaN;
                datavar=prof_data{kprof}.(varname)(prof_data{kprof}.(flagname)==1);
                datadepth=prof_data{kprof}.depth(prof_data{kprof}.(flagname)==1);
            else
                datavar=prof_data{kprof}.(varname);
                datadepth=prof_data{kprof}.depth;
            end
            indprof=indprof+1;
            if absolute
                %scatter(abs(prof_data{kprof}.(varname)),prof_data{kprof}.depth,'o','SizeData',3,'MarkerFaceColor',col(indprof,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
                %hold on
                hp(end+1)=plot(movmean(abs(datavar),smoothing_window),datadepth,'-o','Color',col(indprof,:),'MarkerSize',2,'MarkerFaceColor',col(indprof,:),'LineWidth',1);
                hold on
            else
                hp(end+1)=plot(movmean(datavar,smoothing_window,"omitnan"),datadepth,'-o','Color',col(indprof,:),'MarkerSize',2,'MarkerFaceColor',col(indprof,:),'LineWidth',1);
                hold on
            end
        end
    else  % Plot the longest profile only
        depth_range=NaN(1,length(prof_data)); % [m]
        for kprof=1:length(prof_data)
            depth_range(kprof)=max(prof_data{kprof}.depth)-min(prof_data{kprof}.depth);
        end
        [~,indmax]=max(depth_range);
        % datavar=prof_data{indmax}.(varname);
        if ~isempty(flagname)
            % datavar(prof_data{indmax}.(flagname)==0)=NaN;
            datavar=prof_data{indmax}.(varname)(prof_data{indmax}.(flagname)==1);
            datadepth=prof_data{indmax}.depth(prof_data{indmax}.(flagname)==1);
        else
            datavar=prof_data{indmax}.(varname);
            datadepth=prof_data{indmax}.depth;
        end
        if absolute
            %scatter(abs(prof_data{indmax}.(varname)),prof_data{indmax}.depth,'o','SizeData',3,'MarkerFaceColor',col(kf,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.1);
            %hold on
            hp(end+1)=plot(movmean(abs(datavar),smoothing_window),datadepth,'-o','Color',col(kf,:),'MarkerSize',2,'MarkerFaceColor',col(kf,:),'LineWidth',1);
            hold on
        else
            hp(end+1)=plot(movmean(datavar,smoothing_window,"omitnan"),datadepth,'-o','Color',col(kf,:),'MarkerSize',2,'MarkerFaceColor',col(kf,:),'LineWidth',1);
            hold on
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hp]=plot_prof_ql(DATA,filenames_all,filenames,col,select_profile,varname,indprobe,smoothing_window,absolute)
% Plot all the profiles of a binned variable

hp=[];

if nargin<9
    absolute=false;
end

indprof=0;
for kf=1:length(filenames)
    ind_file=find(strcmp(filenames_all,filenames{kf}),1);
    prof_data=DATA(ind_file).DISS_QL;
    if strcmp(select_profile,'all')
        for kprof=1:length(prof_data)
            indprof=indprof+1;
            if absolute
                %scatter(abs(prof_data{kprof}.(varname)),prof_data{kprof}.depth,'o','SizeData',3,'MarkerFaceColor',col(indprof,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
                %hold on
                hp(end+1)=plot(movmean(abs(prof_data{kprof}.(varname)(indprobe,:)),smoothing_window),prof_data{kprof}.P,'Color',col(indprof,:),'LineWidth',1);
                hold on
            else
                hp(end+1)=plot(movmean(prof_data{kprof}.(varname)(indprobe,:),smoothing_window),prof_data{kprof}.P,'Color',col(indprof,:),'LineWidth',1);
                hold on
            end
        end
    else  % Plot the longest profile only
        depth_range=NaN(1,length(prof_data)); % [m]
        for kprof=1:length(prof_data)
            depth_range(kprof)=max(prof_data{kprof}.P)-min(prof_data{kprof}.P);
        end
        [~,indmax]=max(depth_range);
        if absolute
            %scatter(abs(prof_data{indmax}.(varname)),prof_data{indmax}.depth,'o','SizeData',3,'MarkerFaceColor',col(kf,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.1);
            %hold on
            hp(end+1)=plot(movmean(abs(prof_data{indmax}.(varname)(indprobe,:)),smoothing_window),prof_data{indmax}.P,'Color',col(kf,:),'LineWidth',1);
            hold on
        else
            hp(end+1)=plot(movmean(prof_data{indmax}.(varname)(indprobe,:),smoothing_window),prof_data{indmax}.P,'Color',col(kf,:),'LineWidth',1);
            hold on
        end
    end
end
end