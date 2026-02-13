function [] = check_data(param,data_prof,savedata,savefig,simple_prof,change_axis_top)
%CHECK_DATA Plot the profile from a .P file.
%
%   INPUTS:
%   param (structure): parameters produced by load_parameters_lakename
%   function.
%   data_prof (structure): data from a .P file imported with odas_p2mat.
%   savedata (boolean) [optional]: if =true, save the data. Default:
%   true.
%   savefig (boolean) [optional]: if =true, save the figures. Default:
%   true.
%   simple_prof (boolean) [optional]: if =true, simple version of the plot
%   (five panels). Default: false.
%   change_axis_top (boolean) [optional]: if =true, change the axis order for interactivity. Default: true.
%
%   OUTPUTS:
%   None
%
% T.Doda, 30.01.2026


%% Parameters 

if nargin<3
    savedata=true;
end

if nargin<4
    savefig=true;
end

if nargin<5
    simple_prof=false;
end

if nargin<6
    change_axis_top=true;
end

axlim_shear=[-10,10];
axlim_dT=[-3,3];

%% Analyze each data file

if savedata
    data_folder=[param.folder,'..\Check_data'];
    if ~exist(data_folder, 'dir')
        mkdir(data_folder)
    end

    
end

output_data_file=[data_folder,'\data_',data_prof.filename,'.mat'];
    
if exist(output_data_file,'file')
    analyze_raw=input(sprintf('Mat file already exists for %s: load it (1) or re-analyze the data (2)?',data_prof.filename));
else
    analyze_raw=2;
end


if analyze_raw==2 % Use the structure data_prof
    ind_prof_slow=data_prof.ind_prof_slow_initial;
    ind_prof_fast=data_prof.ind_prof_fast_initial;

    [~,Nprf] = size(ind_prof_slow); % Number of profiles
    if Nprf==0
        error('>>> No profile for %s: data not saved',data_prof.filename)
    end

    data_prof.depth_slow=data_prof.P_slow-min(data_prof.P_slow);
    data_prof.depth_fast=data_prof.P_fast-min(data_prof.P_fast);

    data_prof.ind_prof_slow=ind_prof_slow;
    data_prof.ind_prof_fast=ind_prof_fast;
    data_prof.Nprf=Nprf; 

    varnames_slow={'t_slow','tnum_slow','tdate_slow','depth_slow','P_slow','W_slow',param.CTD_T,param.CTD_C,param.CTD_Chl,param.CTD_Turb,'Salin','Salin_smooth','rho','T1_slow','T2_slow'};
    varnames_fast={'t_fast','tnum_fast','tdate_fast','depth_fast','P_fast','W_fast',...
        'T1_fast','T2_fast','gradT1','gradT2','sh1','sh2','C1_fast','C2_fast',...
        'gradC1','gradC2','Ax','Ay','Az'};
    prof_slow=struct();
    prof_fast=struct();
    
    for kvar=1:length(varnames_slow)
        if isfield(data_prof,varnames_slow{kvar})
            for kprof=1:data_prof.Nprf
                ind_slow=data_prof.ind_prof_slow(1,kprof):data_prof.ind_prof_slow(2,kprof);
                prof_slow(kprof).(varnames_slow{kvar})=data_prof.(varnames_slow{kvar})(ind_slow);
            end
        end
    end

    for kvar=1:length(varnames_fast)
        if isfield(data_prof,varnames_fast{kvar})
            for kprof=1:data_prof.Nprf
                ind_fast=data_prof.ind_prof_fast(1,kprof):data_prof.ind_prof_fast(2,kprof);
                prof_fast(kprof).(varnames_fast{kvar})=data_prof.(varnames_fast{kvar})(ind_fast);
            end
        end
    end

    data_prof.prof_slow=prof_slow;
    data_prof.prof_fast=prof_fast;

    % 3. Save the data in a structure

    DATA=data_prof;
    if savedata
        save(output_data_file,'DATA')
        disp('Data saved')
    end
else
    load(output_data_file,'DATA');
end



%% Plot the profile data with 10 subplots
fig_folder=[param.folder,'..\Check_figures'];
if savefig
    if ~exist(fig_folder, 'dir')
        %rmdir([param.folder,'Check_figures'],'s')
        mkdir(fig_folder)
    end
    
end

if ~simple_prof
    % 1. Plot of pressure & velocity
    figure('Units','centimeters','Position',[1 1 20 10])

    ax1=subplot(2,1,1);
    plot(DATA.tdate_slow,DATA.P_slow,'k')
    hold on
    for kprof=1:DATA.Nprf
        plot(DATA.prof_slow(kprof).tdate_slow,DATA.prof_slow(kprof).P_slow,'r','linewidth',1)
    end
    ylabel('Pressure [dbar]')
    title(sprintf('File %s',erase(DATA.filename,'_')))

    ax2=subplot(2,1,2);
    Psmooth=movmean(DATA.P_slow,5);
    sink_speed1=vertcat(NaN,diff(Psmooth)./diff(DATA.t_slow)); % [dbar/s]
    %sink_speed2=vertcat(NaN,NaN,(DATA.P_slow(5:end)-DATA.P_slow(1:end-4))./(DATA.t_slow(5:end)-DATA.t_slow(1:end-4)),NaN,NaN);
    plot(DATA.tdate_slow,DATA.W_slow)
    hold on
    plot(DATA.tdate_slow,sink_speed1)
    %plot(DATA.tdate_slow,sink_speed2)
    ylabel('Speed [dbar/s]')

    linkaxes([ax1,ax2],'x')

    if savefig
        saveas(gcf,[param.folder,'..\Check_figures\pressure_',erase(DATA.filename,'_'),'.fig'])
        exportgraphics(gcf,[param.folder,'..\Check_figures\pressure_',erase(DATA.filename,'_'),'_P',num2str(kprof),'.png'],'Resolution',400)
    end

    % 2. Profiles
    colmat=lines;
    for kprof=1:DATA.Nprf

        figure('Units','centimeters','Position',[1 1 20 20])
        data_slow=DATA.prof_slow(kprof);
        data_fast=DATA.prof_fast(kprof);

        % Depth & velocity
        axtop=subplot(2,5,1); 
        % Dsmooth_prof=movmean(DATA.depth_slow(ind_plot),5);
        %sink_speed_prof=vertcat(NaN,diff(Dsmooth_prof)./diff(DATA.t_slow(ind_plot))); % [dbar/s]
        sink_speed_prof=vertcat(NaN,diff(data_slow.depth_slow)./diff(data_slow.t_slow)); % [dbar/s]
        plot(sink_speed_prof,data_slow.depth_slow,'k');
        median_speed=median(sink_speed_prof,'omitnan');
        hold on
        plot([median_speed,median_speed],[data_slow.depth_slow(1),data_slow.depth_slow(end)],'--r');
        text(median_speed-0.01,data_slow.depth_slow(1)+0.5*(data_slow.depth_slow(end)-data_slow.depth_slow(1)), ...
            sprintf('%0.2f m/s',median_speed),'color','r','HorizontalAlignment','right')
        xlabel('Speed [m/s]')
        ylabel('Depth [m]') 
        set(axtop,'ydir','reverse')

        % Temp & Cond
        axprof1=subplot(2,5,2);
        plot(data_slow.(param.CTD_T),data_slow.depth_slow,'k')
        xlabel('Temp [°C]') 
        axprof1_top = add_xaxis(axprof1,'on');
        plot(data_slow.(param.CTD_C),data_slow.depth_slow,'color',[0 0 0.8])
        %plot(data_slow.Salin_smooth,data_slow.depth_slow,'color',[0 0 0.8])
        set(axprof1_top,'xcolor',[0 0 0.8])
        xlabel('Cond [mS/cm]','color',[0 0 0.8])
        %xlabel('Salin [PSU]','color',[0 0 0.8])
        set(axprof1,'ydir','reverse')
        set(axprof1_top,'ydir','reverse')
        

        % Ch-a & Turb
        axprof2=subplot(2,5,3);
        plot(data_slow.(param.CTD_Chl),data_slow.depth_slow,'Color',[0 0.8 0])
        xlabel('Chl [ppb]','Color',[0 0.8 0])
        set(axprof2,'xcolor',[0 0.8 0])
        axprof2_top = add_xaxis(axprof2,'on');
        plot(data_slow.(param.CTD_Turb),data_slow.depth_slow,'Color',[0.8 0 0])
        set(axprof2_top,'xcolor',[0.8 0 0])
        xlabel('Turb [FTU]','Color',[0.8 0 0])
        set(axprof2,'ydir','reverse')
        set(axprof2_top,'ydir','reverse')
        title(sprintf('File %s: profile %d/%d',erase(DATA.filename,'_'),kprof,DATA.Nprf))


        % microT
        axprof3=subplot(2,5,4);
        if isfield(data_fast,'T1_fast')
            plot(data_fast.T1_fast,data_fast.depth_fast,'color',colmat(1,:))
        end
        xlabel('T1 [°C]','color',colmat(1,:))
        set(axprof3,'xcolor',colmat(1,:))
        axprof3_top = add_xaxis(axprof3,'on');
        if isfield(data_fast,'T2_fast')
            plot(data_fast.T2_fast,data_fast.depth_fast,'color',colmat(2,:))
        end
        set(axprof3_top,'xcolor',colmat(2,:))
        xlabel('T2 [°C]','color',colmat(2,:))
        [axprof3.XLim,axprof3_top.XLim] = scale_limits(axprof3.XLim,axprof3_top.XLim);
        set(axprof3,'ydir','reverse')
        set(axprof3_top,'ydir','reverse')
        
        % gradT
        axprof4=subplot(2,5,5);
        if isfield(data_fast,'gradT1')
            plot(data_fast.gradT1,data_fast.depth_fast,'color',colmat(1,:))
        end
        xlabel('dT1/dz [°C/m]','color',colmat(1,:))
        set(axprof4,'xcolor',colmat(1,:))
        axprof4_top = add_xaxis(axprof4,'on');
        if isfield(data_fast,'gradT2')
            plot(data_fast.gradT2,data_fast.depth_fast,'color',colmat(2,:))
        end
        set(axprof4_top,'xcolor',colmat(2,:))
        xlabel('dT2/dz [°C/m]','color',colmat(2,:))
        %[axprof4.XLim,axprof4_top.XLim] = scale_limits(axprof4.XLim,axprof4_top.XLim);
        axprof4.XLim=axlim_dT;axprof4_top.XLim=axlim_dT;
        set(axprof4,'ydir','reverse')
        set(axprof4_top,'ydir','reverse')

        % Shear
        axprof5=subplot(2,5,6);
        if isfield(data_fast,'sh1')
            plot(data_fast.sh1,data_fast.depth_fast,'color',colmat(1,:))
        end
        xlabel('S1 [s^{-1}]','color',colmat(1,:))
        ylabel('Depth [m]')
        set(axprof5,'xcolor',colmat(1,:))
        axprof5_top = add_xaxis(axprof5,'on');
        if isfield(data_fast,'sh2')
            plot(data_fast.sh2,data_fast.depth_fast,'color',colmat(2,:))
        end
        set(axprof5_top,'xcolor',colmat(2,:))
        xlabel('S2 [s^{-1}]','color',colmat(2,:))
        %[axprof5.XLim,axprof5_top.XLim] = scale_limits(axprof5.XLim,axprof5_top.XLim);
        axprof5.XLim=axlim_shear;axprof5_top.XLim=axlim_shear;
        set(axprof5,'ydir','reverse')
        set(axprof5_top,'ydir','reverse')

        % microC
        axprof6=subplot(2,5,7);
        if isfield(data_fast,'C1_fast')
            plot(data_fast.C1_fast,data_fast.depth_fast,'color',colmat(1,:))
        end
        xlabel('C1 [mS/cm]','color',colmat(1,:))
        set(axprof6,'xcolor',colmat(1,:))
        axprof6_top = add_xaxis(axprof6,'on');
        if isfield(data_fast,'C2_fast')
            plot(data_fast.C2_fast,data_fast.depth_fast,'color',colmat(2,:))
        end
        set(axprof6_top,'xcolor',colmat(2,:))
        xlabel('C2 [mS/cm]','color',colmat(2,:))
        [axprof6.XLim,axprof6_top.XLim] = scale_limits(axprof6.XLim,axprof6_top.XLim);
        set(axprof6,'ydir','reverse')
        set(axprof6_top,'ydir','reverse')

        % gradC
        axprof7=subplot(2,5,8);
        if isfield(data_fast,'gradC1')
            plot(data_fast.gradC1,data_fast.depth_fast,'color',colmat(1,:))
        end
        xlabel('dC1/dz [mS.cm^{-1}.m{-1}]','color',colmat(1,:))
        set(axprof7,'xcolor',colmat(1,:))
        axprof7_top = add_xaxis(axprof7,'on');
        if isfield(data_fast,'gradC2')
            plot(data_fast.gradC2,data_fast.depth_fast,'color',colmat(2,:))
        end
        set(axprof7_top,'xcolor',colmat(2,:))
        xlabel('dC2/dz [mS.cm^{-1}.m{-1}]','color',colmat(2,:))
        [axprof7.XLim,axprof7_top.XLim] = scale_limits(axprof7.XLim,axprof7_top.XLim);
        set(axprof7,'ydir','reverse')
        set(axprof7_top,'ydir','reverse')

        % Ax, Ay
        axprof8=subplot(2,5,9);
        plot(data_fast.Ax,data_fast.depth_fast,'color',colmat(3,:))
        xlabel('Ax [m.s^{-2}]','color',colmat(3,:))
        set(axprof8,'xcolor',colmat(3,:))
        axprof8_top = add_xaxis(axprof8,'on');
        plot(data_fast.Ay,data_fast.depth_fast,'color',colmat(4,:))
        set(axprof8_top,'xcolor',colmat(4,:))
        xlabel('Ay [m.s^{-2}]','color',colmat(4,:))
        [axprof8.XLim,axprof8_top.XLim] = scale_limits(axprof8.XLim,axprof8_top.XLim);
        set(axprof8,'ydir','reverse')
        set(axprof8_top,'ydir','reverse')

        % Az
        axprof9=subplot(2,5,10);
        if isfield(data_fast,'Az')   
            plot(data_fast.Az,data_fast.depth_fast,'color','k')
            xlabel('Az [m.s^{-2}]')
            set(axprof9,'ydir','reverse')
        end
        linkaxes([axprof1,axprof1_top,axprof2,axprof2_top,axprof3,axprof3_top,...
            axprof4,axprof4_top,axprof5,axprof5_top,axprof6,axprof6_top,...
            axprof7,axprof7_top,axprof8,axprof8_top,axprof9],'y');

        if savefig
            saveas(gcf,[fig_folder,'\profile_',erase(DATA.filename,'_'),'_P',num2str(kprof),'.fig'])
            exportgraphics(gcf,[fig_folder,'\profile_',erase(DATA.filename,'_'),'_P',num2str(kprof),'.png'],'Resolution',400)
            disp('Figure saved')
        end
    end
end
%% Plot simpler version of profiles (5 panels)

if simple_prof
  
                                        colmat=lines;
    for kprof=1:DATA.Nprf

        figure('Units','centimeters','Position',[1 1 30 10])
        data_slow=DATA.prof_slow(kprof);
        data_fast=DATA.prof_fast(kprof);

        % Depth & velocity
        axtop=subplot(1,5,1); 
        % Dsmooth_prof=movmean(DATA.depth_slow(ind_plot),5);
        %sink_speed_prof=vertcat(NaN,diff(Dsmooth_prof)./diff(DATA.t_slow(ind_plot))); % [dbar/s]
        sink_speed_prof=vertcat(NaN,diff(data_slow.depth_slow)./diff(data_slow.t_slow)); % [dbar/s]
        plot(sink_speed_prof,data_slow.depth_slow,'k');
        median_speed=median(sink_speed_prof,'omitnan');
        hold on
        plot([median_speed,median_speed],[data_slow.depth_slow(1),data_slow.depth_slow(end)],'--r');
        text(median_speed-0.01,data_slow.depth_slow(1)+0.5*(data_slow.depth_slow(end)-data_slow.depth_slow(1)), ...
            sprintf('%0.2f m/s',median_speed),'color','r','HorizontalAlignment','right')
        xlabel('Speed [m/s]')
        ylabel('Depth [m]') 
        set(axtop,'ydir','reverse')

        % Temp, Chl-a, Turb
        axprof1=subplot(1,5,2);
        plot(data_slow.(param.CTD_T),data_slow.depth_slow,'k')
        xlabel('Temp [°C]') 
        axprof1_top = add_xaxis(axprof1,'on');
        plot(data_slow.(param.CTD_Chl),data_slow.depth_slow,'Color',[0 0.8 0])
        plot(data_slow.(param.CTD_Turb),data_slow.depth_slow,'Color',[0.8 0 0])
        set(axprof1_top,'xcolor','k')
        xlabel('Turb [FTU]/Chl [ppb]')
        set(axprof1,'ydir','reverse')
        set(axprof1_top,'ydir','reverse')


        % microT
        axprof3=subplot(1,5,3);
        if isfield(data_fast,'T1_fast')
            plot(data_fast.T1_fast,data_fast.depth_fast,'color',colmat(1,:))
        end
        xlabel('T1 [°C]','color',colmat(1,:))
        set(axprof3,'xcolor',colmat(1,:))
        axprof3_top = add_xaxis(axprof3,'on');
        if isfield(data_fast,'T2_fast')
            plot(data_fast.T2_fast,data_fast.depth_fast,'color',colmat(2,:))
        end
        set(axprof3_top,'xcolor',colmat(2,:))
        xlabel('T2 [°C]','color',colmat(2,:))
        [axprof3.XLim,axprof3_top.XLim] = scale_limits(axprof3.XLim,axprof3_top.XLim);
        set(axprof3,'ydir','reverse')
        set(axprof3_top,'ydir','reverse')
        title(sprintf('File %s: profile %d/%d',erase(DATA.filename,'_'),kprof,DATA.Nprf))
        
        % gradT
        axprof4=subplot(1,5,4);
        if isfield(data_fast,'gradT1')
            plot(data_fast.gradT1,data_fast.depth_fast,'color',colmat(1,:))
        end
        xlabel('dT1/dz [°C/m]','color',colmat(1,:))
        set(axprof4,'xcolor',colmat(1,:))
        axprof4_top = add_xaxis(axprof4,'on');
        if isfield(data_fast,'gradT2')
            plot(data_fast.gradT2,data_fast.depth_fast,'color',colmat(2,:))
        end
        set(axprof4_top,'xcolor',colmat(2,:))
        xlabel('dT2/dz [°C/m]','color',colmat(2,:))
        [axprof4.XLim,axprof4_top.XLim] = scale_limits(axprof4.XLim,axprof4_top.XLim);
        set(axprof4,'ydir','reverse')
        set(axprof4_top,'ydir','reverse')

        % Shear
        axprof5=subplot(1,5,5);
        if isfield(data_fast,'sh1')
            plot(data_fast.sh1,data_fast.depth_fast,'color',colmat(1,:))
        end
        xlabel('S1 [s^{-1}]','color',colmat(1,:))
        ylabel('Depth [m]')
        set(axprof5,'xcolor',colmat(1,:))
        axprof5_top = add_xaxis(axprof5,'on');
        if isfield(data_fast,'sh2')
            plot(data_fast.sh2,data_fast.depth_fast,'color',colmat(2,:))
        end
        set(axprof5_top,'xcolor',colmat(2,:))
        xlabel('S2 [s^{-1}]','color',colmat(2,:))
        [axprof5.XLim,axprof5_top.XLim] = scale_limits(axprof5.XLim,axprof5_top.XLim);
        set(axprof5,'ydir','reverse')
        set(axprof5_top,'ydir','reverse')

        
        if savefig
            saveas(gcf,[fig_folder,'\simple_profile_',erase(DATA.filename,'_'),'_P',num2str(kprof),'.fig'])
            exportgraphics(gcf,[fig_folder,'\simple_profile_',erase(DATA.filename,'_'),'_P',num2str(kprof),'.png'],'Resolution',400)
            disp('Figure saved')
        end
    end
    
end

%% Change axis order

if change_axis_top
    uistack(axprof1, 'top')
    set(axprof1,'color','none')
    uistack(axprof2, 'top')
    set(axprof2,'color','none')
    uistack(axprof3, 'top')
    set(axprof3,'color','none')
    uistack(axprof4, 'top')
    set(axprof4,'color','none')
    uistack(axprof5, 'top')
    set(axprof5,'color','none')
    if ~simple_prof
        uistack(axprof6, 'top')
        set(axprof6,'color','none')
        uistack(axprof7, 'top')
        set(axprof7,'color','none')
        uistack(axprof8, 'top')
        set(axprof8,'color','none')
    end
else
    uistack(axprof1_top, 'top')
    uistack(axprof2_top, 'top')
    uistack(axprof3_top, 'top')
    uistack(axprof4_top, 'top')
    uistack(axprof5_top, 'top')
    if ~simple_prof
        uistack(axprof6_top, 'top')
        uistack(axprof7_top, 'top')
        uistack(axprof8_top, 'top')
    end
end



%% Plot the profile data & pressure time series
% for indfile=1:length(param.filename_list)
%     output_data_file=[param.folder,'data_',param.filename_list{indfile},'.mat'];
%     load(output_data_file,'DATA');
% 
%     % 1. Plot of pressure & velocity
%     figure('Units','centimeters','Position',[5 5 20 10])
% 
%     ax1=subplot(2,1,1);
%     plot(DATA.tdate_slow,DATA.P_slow)
%     ylabel('Pressure [dbar]')
% 
%     ax2=subplot(2,1,2);
%     Psmooth=movmean(DATA.P_slow,5);
%     sink_speed1=vertcat(NaN,diff(Psmooth)./diff(DATA.t_slow)); % [dbar/s]
%     %sink_speed2=vertcat(NaN,NaN,(DATA.P_slow(5:end)-DATA.P_slow(1:end-4))./(DATA.t_slow(5:end)-DATA.t_slow(1:end-4)),NaN,NaN);
%     plot(DATA.tdate_slow,DATA.W_slow)
%     hold on
%     plot(DATA.tdate_slow,sink_speed1)
%     %plot(DATA.tdate_slow,sink_speed2)
%     ylabel('Speed [dbar/s]')
% 
%     linkaxes([ax1,ax2],'x')
% 
%     % 2. Profiles
%     colmat=lines;
%     for kprof=1:DATA.Nprf
% 
%         figure('Units','centimeters','Position',[5 5 20 10])
%         data_slow=DATA.prof_slow(kprof);
%         data_fast=DATA.prof_fast(kprof);
% 
%         % Depth & velocity
%         ax1=subplot(3,5,1);
%         vertborder=(1-ax1.Position(2)-ax1.Position(4))*1.1;
%         horborder=0.5*ax1.Position(1);
%         vertbetween=(1-2*vertborder-3*ax1.Position(4))/2;
%         horbetween=(1-2*horborder-5*ax1.Position(3))/4;
%         axtop=axes('Position',[horborder,ax1.Position(2),1-2*horborder,ax1.Position(4)]);
%         set(ax1,'visible','off')
%         plot(axtop,data_slow.tdate_slow,data_slow.depth_slow,'k');
%         ylabel('Depth [m]') 
%         hold on
%         yyaxis right
%         % Dsmooth_prof=movmean(DATA.depth_slow(ind_plot),5);
%         %sink_speed_prof=vertcat(NaN,diff(Dsmooth_prof)./diff(DATA.t_slow(ind_plot))); % [dbar/s]
%         sink_speed_prof=vertcat(NaN,diff(data_slow.depth_slow)./diff(data_slow.t_slow)); % [dbar/s]
%         plot(data_slow.tdate_slow,sink_speed_prof,'m');
%         set(gca,'ycolor','m')
%         ylabel('Speed [m/s]','color','m')
%         title(sprintf('File %s: profile %d/%d',erase(DATA.filename,'_'),kprof,DATA.Nprf))
%         yyaxis left
%         set(axtop,'ydir','reverse')
% 
% 
%         % Temp & Cond
%         heightprof=0.9*(1-2*vertborder-vertbetween-ax1.Position(4));
%         axprof1=axes('Position',[horborder,vertborder,ax1.Position(3),heightprof]);
%         plot(data_slow.JAC_T,data_slow.depth_slow,'k')
%         xlabel('Temp [°C]')
%         ylabel('Depth [m]') 
%         axprof1_top = add_xaxis(axprof1,'on');
%         %plot(data_slow.JAC_C,data_slow.depth_slow,'color',[0 0 0.8])
%         plot(data_slow.JAC_S_smooth,data_slow.depth_slow,'color',[0 0 0.8])
%         set(axprof1_top,'xcolor',[0 0 0.8])
%         %xlabel('Cond [mS/cm]','color',[0 0 0.8])
%         xlabel('Salin [PSU]','color',[0 0 0.8])
%         set(axprof1,'ydir','reverse')
%         set(axprof1_top,'ydir','reverse')
% 
%         % Ch-a & Turb
%         axprof2=axes('Position',[ax1.Position(3)+horborder+horbetween,vertborder,ax1.Position(3),heightprof]);
%         plot(data_slow.Chlorophyll,data_slow.depth_slow,'Color',[0 0.8 0])
%         xlabel('Chl [ppb]','Color',[0 0.8 0])
%         set(axprof2,'xcolor',[0 0.8 0])
%         axprof2_top = add_xaxis(axprof2,'on');
%         plot(data_slow.Turbidity,data_slow.depth_slow,'Color',[0.8 0 0])
%         set(axprof2_top,'xcolor',[0.8 0 0])
%         xlabel('Turb [FTU]','Color',[0.8 0 0])
%         set(axprof2,'ydir','reverse')
%         set(axprof2_top,'ydir','reverse')
% 
%         axprof3=axes('Position',[2*ax1.Position(3)+horborder+2*horbetween,vertborder,ax1.Position(3),heightprof]);
%         if isfield(data_fast,'T1_fast')
%             plot(data_fast.T1_fast,data_fast.depth_fast,'color',colmat(1,:))
%         end
%         xlabel('T1 [°C]','color',colmat(1,:))
%         set(axprof3,'xcolor',colmat(1,:))
%         axprof3_top = add_xaxis(axprof3,'on');
%         if isfield(data_fast,'T2_fast')
%             plot(data_fast.T2_fast,data_fast.depth_fast,'color',colmat(2,:))
%         end
%         set(axprof3_top,'xcolor',colmat(2,:))
%         xlabel('T2 [°C]','color',colmat(2,:))
%         [axprof3.XLim,axprof3_top.XLim] = scale_limits(axprof3.XLim,axprof3_top.XLim);
%         set(axprof3,'ydir','reverse')
%         set(axprof3_top,'ydir','reverse')
% 
%         axprof4=axes('Position',[3*ax1.Position(3)+horborder+3*horbetween,vertborder,ax1.Position(3),heightprof]);
%         if isfield(data_fast,'sh1')
%             plot(data_fast.sh1,data_fast.depth_fast,'color',colmat(1,:))
%         end
%         xlabel('S1 [s^{-1}]','color',colmat(1,:))
%         set(axprof4,'xcolor',colmat(1,:))
%         axprof4_top = add_xaxis(axprof4,'on');
%         if isfield(data_fast,'sh2')
%             plot(data_fast.sh2,data_fast.depth_fast,'color',colmat(2,:))
%         end
%         set(axprof4_top,'xcolor',colmat(2,:))
%         xlabel('S2 [s^{-1}]','color',colmat(2,:))
%         [axprof4.XLim,axprof4_top.XLim] = scale_limits(axprof4.XLim,axprof4_top.XLim);
%         set(axprof4,'ydir','reverse')
%         set(axprof4_top,'ydir','reverse')
% 
%         axprof5=axes('Position',[4*ax1.Position(3)+horborder+4*horbetween,vertborder,ax1.Position(3),heightprof]);
%         if isfield(data_fast,'C1')
%             plot(data_fast.C1,data_fast.depth_fast,'color',colmat(1,:))
%         end
%         xlabel('C1 [mS/cm]','color',colmat(1,:))
%         set(axprof5,'xcolor',colmat(1,:))
%         axprof5_top = add_xaxis(axprof5,'on');
%         if isfield(data_fast,'C2')
%             plot(data_fast.C2,data_fast.depth_fast,'color',colmat(2,:))
%         end
%         set(axprof5_top,'xcolor',colmat(2,:))
%         xlabel('C2 [mS/cm]','color',colmat(2,:))
%         [axprof5.XLim,axprof5_top.XLim] = scale_limits(axprof3.XLim,axprof5_top.XLim);
%         set(axprof5,'ydir','reverse')
%         set(axprof5_top,'ydir','reverse')
% 
%         linkaxes([axprof1,axprof1_top,axprof2,axprof2_top,axprof3,axprof3_top,axprof4,axprof4_top,axprof5,axprof5_top],'y');
% 
% 
%     end
% 
% end
% 
% 

end
