function [data_prof,press_atm]=correct_pressure(data_prof,param,ind_prof_slow,ind_prof_fast,kprof,make_plot,folder_out)
%CORRECT_PRESSURE Compute pressure with respect to the atmospheric
%pressure.
%
%   INPUTS:
%   data_prof (structure): profiling data in a given data file.
%   param (structure): parameters specific to the campaign
%   ind_prof_slow (2*Nprf int array): starting and ending indices of the
%   Nprf different profiles in the data file, computed from the low
%   frequency pressure data.
%   ind_prof_fast (2*Nprf int array): starting and ending indices of the
%   Nprf different profiles in the data file, computed from the high
%   frequency pressure data.
%   kprof (int): profile index between 1 and Nprf.
%   make_plot (boolean) [optional, default = false]: = true to plot data near the interface.
%   folder_out (str) [optional, default = '']: name of the foler where the
%   figure is saved.
%
%   OUTPUTS:
%   data_prof (structure): updated profiling data with corrected pressure.
%   press_atm (double): atmospheric pressure subtracted to the raw pressure
%   to get the pressure relative to the air [dbar].
%
% T. Doda based on S. Piccolroaz, last version: 12.01.2026
%% Compute pressure relative to the atmosphere

% Computed variables:
%   - ind_prof_slow_x (n1*1 int array): indices of the data points near the
%   interface (surface if upward, bottom if downward), based on the low
%   frequency pressure data.
%   - ind_prof_fast_x (n2*1 int array): indices of the data points near the
%   interface (surface if upward, bottom if downward), based on the high
%   frequency pressure data.
%   - ind_prof_slow_plot (n3*1 int array): indices of the data points near the
%   interface for plotting, based on the low
%   frequency pressure data.
%   - ind_prof_fast_plot (n4*1 int array): indices of the data points near the
%   interface for plotting, based on the high
%   frequency pressure data.
%   - istop (n2*1 boolean array): = true where a peak in fast cond or fast
%   temperature is detected (i.e., surface); is empty for downward
%   profiles.
g=9.81; % gravitational acceleration [m/s2]
if nargin<7
    folder_out='';
end
if nargin<6
    make_plot=false;
end
if strcmp(param.info.prof_dir,'up')
    % Correction for precise upward profiles:
    % Identify when the sensor exits the lake looking at the
    % discontinuity of the fast conductimeter or, when not available,
    % of the fast thermistor. Use this information to correct the
    % pressure. In this way, the pressure will be relative to the depth
    % of the fast sensors.
    irange=0.4; % look for the discontinuity within a range of 0.4 seconds around the current start of the profile
    ind_prof_slow_x=ind_prof_slow(2,kprof)-round(irange*data_prof.fs_slow):ind_prof_slow(2,kprof)+round(irange*data_prof.fs_slow);
    ind_prof_fast_x=ind_prof_fast(2,kprof)-round(irange*data_prof.fs_fast):ind_prof_fast(2,kprof)+round(irange*data_prof.fs_fast);

    if strcmp(param.atm_press_method,'cond') && isfield(data_prof,'C1_fast') % Use the microC data
        fprintf(">>> microcond used to detect the surface")
        istop=ischange(data_prof.gradC1(ind_prof_fast_x),'linear'); % Abrupt changes in the slope of gradC1 (=peak)
        istop=find(istop==1,1,'first'); istop=istop-2-ceil(abs(0.001/data_prof.W_fast(ind_prof_fast_x(istop))*data_prof.fs_fast));   % shear probes are 1 mm apart from micro cond. Account also for time response (2 counts to be safe)
    elseif strcmp(param.atm_press_method,'FP07')% Use the FP07 data
        fprintf(">>> FP07 used to detect the surface")
        istop=ischange(data_prof.gradT1(ind_prof_fast_x),'linear','maxnumchanges',1); % Abrupt changes in the slope of gradT1 (=peak)
        istop=find(istop==1,1,'first'); istop=istop-5-ceil(abs(0.003/data_prof.W_fast(ind_prof_fast_x(istop))*data_prof.fs_fast));   % shear probes are 3 mm apart from micro temp. Account also for the time response: 0.007 s *512 count/s ~ 4 counts (5 counts to be safe)
    elseif strcmp(param.atm_press_method,'min')
        press_atm=min(data_prof.P_slow_raw); % Miminum pressure of the file = atmospheric pressure
        istop=[];
    elseif strcmp(param.atm_press_method,'offset')
        istop=[]; 
    else
        error("Wrong atm pressure method!")
    end

    % Compute press_atm with 'cond', 'FP07' and 'offset' methods:
    if ~strcmp(param.atm_press_method,'min')
        if isempty(istop)
            if ~strcmp(param.atm_press_method,'offset') % Display warning only if offset was not the chosen method
                warning("No surface detection! Assume press_atm=0.")
            end
            press_atm=0; % Assumes that pressure was already corrected with coef0?
        else
            press_atm = data_prof.P_fast_raw(ind_prof_fast_x(istop));   % correction
        end
    end

    % Correct pressure only for the specific profile
    P_slow_corr = data_prof.P_slow_raw-press_atm;
    P_fast_corr = data_prof.P_fast_raw-press_atm;
    ind_prof_fast_plot=ind_prof_fast_x;
    ind_prof_slow_plot=ind_prof_slow_x*data_prof.fs_fast/data_prof.fs_slow;
else % Downward profile

    if strcmp(param.atm_press_method,'min')
        press_atm=min(data_prof.P_slow_raw); % Miminum pressure of the file = atmospheric pressure
        istop=[];
    elseif strcmp(param.atm_press_method,'offset')
        istop=[]; 
        press_atm=0;
    else % Methods 'cond' and 'FP07' do not apply to downward profiles
        error("Wrong atm pressure method!") 
    end

    P_slow_corr = data_prof.P_slow_raw-press_atm;
    P_fast_corr = data_prof.P_fast_raw-press_atm;
    irange=0.0; % [s], keep the first profiling index
    ind_prof_slow_x=ind_prof_slow(1,kprof)-round(irange*data_prof.fs_slow):ind_prof_slow(1,kprof)+round(irange*data_prof.fs_slow);
    ind_prof_fast_x=ind_prof_fast(1,kprof)-round(irange*data_prof.fs_fast):ind_prof_fast(1,kprof)+round(irange*data_prof.fs_fast);
    ind_prof_fast_plot=ind_prof_fast_x;
    ind_prof_slow_plot=ind_prof_slow_x*data_prof.fs_fast/data_prof.fs_slow;
    istop=[];
end

fprintf('>>>>> Atmospheric pressure: %0.3f dbar\n',press_atm)
%% Re-load the profile with pmin=0 (T. Doda: this part was originally in resolve_profile_all)
data_prof.ind_prof_slow = get_profile(P_slow_corr,data_prof.W_slow,0,...
    param.info.minvel_detect,param.info.prof_dir,param.info.mindur_detect,data_prof.fs_slow);
data_prof.ind_prof_fast  = get_profile(P_fast_corr,data_prof.W_fast,0,...
    param.info.minvel_detect,param.info.prof_dir,param.info.mindur_detect,data_prof.fs_fast);
[~,data_prof.Nprf] = size(data_prof.ind_prof_slow); % Number of profiles

if size(data_prof.ind_prof_slow,2)~=data_prof.Nprf
    error('Change in number of profiles!')
end

%% Save the corrected pressure only for the specific profile
data_prof.P_slow=NaN(length(data_prof.P_slow),1);
data_prof.P_slow(data_prof.ind_prof_slow(1,kprof):data_prof.ind_prof_slow(2,kprof))=P_slow_corr(data_prof.ind_prof_slow(1,kprof):data_prof.ind_prof_slow(2,kprof));
data_prof.P_fast=NaN(length(data_prof.P_fast),1);
data_prof.P_fast(data_prof.ind_prof_fast(1,kprof):data_prof.ind_prof_fast(2,kprof))=P_fast_corr(data_prof.ind_prof_fast(1,kprof):data_prof.ind_prof_fast(2,kprof));

%% Compute the depth
mrho = cumsum(data_prof.rhoTS)./(1:length(data_prof.rhoTS))'; % [kg/m3]
data_prof.depth=10000*data_prof.P_slow./(mrho*g); % [m]

%% Plot
if make_plot

    if length(ind_prof_fast_plot)>1

        % Plot from S. Piccolroaz
        figure
        set(gcf, 'Units', 'centimeters', 'Position', [1 1 20 8]);
        plot(ind_prof_fast_plot,data_prof.P_fast(ind_prof_fast_plot),'.-k'); hold on
        if isfield(data_prof,'C1_fast')
            plot(ind_prof_fast_plot,normalize(data_prof.C1_fast(ind_prof_fast_plot))/2,'-g');
        end
        plot(ind_prof_fast_plot,normalize(data_prof.T1_fast(ind_prof_fast_plot))/2,'.-r');
        plot(ind_prof_fast_plot,data_prof.W_fast(ind_prof_fast_plot),'.-k');
        plot(ind_prof_fast_plot,data_prof.sh1(ind_prof_fast_plot)/50,'.-m'); hold on
        plot(ind_prof_slow_plot,normalize(data_prof.(CTD_C)(ind_prof_slow_x))/2,'.-b'); hold on
        plot([ind_prof_fast_plot(1) ind_prof_fast_plot(end)],[0 0],'--k');
        if ~isempty(istop)
            plot([ind_prof_fast_plot(istop) ind_prof_fast_plot(istop)],[-1 1],'--g');
        end
        ylim([-1 1])
    end

    % Plot to compare before-after pressure correction
    indslow=data_prof.ind_prof_slow(1,kprof):data_prof.ind_prof_slow(2,kprof);
    indfast=data_prof.ind_prof_fast(1,kprof):data_prof.ind_prof_fast(2,kprof);

    figure

    ax1=subplot(1,2,1);
    plot(data_prof.(param.CTD_T)(indslow),data_prof.P_slow_raw(indslow))
    hold on
    plot(data_prof.(param.CTD_T)(indslow),data_prof.P_slow(indslow))
    plot(data_prof.(param.CTD_T)(indslow),data_prof.depth(indslow),'k')
    xlabel('Temperature [°C]')
    ylabel('Pressure [dbar]/Depth [m]')
    legend('Praw','Pcorr','Depth')
    set(gca,'ydir','reverse')
    title(gca,'Slow data')

    ax2=subplot(1,2,2);
    plot(data_prof.T1_fast(indfast),data_prof.P_fast_raw(indfast))
    hold on
    plot(data_prof.T1_fast(indfast),data_prof.P_fast(indfast))
    xlabel('Temperature [°C]')
    legend('Praw','Pcorr')
    set(gca,'ydir','reverse')
    title(gca,'Fast data')

    linkaxes([ax1,ax2],'y')

    if ~exist([folder_out,'Figures'],'dir')
        mkdir([folder_out,'Figures'])
    end
    saveas(gcf,[folder_out,'Figures\Pcorrection',num2str(kprof),'.fig'])
    exportgraphics(gcf,[folder_out,'Figures\Pcorrection',num2str(kprof),'.png'],'Resolution',400)

end


end