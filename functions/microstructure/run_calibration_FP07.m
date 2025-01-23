function [data_prof,cfgfile_cal] = run_calibration_FP07(param,data_prof,CTD_T,ind_prof,Nprf,modified_data_file,cfgfile_mod,make_plot,folder_out)
%RUN_CALIBRATION_FP07 Calibrate FP07
%   Detailed explanation goes here
% make_plot: optional (=False by default)
% folder_out: optional
%%
if nargin<9
    folder_out='';
end
if nargin<8
    make_plot=false;
end
%% Calibration
fprintf(">>> FP07 calibration\n")
m=[]; % position of values along a profile
for i=1:Nprf
    m = [m,ind_prof(1,i):ind_prof(2,i)];
end

temp_info = cal_FP07_in_situ; % Get the calibration parameters
if max(data_prof.(CTD_T)(m))-min(data_prof.(CTD_T)(m))<=8 % Less than 8°C of difference between the min and max temperature
    temp_info.order=1;
else
    temp_info.order=2;
end





% Run calibration and change parameters in the config file
cfgfile_cal=[cfgfile_mod '_calibrated'];
copyfile([cfgfile_mod '.cfg'],[cfgfile_cal '.cfg']);
if param.config.T1
    % Save precalibration data:
    T1_precal=data_prof.T1_fast;
    [T_01,beta1,lag1] = cal_FP07_in_situ_EPFL(data_prof,m,CTD_T,'T1',cfgfile_cal,temp_info);
end

if param.config.T2
    T2_precal=data_prof.T2_fast;
    [T_02,beta2,lag2] = cal_FP07_in_situ_EPFL(data_prof,m,CTD_T,'T2',cfgfile_cal,temp_info);
end


%% Re-pacth the P-file
% if exist(modified_data_file,'file')
%     delete(modified_data_file)
% end
if ~exist(modified_data_file,'file')
    error('No patched data file available')
end
%copyfile(original_data_file,modified_data_file); % Create a copy of the P file where the configuration file will be modified
patch_setupstr(modified_data_file,cfgfile_cal); % patch the new cfg file

%% Re-convert the data to physical units
clear data_prof
default_parameters=odas_p2mat;
data_prof = odas_p2mat(modified_data_file,default_parameters);          % re-convert data to physical units
if ~strcmp(data_prof.input_parameters.gradT_method,'high_pass')
    error('Error: the gradT_method should be high_pass (if first_difference, pass the info to get_scalar_spectra_odas)')
end

%% Make plot
if make_plot
    figure('Units','centimeters','Position',[1 1 18 10])
    leg={};

    subplot(1,3,1) % Before calibration
    if param.config.T1
        plot(T1_precal,data_prof.P_fast,'LineWidth',1)
        hold on
        leg{end+1}='T1';
    end
    if param.config.T2
        plot(T2_precal,data_prof.P_fast,'LineWidth',1)
        hold on
        leg{end+1}='T2';
    end
    plot(data_prof.(param.CTD_T),data_prof.P_slow,'k','LineWidth',1)
    leg{end+1}='CTD-T';
    legend(leg);
    set(gca,'ydir','reverse')
    xlabel('Temperature [°C]')
    ylabel('Pressure [dbar]')
    title(gca,'Before calibration')

    subplot(1,3,2) % After calibration
    if param.config.T1
        plot(data_prof.T1_fast,data_prof.P_fast,'LineWidth',1)
        hold on
    end
    if param.config.T2
        plot(data_prof.T2_fast,data_prof.P_fast,'LineWidth',1)
        hold on
    end
    plot(data_prof.(param.CTD_T),data_prof.P_slow,'k','LineWidth',1)
    set(gca,'ydir','reverse')
    xlabel('Temperature [°C]')
    title(gca,'After calibration')

    subplot(1,3,3) % x-x comparison after calibration
    leg3={};
    p3=[];
    if param.config.T1
        p3(end+1)=plot(data_prof.(param.CTD_T),interp1(data_prof.P_fast,data_prof.T1_fast,data_prof.P_slow),'.');
        hold on
        leg3{end+1}='T1';
    end
    if param.config.T2
        p3(end+1)=plot(data_prof.(param.CTD_T),interp1(data_prof.P_fast,data_prof.T2_fast,data_prof.P_slow),'.');
        leg3{end+1}='T2';
    end 
    plot(get(gca,'XLim'),get(gca,'XLim'),'--k','LineWidth',1)
    xlabel('CTD-T [°C]')
    ylabel('FP07-T [°C]')
    legend(p3,leg3)

    
    if ~exist([folder_out,'Figures'],'dir')
        mkdir([folder_out,'Figures'])
    end
    saveas(gcf,[folder_out,'Figures\FP07_calibration.fig'])
    exportgraphics(gcf,[folder_out,'Figures\FP07_calibration.png'],'Resolution',400)

end

end