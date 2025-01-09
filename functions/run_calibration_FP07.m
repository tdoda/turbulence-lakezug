function [data_prof,filename,cfgfile_mod] = run_calibration_FP07(param,data_prof,CTD_T,ind_prof,Nprf,kf,cfgfile_mod)
%RUN_CALIBRATION_FP07 Calibrate FP07
%   Detailed explanation goes here

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
if param.config.T1
    [T_01,beta1,lag1] = cal_FP07_in_situ_EPFL(data_prof,m,CTD_T,'T1',cfgfile_mod,temp_info);
end

if param.config.T2
    [T_02,beta2,lag2] = cal_FP07_in_situ_EPFL(data_prof,m,CTD_T,'T2',cfgfile_mod,temp_info);
end


%% Re-pacth the P-file 
filename0 = param.filename_list{kf};
filename = [filename0,'_patched'];
original_data_file=[param.folder,filename0,'.P'];
modified_data_file=[param.folder,filename,'.P'];
if exist(modified_data_file,'file')
    delete(modified_data_file)
end
copyfile(original_data_file,modified_data_file); % Create a copy of the P file where the configuration file will be modified

patch_setupstr([param.folder,filename],cfgfile_mod); % patch the new cfg file
%delete([cfgfile_tmp,'.cfg']) % remove the new cfg file

%% Re-convert the data to physical units
clear data_prof
default_parameters=odas_p2mat;
data_prof = odas_p2mat([param.folder,filename,'.P'],default_parameters);          % re-convert data to physical units
if ~strcmp(data_prof.input_parameters.gradT_method,'high_pass')
    error('Error: the gradT_method should be high_pass (if first_difference, pass the info to get_scalar_spectra_odas)')
end

end