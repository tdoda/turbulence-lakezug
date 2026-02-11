
function [param] = load_parameters_Zug(lakename,date,general_data_folder,direction,instrument)
%LOAD_PARAMETERS Give parameters specific to a field campaign.
%
%   INPUTS:
%   lakename (str): name of the lake (options: "Geneva", "Zug", "default")
%   date (str): date of the campaign with the format "yyyymmdd"
%   general_data_folder (str): path to the data folder
%   direction (str) [optional]: profile direction (options: "up", "down")
%   instrument (str) [optional]: name of the instrument (options: "microCTD", "VMP")
%
%   OUTPUTS:
%   param (structure): parameters specific to the campaign
%
% T. Doda, 06.12.2024
%% Default values

param.info.mindur_detect = 30; % Minimim duration of a profile [s]
param.info.minvel_detect=0.1; % Minimum speed for profile detection [m/s]
param.info.minKT = 1; % Minimum wavenumber for temperature spectral integration [cpm]
param.info.minKS = 0.1; % Minimum wavenumber for shear spectral integration [cpm]
param.info.maxKS = 14; % Maximum wavenumber for shear spectral integration [cpm]
param.info.fAA = 90; % ~90%*f_AA, where f_AA = 98 Hz
%param.info.Tmethod = 'B'; % Options: 'B'=Batchelor; 'K'=Kraichnan
param.info.Tspec = 'K';% Options: 'B'=Batchelor; 'K'=Kraichnan
param.info.q = 5.26; %3.7; %5.26; %5.26; %2; %1.5; %3.9; %%% q turbulent parameter
param.info.num_fft = 3; % number of fft segments (with 50% overlapping by default): typically 3 or 5
param.info.overlap_pct = 50; % percentage of overlap between fft segments: typically 50%
param.info.int_range = 'L'; % integration range, options: 'S'=Steinbuck 2009; 'L'=Luketina and Imberger 2001
param.info.time_corr = 'NAS'; % Options: 'RSI'; 'KOC'= Kocsis (tau=tau0*W^{-0.5}), 'NAS'=Nash et al., 1999 (tau=tau0*W^{-0.12})
param.info.time_res = 0.0058; %0.0035;  0.0058   %%% 0.0-> no time response correction. Used only if time_corr='KOC' or 'NAS'
param.info.Nasmyth_spec = 'EPFL';% Options: 'ODAS'=default by RSI; 'EPFL'=Bieito's version
param.info.noise_corr = 'Goodman'; % Type of denoising correction. Options: none or 'Goodman'
param.info.npoles = 'single';% Options: 'single' or 'double' -> single-pole or double-pole time response correction of FP07
param.info.kmax_factor = 1/1.66;
param.info.despike_sh  = [ 8  0.5 0.04]; % Parameters to remove spikes in shear data (see odas function "despike")
param.info.despike_A = [8 0.5000 0.0400]; % Parameters to remove spikes in accelerometer data (see odas function "despike")
param.info.ksfact=0.04; % Upper bound of the inertial-convective subrange, as in Steinbuck et al. (2009)
param.info.Snfact=1.55; % Maximum acceptable signal to noise ratio 1.55, as in Goto et al. (2016)

%% Campaign-specific values (can overwrite default values)
if strcmp(lakename,"Zug")
    param=load_campaign_Zug(param,date,general_data_folder);
elseif strcmp(lakename,"default")
    if nargin<4 % Direction was not specified
        direction='down'; % Downward profile by default
        warning('Direction was not specified: use downward direction as default.')
    end
    if nargin<5 % Instrument was not specified
        param.instrument='microCTD';
        warning('Instrument was not specified: use microCTD as default.')
    else
        param.instrument=instrument;
    end
    param=load_campaign_manual(param,direction,general_data_folder);
else
    error("Invalid lake name")
end

%% Instrument-specific values

% About k_HP_cut: low-frequency motions of the free-falling profiler 
% (fluctuations with the scales corresponding
% to half of the instrument length) can be a source of low-frequency noise in the microstructure
% shear signal (see DOI: 10.1016/j.pocean.2006.07.003). 
% MicroCTD length = 1 m --> k_HP_cut = 0.5 cpm (or m-1)
% VMP length = 1.7 m --> k_HP_cut = 0.85 cpm (or m-1)

if strcmp(param.instrument,'VMP')
    param.SNname='028';
    param.CTD_T='SBT';
    param.CTD_C='SBC';
    param.CTD_Chl='FL';
    param.CTD_Turb='NTU';
    param.unit_Chl='ppb';
    param.unit_Turb='FTU';
    param.info.k_HP_cut = 0.85;
    param.space_cfg=false;
elseif strcmp(param.instrument,'microCTD')
    param.SNname='310'; 
    param.CTD_T='JAC_T';
    param.CTD_C='JAC_C';
    param.CTD_Chl='Chlorophyll';
    param.CTD_Turb='Turbidity';
    param.unit_Chl='ppb';
    param.unit_Turb='FTU';
    param.info.k_HP_cut = 0.5;
    param.space_cfg=true;
else
    error('Wrong instrument name')
end

%% T1, T2 parameters
if param.config.T1
    param.info.peak_rem_T1 = [0,0];
    param.info.noisep_T1 = [-10.24,-0.89,param.info.fAA];
    %add some more noise
    info.noisep_T1 = param.info.noisep_T1 + [0.25,0,0];
end

if param.config.T2
    param.info.peak_rem_T2 = [0,0];
    param.info.noisep_T2 = [-10.08,-0.97,param.info.fAA];
    %add some more noise
    info.noisep_T2 = param.info.noisep_T2 + [0.25,0,0];
end
    
    

end


function [param] = load_campaign_manual(param,direction,general_data_folder)
%LOAD_CAMPAIGN_MANUAL Select the data files manually.
%
%   INPUTS:
%   param (structure): default parameters 
%   direction (str): profile direction (options: "up", "down")
%   general_data_folder (str): path to the data folder
%
%   OUTPUTS:
%   param (structure): default parameters 
%
% T. Doda, 29.11.2024

%% Open file selection dialog
[filename, data_folder] = uigetfile('*.P', 'Select a data file',general_data_folder,'MultiSelect','on');
                              
% Check if user selected a file (not canceled)
if data_folder == 0
    error('No file was selected.');
end

% Store file names in a structure without .P extension
if ischar(filename)
    filename=filename(1:strfind(lower(filename),'.p')-1);
    param.filename_list={filename};
else
    for kf=1:length(filename)
        filename{kf}=filename{kf}(1:strfind(lower(filename{kf}),'.p')-1);
    end
    param.filename_list=filename;
end

%param.offset_P=; % Pressure offset (coef 0) if not specified in the configuration file
%param.S_sh1=; % Sensitivity of the shear probe if not specified in the configuration file
%param.S_sh2=; % Sensitivity of the shear probe if not specified in the configuration file
%param.cfgfile = ''; 
param.folder = data_folder;
param.info.system = 'default';

if strcmp(direction,'up')
    param.info.pmin= 0; % Minimum depth [m]
    param.info.pmax = 30; % Maximum depth [m]
    param.info.dp = 0.5; % Bin overlap = bin interval [m]
    param.info.dpD = 1; % Bin size [m]
elseif strcmp(direction,'down')
    param.info.pmin = 1;
    param.info.pmax = 100;
    param.info.dp = 0.5;
    param.info.dpD = 1;
else
    error("Invalid direction")
end

param.info.prof_dir = direction;

%% Sensor configuration:
if strcmp(param.instrument,'VMP')
    param.config.T1=true;
    param.config.T2=true;
    param.config.S1=true;
    param.config.S2=true;
    param.config.uC1=true;
    param.config.uC2=true;
elseif strcmp(param.instrument,'microCTD')
    param.config.T1=true;
    param.config.T2=true;
    param.config.S1=true;
    param.config.S2=false;
    param.config.uC1=true;
    param.config.uC2=false;
else
    error('Wrong instrument name')
end

end

function [param] = load_campaign_Zug(param,date,general_data_folder)
%LOAD_CAMPAIGN_ZUG Give parameters specific to field campaigns in Lake Zug (VMP).
%
%   INPUTS:
%   param (structure): default parameters 
%   date (str): date of the campaign with the format "yyyymmdd"
%   direction (str): profile direction (options: "up", "down")
%   general_data_folder (str): path to the data folder
%
%   OUTPUTS:
%   param (structure): parameters specific to the campaign
%
% T. Doda, 29.11.2024
%% Values for all campaigns
param.info.system = 'Zug';
param.info.pmin = 1; % Minimum depth [m] to start profiles (used in function get_profile). It is only used to perform a first detection of the profiles; a more accurate detection is done in correct_pressure().
param.info.pmax = 180; % Maximum depth [m] of the bin-profiles (not used to detect the profiles, only velocity criterion)
param.info.dp = 0.5; % Bin overlap [m]
param.info.dpD = 1; % Bin size [m]

%% Campaign-specific values
%**************************************************************************
if strcmp(date,'20211110') % Oscar's campaign #1
    param.folder = [general_data_folder,'20211110\Level0\'];
    param.filename_list={'VMP001','VMP002','VMP003','VMP005'}; % Several files can be listed here
    % param.filename_list={'VMP001'};

    % Set P offset and sh probe sensitivity
    % param.offset_P=-0.33;
    % Use shear sensitivities specified in config file: 
    param.cfgfile = 'setup_EAWAG_2018_07_18_OS';

    param.instrument='VMP';
    param.info.prof_dir = 'down';

    param.config.T1=true;
    param.config.T2=true;
    param.config.S1=true;
    param.config.S2=true;
    param.config.uC1=true;
    param.config.uC2=true;
%**************************************************************************
elseif strcmp(date,'20241204') 
    param.folder = [general_data_folder,'20241204\Level0\'];
    param.filename_list={'VMP002','VMP003','VMP004','VMP005','VMP006','VMP008'}; % Several files can be listed here

    % Set P offset and sh probe sensitivity
    % param.offset_P=-0.33;
    % Use shear sensitivities specified in config file: 
    param.cfgfile = 'setup_EAWAG_Zug_2024_12_04';
    param.instrument='VMP';
    param.info.prof_dir = 'down';

    param.config.T1=true;
    param.config.T2=true;
    param.config.S1=true;
    param.config.S2=true;
    param.config.uC1=true;
    param.config.uC2=true;
%**************************************************************************
elseif strcmp(date,'20241205') 
    param.folder = [general_data_folder,'20241205\Level0\'];
    param.filename_list={'VMP001','VMP002','VMP003','VMP004','VMP005','VMP006','VMP007','VMP008','VMP010','VMP011'}; % Several files can be listed here
    


    % Set P offset and sh probe sensitivity
    % param.offset_P=-0.33;
    % Use shear sensitivities specified in config file: 
    param.cfgfile = 'setup_EAWAG_Zug_2024_12_05';
    param.instrument='VMP';
    param.info.prof_dir = 'down';

    param.config.T1=true;
    param.config.T2=true;
    param.config.S1=true;
    param.config.S2=true;
    param.config.uC1=true;
    param.config.uC2=true;
%**************************************************************************
elseif strcmp(date,'20250211') 
    param.folder = [general_data_folder,'20250211\Level0\'];
    param.filename_list={'VMP003','VMP004','VMP007','VMP008','VMP009','VMP010','VMP012','VMP013','VMP014'}; % Several files can be listed here

    % Set P offset and sh probe sensitivity
    % param.offset_P=-0.33;
    % Use shear sensitivities specified in config file: 
    param.cfgfile = 'setup_EAWAG_Zug_2025_02_11';
    param.instrument='VMP';
    param.info.prof_dir = 'down';

    param.config.T1=true;
    param.config.T2=true;
    param.config.S1=true;
    param.config.S2=true;
    param.config.uC1=false;
    param.config.uC2=false;
%**************************************************************************
elseif strcmp(date,'20250626') 
    param.folder = [general_data_folder,'20250626\Level0\'];
    param.filename_list={'VMP002','VMP003','VMP004','VMP005','VMP006'}; % Several files can be listed here

    % Set P offset and sh probe sensitivity
    % param.offset_P=-0.33;
    % Use shear sensitivities specified in config file: 
    param.cfgfile = 'setup_EAWAG_Zug_2025_06_26';
    param.instrument='VMP';
    param.info.prof_dir = 'down';

    param.config.T1=true;
    param.config.T2=true;
    param.config.S1=true;
    param.config.S2=true;
    param.config.uC1=false;
    param.config.uC2=false;
%**************************************************************************
elseif strcmp(date,'20251126') 
    param.folder = [general_data_folder,'20251126\Level0\'];
    param.filename_list={'VMP003','VMP004','VMP005','VMP009','VMP010'}; % Several files can be listed here

    % Set P offset and sh probe sensitivity
    param.offset_P=-0.15;
    % Use shear sensitivities specified in config file: 
    param.cfgfile = 'setup_EAWAG_Zug_2025_11_26_updated';

    param.atm_press_method='offset'; % Options: 'cond' (only for upward),'FP07' (only for upward),'offset','min'
    
    param.instrument='VMP';
    param.info.prof_dir = 'down';
    param.config.T1=true;
    param.config.T2=true;
    param.config.S1=true;
    param.config.S2=false;
    param.config.uC1=false;
    param.config.uC2=false;
%**************************************************************************
elseif strcmp(date,'20251127') 
    param.folder = [general_data_folder,'20251127\Level0\'];
    param.filename_list={'VMP002','VMP003','VMP004','VMP005','VMP006'}; % Several files can be listed here

    % Set P offset and sh probe sensitivity
    param.offset_P=-0.15;
    % Use shear sensitivities specified in config file: 
    param.cfgfile = 'setup_EAWAG_Zug_2025_11_27';

    param.instrument='VMP';
    param.info.prof_dir = 'down';
    param.config.T1=true;
    param.config.T2=true;
    param.config.S1=true;
    param.config.S2=false;
    param.config.uC1=false;
    param.config.uC2=false;
%**************************************************************************
elseif strcmp(date,'20260113') 
    param.folder = [general_data_folder,'20260113\Level0\'];
    param.filename_list={'DAT_053','DAT_055','DAT_057','DAT_059'}; % Several files can be listed here
    param.logbook='Logbook_Zug_20260112'; % csv file

    % Set P offset and sh probe sensitivity
    %param.offset_P=0;
    % Use shear sensitivities specified in config file: 
    param.cfgfile = 'SETUP_updated';
    param.atm_press_method='min'; % Options: 'cond' (only for upward),'FP07' (only for upward),'offset','min'
    
    param.instrument='microCTD';
    param.info.prof_dir = 'down';
    param.config.T1=true;
    param.config.T2=true;
    param.config.S1=true;
    param.config.S2=true;
    param.config.uC1=false;
    param.config.uC2=false;
    
    
else
    error("Invalid campaign date")
end






end

