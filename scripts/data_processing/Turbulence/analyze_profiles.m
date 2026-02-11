%ANALYZE_PROFILES Main script to process microstructure data .P files (VMP,
%microCTD).
%
% T.Doda, last version: 09.02.2026
%%

close all
fclose('all');
clear 
clc

%% Parameters to adapt

lakename='Zug'; % Options: 'Zug' or 'default' (see load_parameters_Zug function)
general_data_folder='..\..\..\data\microCTD\'; % Where fieldwork data is stored
odas_folder='..\..\functions\odas_v4.4\'; % Where ODAS functions are stored
functions_folder="..\..\functions\microstructure\"; % Where microstructure functions are stored
date_campaign="20260113"; % Should match the date in "load_parameters" function except if "default" is used

turbulence_analysis=true; % If =true, run the full turbulence analysis, if =false just check the profiles
modify_cfg=true; % Modify the configuration file (if "false", configuration from .P file is used)
calibrate_FP07=true; % Calibrate FP07

% If turbulence_analysis=false (check data)
save_checkdata=true; % If =true, save the "checked" data
save_checkfig=true; % If =true, save the "checked" figures

% If turbulence_analysis=true (run turbulence analysis):
run_quick_look=false; % Apply quick_look function from Rockland (shear dissipation only)
run_dissip=true; % Compute dissipation based on Bieito's and Sebastiano's script
make_plot_prof = true; % Make profile-related plots.
ind_plot_spectra = []; % Indices of bins where spectra should be plotted (temperature and shear spectra).
show_progress=true;

% To avoid user input
erase_folder=true;


%% Add paths
addpath(odas_folder)
addpath(functions_folder) % Add microstructure functions

%% Load metadata
param=load_parameters_Zug(lakename,date_campaign,general_data_folder);
%param=load_parameters_Geneva(lakename,date_campaign,general_data_folder);
param.filename_list={'DAT_053'};

if modify_cfg 
    if (~isfield(param,'cfgfile') || strcmp(param.cfgfile,''))
        if exist('cfg_file','var') && ~strcmp(cfg_file,'')
            param.cfgfile=[param.folder,cfg_file];
        else
            param.cfg_file=''; 
            warning("No configuration file specified: shear sentivities will be imported from the parameters")
            %error('A configuration file must be specified')
        end
    end
end

if ~isfield(param,'space_cfg')
    param.space_cfg=true;
end

if ~isfield(param,'atm_press_method')
    warning("No atm pressure method specified: use offset")
    param.atm_press_method='offset';
end

% Read logbook
if isfield(param,'logbook')
    data_logbook = readtable([param.folder,param.logbook,'.csv'], 'PreserveVariableNames', true);
    add_coord=true;
else
    add_coord=false;
end

%% Analyze each data file

for kf=1:length(param.filename_list)
    disp('')
    disp('*******************************************')
    fprintf("File %s (%d/%d)\n",param.filename_list{kf},kf,length(param.filename_list))
    
    
    %% Create an output folder for each data file and with metadata in name

    if strcmp(param.info.time_corr,'RSI')
        tmp = param.info.time_corr;
    else
        tmp = [param.info.time_corr num2str(param.info.time_res,'%6.4f') ];
    end
    folder_out = [param.filename_list{kf} '_' param.info.prof_dir num2str(param.info.dpD,'%3.1f') '_'  param.info.Tspec num2str(param.info.q,'%3.1f') ...
        '_' param.info.int_range '_' tmp '_' param.info.npoles 'pole_nfft' num2str(param.info.num_fft)   ...
        '_' param.info.Nasmyth_spec '_' param.info.noise_corr '\'];
    folder_L1 = [param.folder '..\Level1\' folder_out];
    folder_L2 = [param.folder '..\Level2\' folder_out];
    if exist(folder_L1, 'dir') || exist(folder_L2, 'dir')
        if ~erase_folder
            gohead=input('>>> Warning: the folders already exist, do you want to remove them and proceed (y/n): ','s');
        else
            gohead='y';
        end
        if strcmpi(gohead,'y')
            rehash()
            if exist(folder_L1, 'dir')
                rmdir(folder_L1,'s')
            end
            if exist(folder_L2, 'dir')
                rmdir(folder_L2,'s')
            end
        elseif strcmpi(gohead,'n')
            disp('>>> Stop!')
            return
        else
            error('>>> Wrong input. Stop')
        end
    end
    mkdir(folder_L1)
    mkdir(folder_L2)
    
    %% 1st conversion to physical units
    default_parameters=odas_p2mat;
    %default_parameters.speed_tau=0.68/0.99999*2/64; % To avoid smoothing W
    data_prof=odas_p2mat_print([param.folder,param.filename_list{kf},'.P'],false,default_parameters);
   
    if turbulence_analysis
        % Save Level 1 data:
        save([folder_L1,'L1_',param.filename_list{kf},'_',param.info.prof_dir,'.mat'],'data_prof','param')
        disp('>>> Level 1 data saved!')
    end
    
    %% Modify and patch the config file (data_prof is reloaded)
    [data_prof,modified_data_file,cfgfile_mod] = get_config(param,data_prof,folder_L2,modify_cfg,kf,param.space_cfg);
    if ~modify_cfg
        param.cfgfile=[param.folder,cfgfile_mod];
    end
    %% Extract the profiles
    % Start and end indices of the profiles
    ind_prof_slow = get_profile(data_prof.P_slow,data_prof.W_slow,param.info.pmin,...
    param.info.minvel_detect,param.info.prof_dir,param.info.mindur_detect,data_prof.fs_slow);
    ind_prof_fast  = get_profile(data_prof.P_fast,data_prof.W_fast,param.info.pmin,...
        param.info.minvel_detect,param.info.prof_dir,param.info.mindur_detect,data_prof.fs_fast);
    [~,Nprf] = size(ind_prof_slow); % Number of profiles
    if Nprf==0
        warning('>>> No profile for %s: data not saved',filename0)
        continue
    end
    
    %% Calibration of the fast thermistors (data_prof is reloaded)

    if calibrate_FP07
        [data_prof,cfgfile_cal]=run_calibration_FP07(param,data_prof,param.CTD_T,ind_prof_slow,Nprf,modified_data_file,cfgfile_mod,make_plot_prof,folder_L2);
    end

    DATA(kf)=data_prof;
    
    %% Default data analysis from RSI on the modified .P file
    DISS_QL={};

    if turbulence_analysis && run_quick_look
        fprintf(">>> Quick look at the data\n")
        % Pressure limits to show spectra with quick_look
        pmin_ql=param.info.pmin; % min pressure [dbar]
        pmax_ql=param.info.pmax; % max pressure [dbar]
        ql_info=quick_look; % Default parameters
        ql_info.despike_A=param.info.despike_A;
        ql_info.despike_sh=param.info.despike_sh;
        ql_info.profile_min_P=param.info.pmin;
        ql_info.profile_min_W=param.info.minvel_detect;
        ql_info.profile_min_duration=param.info.mindur_detect;  
        if strcmp(param.info.prof_dir,'up')
            ql_info.vehicle='rvmp'; % Upward profiler 
        else
            ql_info.vehicle='vmp'; % Downward profiler
        end
        % Other parameters related to fft are kept as default.
        
        if ~make_plot_prof
            ql_info.make_figures=false; % No plot
        end
        
        for kprof=1:Nprf
            fprintf(">>>>> Profile %d/%d\n",kprof,Nprf)
            ql_info.profile_num=kprof;
            DISS_QL{kprof}=quick_look(modified_data_file,pmin_ql,pmax_ql,ql_info);
            if isempty(DISS_QL{kprof}) % The profile number specified is larger than the number profiles detected by QUICK_LOOK
                continue % Leave empty cell
            end
            DISS_QL{kprof}.approx_bin_size =  DISS_QL{kprof}.diss_length/DISS_QL{kprof}.fs_fast * nanmean (DISS_QL{kprof}.speed);  % in m
 
            if make_plot_prof
                if ~exist([folder_L2,'Figures'],'dir')
                    mkdir([folder_L2,'Figures'])
                end
                figHandles = findobj('Type', 'figure');
                for kfig=1:length(figHandles)
                    saveas(figHandles(kfig),[folder_L2,'Figures\ql_prof',num2str(kprof),'_fig',num2str(kfig),'.fig'])
                    exportgraphics(figHandles(kfig),[folder_L2,'Figures\ql_prof',num2str(kprof),'_fig',num2str(kfig),'.png'],'Resolution',400)
                end
            end
            close all
        end

        
        %save([folder_L2,'quick_look_data.mat'],'DISS_QL')

        clear pmin_ql pmax_ql diss_prof
    end
        
    
    %% Add variables to data_prof
    data_prof.tnum_slow=datenum([data_prof.date,' ',data_prof.time],'yyyy-mm-dd HH:MM:SS.FFF')+data_prof.t_slow/86400;
    data_prof.tdate_slow=datetime(data_prof.tnum_slow,'ConvertFrom','datenum');
    data_prof.tnum_fast=datenum([data_prof.date,' ',data_prof.time],'yyyy-mm-dd HH:MM:SS.FFF')+data_prof.t_fast/86400;
    data_prof.tdate_fast=datetime(data_prof.tnum_fast,'ConvertFrom','datenum');

    data_prof.filename=param.filename_list{kf};

    data_prof.ind_prof_slow_initial=ind_prof_slow;
    data_prof.ind_prof_fast_initial=ind_prof_fast;
    
    %% Compute salinity and density 
    [data_prof.rhoTS,data_prof.Cond_corr,data_prof.Cond_20,data_prof.Sal,~] = compute_rho_salinity(lakename,data_prof.(param.CTD_T),...
        data_prof.(param.CTD_C),data_prof.P_slow,true);

    if param.config.T1
        [data_prof.rhoT1S,data_prof.CondT1_corr,data_prof.CondT1_20,data_prof.SalT1,~] = compute_rho_salinity(lakename,data_prof.T1_fast,...
        interp1(data_prof.P_slow,data_prof.(param.CTD_C),data_prof.P_fast,'linear','extrap'),data_prof.P_fast,true);
    end
    
    if param.config.T2
        [data_prof.rhoT2S,data_prof.CondT2_corr,data_prof.CondT2_20,data_prof.SalT2,~] = compute_rho_salinity(lakename,data_prof.T2_fast,...
        interp1(data_prof.P_slow,data_prof.(param.CTD_C),data_prof.P_fast,'linear','extrap'),data_prof.P_fast,true);
    end
    
    % Save raw pressure series (not corrected)
    data_prof.P_slow_raw=data_prof.P_slow;
    data_prof.P_fast_raw=data_prof.P_fast;
        
    %% Core of the analyis for each profile
    if turbulence_analysis
        counter=1;
        indremove=[];

        % Extract coordinates of the profiles
        indprof_log=find(strcmp(data_logbook.filename,param.filename_list{kf}));
        if length(indprof_log)~=Nprf
            warning('Not same number of profiles than in logbook: coordinates not extracted')
            add_coord=false;
        end


        for kprof = 1:Nprf
            fprintf('>>> Analysis of profile %d/%d\n',kprof, Nprf)
            
            % Remove profiles that are too short:
            if diff(ind_prof_fast(:,kprof))<2*1024 % Lower than 2*nfft (minimum length for function csd_odas)
                warning('Not enough samples in the profile %d of %s: profile not considered',...
                    kprof,filename0)
                indremove(end+1)=kprof;
                if ~isempty(DISS_QL)
                    DISS_QL(kprof)=[];
                end
                continue
            end
    
            % Correct pressure with respect to the atmospheric pressure for the
            % given profile
            [data_prof,press_atm]=correct_pressure(data_prof,param,ind_prof_slow,ind_prof_fast,kprof,make_plot_prof,folder_L2);
    
    
            % Turbulence analysis
            tic
            [BINNED0{counter},SLOW0{counter}, FAST0{counter}] = ...
                resolve_turbulence(data_prof,kprof,param,folder_L2,counter,make_plot_prof,ind_plot_spectra,run_dissip,show_progress);
            toc  
    
            
    
            % ------------------------------------------------------------------
            % Additional script from Sebastiano:
    
            % fname_sp=[folder_main,'/',filename(1:7),'_profID',num2str(i)];
            % [~]=save_single_profiles( BINNED0{i}, SLOW0{i},  FAST0{i},fname_sp);
            % 
            % %creates a matrix with processed results
            % elements = fieldnames(BINNED0{i});
            % for j = 1:length(elements)
            %     varb = elements{j};
            %     if i == 1
            %         eval(['BINNED.',varb,' = transpose(BINNED0{i}.',varb,');'])
            %     else
            %         eval(['BINNED.',varb,' = [BINNED.',varb,', transpose( BINNED0{i}.',varb,')];'])
            %     end
            % end
            % 
            % elements = fieldnames(SLOW0{i});
            % % Reasonable number of SLOW data. Required to have a homogeneous matrix.
            % % nSLOW_ref = max_depth * speed * Freq
            % % 0.5 m/s = ref profiling speed. 64 Hz = freq slow channel
            % if strcmp(direction,'up')
            %     nSLOW_ref = 40/0.5*64;
            % else
            %     nSLOW_ref = 110/0.5*64;
            % end
            % for j = 1:length(elements)
            %     SLOW_ref.(elements{j})=NaN(nSLOW_ref,1);
            % end
            % for j = 1:length(elements)
            %     varb = elements{j};
            %     if j==4  % first vector is pres (j=4)
            %         nSLOW=length(SLOW0{i}.(varb));
            %         if nSLOW>nSLOW_ref
            %             error('length SLOW > maximum plausible length')
            %             pause
            %         else
            %             for jj=4:(length(elements)-3) % -3 : can't convert cells 'date', 'direction, and 'time'
            %                 SLOW_ref.(elements{jj})(1:nSLOW) = SLOW0{i}.(elements{jj});
            %             end
            %             SLOW0{i} = SLOW_ref;
            %         end
            %     end
            %     if i == 1
            %         eval(['SLOW.',varb,' = (SLOW0{i}.',varb,');'])  % not transposed, because it is defined differently!
            %     else
            %         eval(['SLOW.',varb,' = [SLOW.',varb,',( SLOW0{i}.',varb,')];'])  % not transposed, because it is defined differently!
            %     end
            % end
            % elements = fieldnames(FAST0{i});
            % % Reasonable number of FAST data. Required to have a homogeneous matrix.
            % % nSLOW_ref = max_depth * speed * Freq
            % % 0.5 m/s = ref profiling speed. 512 Hz = freq slow channel
            % if strcmp(direction,'up')
            %     nFAST_ref = 40/0.5*512;
            % else
            %     nFAST_ref = 110/0.5*512;
            % end
            % for j = 1:length(elements)
            %     FAST_ref.(elements{j})=NaN(nFAST_ref,1);
            % end
            % for j = 1:length(elements)
            %     varb = elements{j};
            %     if j==4  % first vector is pres (j=4)
            %         nFAST=length(FAST0{i}.(varb));
            %         if nFAST>nFAST_ref
            %             error('length FAST > maximum plausible length')
            %             pause
            %         else
            %             for jj=4:(length(elements)-3) % -3 : can't convert cells 'date', 'direction, and 'time'
            %                 FAST_ref.(elements{jj})(1:nFAST) = FAST0{i}.(elements{jj});
            %             end
            %             FAST0{i} = FAST_ref;
            %         end
            %     end
            %     if i == 1
            %         eval(['FAST.',varb,' = (FAST0{i}.',varb,');'])  % not transposed, because it is defined differently!
            %     else
            %         eval(['FAST.',varb,' = [FAST.',varb,',( FAST0{i}.',varb,')];'])  % not transposed, because it is defined differently!
            %     end
            % end
            
            % % Comparison epsilon and diffusivity from sh and T
            % plot_comparison(BINNED,i,[filename '_' num2str(inp,'%02d')],folder_main,folder_L2)
            % % Close figures
            % all_figs = findobj(0, 'type', 'figure');
            % close(setdiff(all_figs,99));
            % i= i+1;
        
            % Saving each profile as netCDF file:
            DATA_NC.BINNED=BINNED0{counter};
            DATA_NC.SLOW=SLOW0{counter};
            DATA_NC.FAST=FAST0{counter};
            if ~isempty(DISS_QL)
                DATA_NC.DISS_QL=DISS_QL{counter};
            end

            param_prof=param;
            if add_coord
                param_prof.x_coord=data_logbook.X_m(indprof_log(kprof));
                param_prof.y_coord=data_logbook.Y_m(indprof_log(kprof));
            end
            export_to_netcdf([folder_L2,'..\L2_',param.filename_list{kf},'_',param.info.prof_dir,'_prof',num2str(counter),'.nc'],DATA_NC,param_prof,'L2')
            
            counter=counter+1;
        end
    
       
        % Saving the data as .mat file:
        if Nprf>0
            if exist('BINNED0','var')
                BINNED=BINNED0;
                SLOW=SLOW0;
                FAST=FAST0;
                clear BINNED0 SLOW0 FAST0
            else
                BINNED={};
                SLOW={};
                FAST={};
            end
            save([folder_L2,'L2_',param.filename_list{kf},'_',param.info.prof_dir,'.mat'],'BINNED','SLOW','FAST','DISS_QL','param')
            disp('>>> Level 2 data saved!')
            % Comparison epsilon and diffusivity from sh and T
            % plot_comparison(BINNED,[1:length(inPall)],[filename '_all'],folder_main,folder_L2);close all;
        else
            warning('No profile for %s: data not saved',param.filename_list{kf})
        end

    else % Quick check of the profiles
        check_data(param,data_prof,save_checkdata,save_checkfig)

    end
        
end


