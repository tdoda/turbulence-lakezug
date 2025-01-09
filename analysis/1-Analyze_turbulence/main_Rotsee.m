tic
%%% THIS File was made on 3 Feb 2020, the goal is to process all
%%% the data with the same method
%%%
%%% The file was edited in May-June 2021 by Francesco Saturnino
close all
clear
clc

% add path of odas functions and supplementary functions if not in matlab
% path
addpath('./odas_v4.4/odas /')
addpath(genpath('./my_functions/'))

%specify the period you want to analyze
%%%%% CHANGE HERE %%%%%
period = '20191106-07'
direction = 'down'
%%%%%%%%%%%%%%%%%%%%%%%
param=load_period_Rotsee(period,direction);

% Extract all the variables of param:
param_names=fieldnames(param);
for f=1:length(param_names)
    eval([param_names{f},'=param.(param_names{f});']);
end

nfl=length(filename_list);
for i_fl=1:nfl
    %% Setup and coefficients
    info.minvel_detect = 0.25;  % This is only used for FP07 calibration. For the analysis the 90% of the median speed will be used
    info.pmin = minimum_depth;
    info.pmax = maximum_depth;              %%% CHANGE: maximum depth
    info.prof_dir = direction;              %%% CHANGE: up or down
    info.fAA = 90;                          % ~90%*f_AA, where f_AA = 98 Hz
    info.Tmethod = 'B';
    info.Tspec = 'K';                       %%% 'B'=Batchelor; 'K'=Kraichnan
    info.q = 5.26; %3.7; %5.26; %5.26; %2; %1.5; %3.9;          %%% q turbulent parameter
    info.num_fft = 3;                       %%% number of fft lengths (with 50% overlapping): typically 3 or 5
    info.system = 'Rot';
    info.int_range = 'L';                   %%% integration range: 'S' Steinbuck 2009 'L' Luketina and Imberger 2001
    info.time_corr = 'NAS';                 %%% RSI, KOC as Kocsis (tau=tau0*W^{-0.5}), NAS as Nash et al., 1999 (tau=tau0*W^{-0.12})
    info.time_res = 0.0058; %0.0035;  0.0058    %%% 0.0-> no time response correction. Used only if time_corr='KOC' or 'NAS'
    info.Nasmyth_spec = 'EPFL';             %%% ODAS-> default by RSI; EPFL-> Bieito's version
    info.noise_corr = 'Goodman';            %%% none or Goodman
    info.npoles = 'single';                 %%% single or double -> single-pole or double-pole correction of FP07
    suffix = 'def__';                            %%% (optional) suffix to add at the end of folder_out
    PLOT = 0;                           %%% flag: spectra T and sh
    
    %% Create the main output folder
    if (i_fl==1)
        if strcmp(info.time_corr,'RSI')
            tmp = info.time_corr;
        else
            tmp = [info.time_corr num2str(info.time_res,'%6.4f') ];
        end
        folder_main = [folder info.prof_dir num2str(info.dpD,'%3.1f') '_'  info.Tspec num2str(info.q,'%3.1f') ...
            '_' info.int_range '_' tmp '_' info.npoles 'pole_nfft' num2str(info.num_fft)   ...
            '_' info.Nasmyth_spec '_' info.noise_corr suffix];
        if exist(folder_main, 'dir')
            gohead=input('Warning: the folder already exists, do you want to remove it and proceed (y/n): ','s');
            if strcmpi(gohead,'y')
                rmdir(folder_main,'s')
            elseif strcmpi(gohead,'n')
                gohead=input('Warning: do you want to proceed anyway without removing the folder (y/n): ','s');
                if strcmpi(gohead,'n')
                    error('Stop')
                elseif ~strcmpi(gohead,'y')
                    error('Wrong input. Stop')
                end
            else
                error('Wrong input. Stop')
            end
        end
        mkdir(folder_main)
    end
    
    %% Open the P-file and patch the config file
    filename0 = filename_list{i_fl}; % 'DAT_076';
    filename = [filename0,'_patched'];
    a=[folder,filename0,'.P']; b=[folder,filename,'.P']; copyfile(a,b);
    cfgfile_tmp=[folder_main '/setup_tmp'];
    if ~isempty(cfgfile)
        patch_setupstr([folder,filename],cfgfile);
        copyfile([cfgfile '.cfg'],[cfgfile_tmp '.cfg']); fileattrib([cfgfile_tmp '.cfg'],'+w');
    else
        if exist([cfgfile_tmp '.cfg'])
            delete([cfgfile_tmp '.cfg'])
        end
        if strcmp(period,'test_response') && ismember(filename0,{'DAT_082','DAT_084','DAT_041'})
            patch_setupstr([folder,filename],'setup_files/setup_SN310_shM1731_Cis');
            copyfile('setup_files/setup_SN310_shM1731_Cis.cfg',[cfgfile_tmp '.cfg']); fileattrib([cfgfile_tmp '.cfg'],'+w');
        else % Use the configuration file from the P file
            disp('Extracting and modifying configuration file')
            extract_setupstr([folder filename '.P'], [cfgfile_tmp '.cfg']);
            % Fix the Sensitivity of the sh probe and P offset
            A = regexp(fileread([cfgfile_tmp '.cfg']),'\n','split');
            idx_sh = find(contains(A,'sens        ='));
            A(idx_sh(1)) =  {['sens        = ', num2str(S_sh)]};
            if ~isnan(offset_P)
                idx_P = find(contains(A,'coef0       ='));
                A(idx_P(1)) =  {['coef0       = ', num2str(-0.08276 - offset_P)]};
            end
            fid = fopen([cfgfile_tmp '.cfg'], 'w');
            fprintf(fid,'%s\n',A{:});
            fclose(fid);
        end
    end
    
    %% Conversion to Physical units
    default_parameters=odas_p2mat;
    default_parameters.speed_tau=0.68/0.99999*2/64;     % I do no want to smooth W
    DATA=odas_p2mat([folder,filename,'.P'],default_parameters);
    
    %% Calibration of the fast thermistors
    iPs0 = get_profile(DATA.P_slow,DATA.W_slow,info.pmin,...
        info.minvel_detect,info.prof_dir,info.mindur_detect,DATA.fs_slow);
    [~,Nprf] = size(iPs0);
    if Nprf==0
        warning('No profile for %s: data not saved',filename0)
        continue
    end
    m=[]; % position of values along a profile
    
    %     if strcmp(filename0,'VMP_018')
    %         if strcmp(direction,'up')
    %             for i=5:Nprf
    %                 m = [m iPs0(1,i):iPs0(2,i)];
    %             end
    %         else
    %             for i=1:3
    %                 m = [m iPs0(1,i):iPs0(2,i)];
    %             end
    %         end
    %     else
    for i=1:Nprf
        m = [m iPs0(1,i):iPs0(2,i)];
    end
    %     end
    
    temp_info = cal_FP07_in_situ;
    if max(DATA.JAC_T(m))-min(DATA.JAC_T(m))<=8
        temp_info.order=1;
    else
        temp_info.order=2;
    end
    if any(strcmpi(fieldnames(DATA),'JAC_T'))  % MicroCTD
        [T_01,beta1,lag1] = cal_FP07_in_situ_EPFL(DATA,m,'JAC_T','T1',cfgfile_tmp,temp_info);
        [T_02,beta2,lag2] = cal_FP07_in_situ_EPFL(DATA,m,'JAC_T','T2',cfgfile_tmp,temp_info);
    end
    %% Re-pacth the P-file and re-convert the data to physical units
    a=[folder,filename0,'.P']; b=[folder,filename,'.P'];
    delete(b)                                           % remove patched file
    copyfile(a,b);
    patch_setupstr([folder,filename],cfgfile_tmp); % patch the new cfg file
    delete([cfgfile_tmp,'.cfg'])                                 % remove the new cfg file
    clear DATA
    DATA = odas_p2mat([folder,filename,'.P'],default_parameters);          % re-convert data to physical units
    if ~strcmp(DATA.input_parameters.gradT_method,'high_pass')
        error('Error: the gradT_method should be high_pass (if first_difference, pass the info to get_scalar_spectra_odas)')
    end
    
    %% Move to the P-file folder where figures will be saved, and remove existing png and mat files for filename0
    fileList = [dir(['*',filename0, '*.png']);dir(['*',filename0, '*.mat'])];
    if ~isempty(fileList)
        delete(fileList.name)
    end
    
    %% Create folder for output files
    folder_out = [folder_main '/' filename0];
    if exist(folder_out, 'dir')
        rmdir(folder_out,'s')
    end
    mkdir(folder_out)
    
    %% Get profiles
    % 50% of the mean speed. Note: m is the index of the profiles found above.
    info.minvel_detect = 0.5*mean(abs(DATA.W_slow(m)));
    iPs0 = get_profile(DATA.P_slow,DATA.W_slow,info.pmin,info.minvel_detect,info.prof_dir,info.mindur_detect,DATA.fs_slow);    % DATA.W_slow
    iPf0 = get_profile(DATA.P_fast,DATA.W_fast,info.pmin,info.minvel_detect,info.prof_dir,info.mindur_detect,DATA.fs_fast);    % DATA.W_slow
    [~,Nprf] = size(iPs0);
    
    %% Here one can list the profiles that will be extracted for each P-file
    inPall = 1:Nprf;
    
    %% Core of the analyis
    i=1;
    indremove=[];
    for inp = inPall
        fprintf('Profile number %d of %d, ',inp, Nprf)
        
        % Remove profiles that are too short:
        if diff(iPf0(:,inp))<2*1024 % Lower than 2*nfft (minimum length for function csd_odas)
            warning('Not enough samples in the profile %d of %s: profile not considered',...
                inp,filename0)
            indremove(end+1)=inp;
            continue
        end
        % Correction for precise upward profiles:
        % Identify when the sensor exits the lake looking at the
        % discontinuity of the fast conductimeter or, when not available,
        % of the fast thermistor. Use this informatio to correct the
        % pressure. In this way, the Pressure will be relative to the depth
        % of the fast sensors.
        if strcmp(info.prof_dir,'up')
            irange=0.4;    % look for the discontinuity within a range of 0.4 seconds around the current P=0
            iPs0_x=iPs0(2,inp)-round(irange*DATA.fs_slow):iPs0(2,inp)+round(irange*DATA.fs_slow);
            iPf0_x=iPf0(2,inp)-round(irange*DATA.fs_fast):iPf0(2,inp)+round(irange*DATA.fs_fast);
            
            if isfield(DATA,'C1_fast')
                istop=ischange(DATA.gradC1(iPf0_x),'linear');
                istop=find(istop==1,1,'first'); istop=istop-2-ceil(abs(0.001/DATA.W_fast(iPf0_x(istop))*DATA.fs_fast));   % shear probes are 1 mm apart from micro cond. Account also for time response (2 counts to be safe)
            else
                istop=ischange(DATA.gradT1(iPf0_x),'linear','maxnumchanges',1);
                istop=find(istop==1,1,'first'); istop=istop-5-ceil(abs(0.003/DATA.W_fast(iPf0_x(istop))*DATA.fs_fast));   % shear probes are 3 mm apart from micro temp. Account also for the time response: 0.007 s *512 count/s ~ 4 counts (5 counts to be safe)
            end
            
            if isempty(istop)
                spres=0;
            else
                spres = DATA.P_fast(iPf0_x(istop));   % correction
            end
            DATA.P_slow = DATA.P_slow-spres;
            DATA.P_fast = DATA.P_fast-spres;
            iPf0_plot=iPf0_x;
            iPs0_plot=iPs0_x*DATA.fs_fast/DATA.fs_slow;
        else
            spres=0; %In the case of Rotsee, this is 0 as the offset_P is fixed. %min(DATA.P_slow(1:iPs0(1,1)));
            DATA.P_slow = DATA.P_slow-spres;
            DATA.P_fast = DATA.P_fast-spres;
            irange=0.0;
            iPs0_x=iPs0(1,inp)-round(irange*DATA.fs_slow):iPs0(1,inp)+round(irange*DATA.fs_slow);
            iPf0_x=iPf0(1,inp)-round(irange*DATA.fs_fast):iPf0(1,inp)+round(irange*DATA.fs_fast);
            iPf0_plot=iPf0_x;
            iPs0_plot=iPs0_x*DATA.fs_fast/DATA.fs_slow;
            istop=[];
        end
        figure(99)
        set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [29 10], 'PaperPositionMode', 'manual', 'PaperPosition', [0 0 29 10]);
        subplot(1,length(inPall),i)
        plot(iPf0_plot,DATA.P_fast(iPf0_plot),'.-k'); hold on
        if isfield(DATA,'C1_fast')
            plot(iPf0_plot,normalize(DATA.C1_fast(iPf0_plot))/2,'-g');
        end
        plot(iPf0_plot,normalize(DATA.T1_fast(iPf0_plot))/2,'.-r');
        plot(iPf0_plot,DATA.W_fast(iPf0_plot),'.-k');
        plot(iPf0_plot,DATA.sh1(iPf0_plot)/50,'.-m'); hold on
        plot(iPs0_plot,normalize(DATA.JAC_C(iPs0_x))/2,'.-b'); hold on
        plot([iPf0_plot(1) iPf0_plot(end)],[0 0],'--k');
        if ~isempty(istop)
            plot([iPf0_plot(istop) iPf0_plot(istop)],[-1 1],'--g');
        end
        ylim([-1 1])
        
        tic
        [BINNED0{i},SLOW0{i}, FAST0{i}] = resolve_microCTD_profile_all(DATA,inp,info,filename,PLOT,folder_out,i);
        toc
        
        %% uncomment longitude/latitude entries once they are programmed
        %        BINNED0{i}.latitude = lat;
        %        BINNED0{i}.longitude = lon;
        BINNED0{i}.bin_size = info.dpD;
        BINNED0{i}.bin_overlap = info.dp;
        %        FAST0{i}.latitude = lat;
        %        FAST0{i}.longitude = lon;
        %        SLOW0{i}.latitude = lat;
        %        SLOW0{i}.longitude = lon;
        BINNED0{i}.profID=i;
        SLOW0{i}.profID=i;
        FAST0{i}.profID=i;
        FAST0{i}.date = BINNED0{i}.date;
        SLOW0{i}.date = BINNED0{i}.date;
        FAST0{i}.direction = BINNED0{i}.direction;
        SLOW0{i}.direction = BINNED0{i}.direction;
        FAST0{i}.time = BINNED0{i}.time;
        SLOW0{i}.time = BINNED0{i}.time;
        
        fname_sp=[folder_main,'/',filename(1:7),'_profID',num2str(i)];
        [~]=save_single_profiles( BINNED0{i}, SLOW0{i},  FAST0{i},fname_sp);
        
        %creates a matrix with processed results
        elements = fieldnames(BINNED0{i});
        for j = 1:length(elements)
            varb = elements{j};
            if i == 1
                eval(['BINNED.',varb,' = transpose(BINNED0{i}.',varb,');'])
            else
                eval(['BINNED.',varb,' = [BINNED.',varb,', transpose( BINNED0{i}.',varb,')];'])
            end
        end
        
        elements = fieldnames(SLOW0{i});
        % Reasonable number of SLOW data. Required to have a homogeneous matrix.
        % nSLOW_ref = max_depth * speed * Freq
        % 0.5 m/s = ref profiling speed. 64 Hz = freq slow channel
        if strcmp(direction,'up')
            nSLOW_ref = 40/0.5*64;
        else
            nSLOW_ref = 110/0.5*64;
        end
        for j = 1:length(elements)
            SLOW_ref.(elements{j})=NaN(nSLOW_ref,1);
        end
        for j = 1:length(elements)
            varb = elements{j};
            if j==4  % first vector is pres (j=4)
                nSLOW=length(SLOW0{i}.(varb));
                if nSLOW>nSLOW_ref
                    error('length SLOW > maximum plausible length')
                    pause
                else
                    for jj=4:(length(elements)-3) % -3 : can't convert cells 'date', 'direction, and 'time'
                        SLOW_ref.(elements{jj})(1:nSLOW) = SLOW0{i}.(elements{jj});
                    end
                    SLOW0{i} = SLOW_ref;
                end
            end
            if i == 1
                eval(['SLOW.',varb,' = (SLOW0{i}.',varb,');'])  % not transposed, because it is defined differently!
            else
                eval(['SLOW.',varb,' = [SLOW.',varb,',( SLOW0{i}.',varb,')];'])  % not transposed, because it is defined differently!
            end
        end
        elements = fieldnames(FAST0{i});
        % Reasonable number of FAST data. Required to have a homogeneous matrix.
        % nSLOW_ref = max_depth * speed * Freq
        % 0.5 m/s = ref profiling speed. 512 Hz = freq slow channel
        if strcmp(direction,'up')
            nFAST_ref = 40/0.5*512;
        else
            nFAST_ref = 110/0.5*512;
        end
        for j = 1:length(elements)
            FAST_ref.(elements{j})=NaN(nFAST_ref,1);
        end
        for j = 1:length(elements)
            varb = elements{j};
            if j==4  % first vector is pres (j=4)
                nFAST=length(FAST0{i}.(varb));
                if nFAST>nFAST_ref
                    error('length FAST > maximum plausible length')
                    pause
                else
                    for jj=4:(length(elements)-3) % -3 : can't convert cells 'date', 'direction, and 'time'
                        FAST_ref.(elements{jj})(1:nFAST) = FAST0{i}.(elements{jj});
                    end
                    FAST0{i} = FAST_ref;
                end
            end
            if i == 1
                eval(['FAST.',varb,' = (FAST0{i}.',varb,');'])  % not transposed, because it is defined differently!
            else
                eval(['FAST.',varb,' = [FAST.',varb,',( FAST0{i}.',varb,')];'])  % not transposed, because it is defined differently!
            end
        end
        
        % Comparison epsilon and diffusivity from sh and T
        plot_comparison(BINNED,i,[filename '_' num2str(inp,'%02d')],folder_main,folder_out)
        % Close figures
        all_figs = findobj(0, 'type', 'figure');
        close(setdiff(all_figs,99));
        i= i+1;
    end
    inPall(ismember(inPall,indremove))=[];
            saveas(99,[folder_out,'/overview_',num2str(inp,'%02d'),info.prof_dir,'.png']);
        close(99);
        
    if ~isempty(inPall)
        save([folder_main,'/microCTD_20200319_',filename,'_',info.prof_dir],'BINNED','SLOW','FAST')

        % Comparison epsilon and diffusivity from sh and T
        plot_comparison(BINNED,[1:length(inPall)],[filename '_all'],folder_main,folder_out);close all;
    else
        warning('No profile for %s: data not saved',filename0)
    end
end