function [data_prof,modified_data_file,cfgfile_mod] = get_config(param,data_prof,folder_out,modify_cfg,kf,space_cfg)
%GET_CONFIG Get the configuration file from the P file and modify it if necessary 
% (pressure offset and shear sensitivities) + patch it to the .P file
%  
%   INPUTS:
%   param (structure): parameters specific to the campaign.
%   data_prof (structure): profiling data from a given data file (output of
%   odas_p2mat).
%   folder_out (string): path to the output folder.
%   modify_cfg (boolean): = True to modify cfg file. 
%   kf (integer): index of the data file among the list of files in param.
%
%   OUTPUTS:
%   data_prof (structure): profiling data from a given data file,
%   re-imported with odas_p2mat using the modified cgf file.
%   modified_data_file (string): path to the .P file with modified
%   configuration
%   cfgfile_mod (string): path to the modified .cfg file 
%   space_cfg [optional]: = True if spaces between name and values in cfg
% file. Default=True.
%
% T. Doda based on S. Piccolroaz, last version: 09.02.2026
%% Set up parameters
if nargin<6
    space_cfg=true;
end

if space_cfg
    equal_charact='        = ';
else
    equal_charact='=';
end

filename0 = param.filename_list{kf};
filename = [filename0,'_patched'];
original_data_file=[param.folder,filename0,'.P'];
modified_data_file=[folder_out,filename,'.P'];

cfgfile_original=[folder_out 'setup_' filename0 '_original' '.cfg'];

extract_setupstr(original_data_file, cfgfile_original); % Save the config file from the original .P file

if exist(modified_data_file,'file')
    delete(modified_data_file)
    warning('Delete already existing patched data file')
end
copyfile(original_data_file,modified_data_file); % Create a copy of the P file where the configuration file will be modified

if modify_cfg
    disp('>>> Modification of the configuration file...'); 
    cfgfile_mod=[folder_out 'setup_' filename0 '_modified'];
    
    
    
    if exist([cfgfile_mod '.cfg'],'file') % Delete any existing modified config file
        delete([cfgfile_mod '.cfg'])
        warning('Delete already existing modified config file')
    end
    if isfield(param,'cfgfile')&&~isempty(param.cfgfile) % Get the configuration from the specified configuration file
        if strcmp(param.atm_press_method,'offset')
            answ=input("WARNING: the pressure offset will only be applied if it has been included in the config file. Continue (Y/N)?",'s');
            if ~strcmp(answ,"Y") && ~strcmp(answ,"y")
                error("Stop the analysis: please update the config file or change the pressure correction method")
            end
        end
        patch_setupstr(modified_data_file,[param.folder,param.cfgfile]);
        copyfile([param.folder param.cfgfile '.cfg'],[cfgfile_mod '.cfg']); fileattrib([cfgfile_mod '.cfg'],'+w');
    else % Use the configuration specified in the parameters
        % Extract the configuration file from the P file and save it as
        % a temporary file:
        extract_setupstr(modified_data_file, [cfgfile_mod '.cfg']);


        
        % Fix the Sensitivity of the shear probe
        A = regexp(fileread([cfgfile_mod '.cfg']),'\n','split'); % Cell array with each cell = a row of the config file
        ind_sh = find(contains(A,['sens', equal_charact])); % Indices of the lines where sensitivity is saved
        if isempty(ind_sh) 
            error("No shear sensitivity found in the configuration file")
        end
        if isfield(param,'S_sh1')
            A(ind_sh(1)) =  {['sens',equal_charact, num2str(param.S_sh1)]};
        else
            error('Shear sentivity S_sh1 must be specified in the parameters')
        end
       
        if length(ind_sh)>1
            if isfield(param,'S_sh2')
                A(ind_sh(2)) =  {['sens',equal_charact, num2str(param.S_sh2)]};
            else
                error('Shear sentivity S_sh2 must be specified in the parameters')
            end
        end

        % Fix the P offset
        if strcmp(param.atm_press_method,'offset')
            ind_Pname = find(contains(A,['name',equal_charact,'P']),1);
            ind_P = find(contains(A,['coef0',equal_charact]));
            ind_P=ind_P(ind_P>ind_Pname);
            row_offset=A{ind_P(1)};
            ind0=strfind(row_offset,'=')+1;
            if length(ind0)>1
                ind0=ind0(1);
            end
            indf=strfind(row_offset,';')-1;
            if length(indf)>1
                indf=indf(1);
            elseif isempty(indf)
                indf=length(row_offset);
            end
            original_offset=str2double(row_offset(ind0:indf));
            if isfield(param,'offset_P') && ~isnan(param.offset_P)
                A(ind_P(1)) =  {['coef0',equal_charact, num2str(original_offset - param.offset_P)]};
            else % Ask the user for the offset value
                % Plot the pressure data
                minP=min(data_prof.P_slow);
    
                figure
                plot(data_prof.tdate_slow,data_prof.P_slow);
                hold on
                plot([data_prof.tdate_slow(1),data_prof.tdate_slow(end)],[minP,minP],'-r');
                ylabel("Pressure [dbar]")
                offset_P=[];
                while isempty(offset_P)
                    offset_P=input(sprintf(">>> Pressure offset to subtract (e.g., minP = %f):",minP));
                end
                A(ind_P(1)) =  {['coef0',equal_charact, num2str(original_offset - offset_P)]};
            end
        end
        fid = fopen([cfgfile_mod '.cfg'], 'w');
        fprintf(fid,'%s\n',A{:});
        fclose(fid);

        patch_setupstr(modified_data_file,cfgfile_mod);
    end

    %% Reconvert to physical units
    default_parameters=odas_p2mat;
    data_prof=odas_p2mat_print(modified_data_file,false,default_parameters); % Do not print results (TD, 20251124)
    data_prof.tnum_slow=datenum([data_prof.date,' ',data_prof.time],'yyyy-mm-dd HH:MM:SS.FFF')+data_prof.t_slow/86400;
    data_prof.tdate_slow=datetime(data_prof.tnum_slow,'ConvertFrom','datenum');
    data_prof.tnum_fast=datenum([data_prof.date,' ',data_prof.time],'yyyy-mm-dd HH:MM:SS.FFF')+data_prof.t_fast/86400;
    data_prof.tdate_fast=datetime(data_prof.tnum_fast,'ConvertFrom','datenum');
else % Use the configuration file from the P file
    % filename = param.filename_list{kf};
    % cfgfile_mod=[folder_out 'setup_' filename];
    % if exist([cfgfile_mod '.cfg'],'file') % Delete any existing modified config file
    %     delete([cfgfile_mod '.cfg'])
    % end
    % extract_setupstr(modified_data_file, [cfgfile_mod '.cfg']);

    cfgfile_mod=cfgfile_original;
end



end