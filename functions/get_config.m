function [data_prof,filename,cfgfile_mod] = get_config(param,data_prof,folder_out,modify_cfg,kf)
%GET_CONFIG Get the configuration file from the P file and modify it if necessary + patch it to the .P file
%   Detailed explanation goes here
if modify_cfg
    disp('>>> Modification of the configuration file...');
    filename0 = param.filename_list{kf};
    filename = [filename0,'_patched'];
    original_data_file=[param.folder,filename0,'.P'];
    modified_data_file=[param.folder,filename,'.P'];
    if exist(modified_data_file)
        delete(modified_data_file)
    end
    copyfile(original_data_file,modified_data_file); % Create a copy of the P file where the configuration file will be modified
    cfgfile_mod=[folder_out '/setup_' filename0 '_modified'];
    if exist([cfgfile_mod '.cfg']) % Delete any existing modified config file
        delete([cfgfile_mod '.cfg'])
    end
    if isfield(param,'cfgfile')&&~isempty(param.cfgfile) % Get the configuration from the specified configuration file
        patch_setupstr([param.folder,filename],param.cfgfile);
        copyfile([param.cfgfile '.cfg'],[cfgfile_mod '.cfg']); fileattrib([cfgfile_mod '.cfg'],'+w');
    else % Use the configuration specified in the parameters
        % Extract the configuration file from the P file and save it as
        % a temporary file:
        extract_setupstr([param.folder filename '.P'], [cfgfile_mod '.cfg']);

        % Fix the Sensitivity of the shear probe
        A = regexp(fileread([cfgfile_mod '.cfg']),'\n','split'); % Cell array with each cell = a row of the config file
        ind_sh = find(contains(A,'sens        =')); % Indices of the lines where sensitivity is saved
        A(ind_sh(1)) =  {['sens        = ', num2str(param.S_sh1)]};
        if length(ind_sh)>1
            A(ind_sh(2)) =  {['sens        = ', num2str(param.S_sh2)]};
        end

        % Fix the P offset
        ind_Pname = find(contains(A,'name        = P'),1);
        ind_P = find(contains(A,'coef0       ='));
        ind_P=ind_P(ind_P>ind_Pname);
        row_offset=A{ind_P(1)};
        ind0=strfind(row_offset,'=')+1;
        if length(ind0)>1
            ind0=ind0(1);
        end
        indf=strfind(row_offset,';')-1;
        if length(indf)>1
            indf=indf(1);
        end
        original_offset=str2double(row_offset(ind0:indf));
        if isfield(param,'offset_P') && ~isnan(param.offset_P)
            A(ind_P(1)) =  {['coef0       = ', num2str(original_offset - param.offset_P)]};
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
            A(ind_P(1)) =  {['coef0       = ', num2str(original_offset - offset_P)]};
        end
        fid = fopen([cfgfile_mod '.cfg'], 'w');
        fprintf(fid,'%s\n',A{:});
        fclose(fid);

        patch_setupstr([param.folder,filename],cfgfile_mod);
    end

    %% Reconvert to physical units
    default_parameters=odas_p2mat;
    data_prof=odas_p2mat([param.folder,filename,'.P'],default_parameters);
    data_prof.tnum_slow=datenum([data_prof.date,' ',data_prof.time],'yyyy-mm-dd HH:MM:SS.FFF')+data_prof.t_slow/86400;
    data_prof.tdate_slow=datetime(data_prof.tnum_slow,'ConvertFrom','datenum');
    data_prof.tnum_fast=datenum([data_prof.date,' ',data_prof.time],'yyyy-mm-dd HH:MM:SS.FFF')+data_prof.t_fast/86400;
    data_prof.tdate_fast=datetime(data_prof.tnum_fast,'ConvertFrom','datenum');
else % Use the configuration file from the P file
    filename = param.filename_list{kf};
    cfgfile_mod=[folder_out '/setup_' filename];
    if exist([cfgfile_mod '.cfg']) % Delete any existing modified config file
        delete([cfgfile_mod '.cfg'])
    end
    extract_setupstr([param.folder filename '.P'], [cfgfile_mod '.cfg']);
end



end