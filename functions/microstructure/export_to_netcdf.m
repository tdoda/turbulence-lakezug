function export_to_netcdf(filename,DATA,ATT,level)
% EXPORT_TO_NETCDF: Export data from a structure array to netCDF using
% low-level functions (similar to NetCDF C library).
%
%   INPUTS:

%   filename (string): path and name of the output file
%   DATA (structure array): data to export to netCDF. Each field is the
%   name of a variable. If the variable is a structure array, all the variables 
%   from that structure array will be extracted separately and renamed.
%   Only the variables defined in create_param_netCDF
%   fill be exported to the netCDF file.
%   ATT (structure array): attributes to add to the general attributes of
%   the netCDF file. Each field is the name of an attribute. 
%   level (string): data output level, either 'L1' or 'L2'
%
% T. Doda, 20.01.2026
%%

% Create parameters of the netCDF file:
PARAM=create_param_netCDF(level);

% Add general attributes:
PARAM=add_att(PARAM,ATT);

% Add data:
PARAM=add_data(PARAM,DATA);


% Create netCDF with nccreate for each variable
% Create file and open it in "define" mode:

%ncid = netcdf.create(filename,'CLOBBER'); % Overwrite existing file
ncid=netcdf.create(filename,'NETCDF4');
c = onCleanup(@() netcdf.close(ncid)); % Close the file if an error occurs


% Define dimensions:
dim_all=fieldnames(PARAM.dimensions);
dim_names=cell(1,length(dim_all));
dimid=NaN(1,length(dim_all));
diminf=netcdf.getConstant('NC_UNLIMITED');
for kdim=1:length(dim_all)
    dim_names{kdim}=PARAM.dimensions.(dim_all{kdim}).dim_name;
    dimsize=PARAM.dimensions.(dim_all{kdim}).dim_size;
    if isinf(dimsize)||isempty(dimsize)||isnan(dimsize)
        dimid(kdim) = netcdf.defDim(ncid,dim_names{kdim},diminf);
    else
        dimid(kdim) = netcdf.defDim(ncid,dim_names{kdim},dimsize);
    end
end

% Define variables and attributes: 
var_all=fieldnames(PARAM.variables);
varid=NaN(1,length(var_all));
for kvar=1:length(var_all)
    varinfo=PARAM.variables.(var_all{kvar});
    dimcell=varinfo.dim;
    if ~iscell(dimcell)
        dimcell={dimcell};
    end
    dimarray=NaN(1,length(dimcell));
    for kdim=1:length(dimcell)
        
        dimarray(kdim)=dimid(strcmp(dim_names,dimcell{kdim}));

    end

    varid(kvar) = netcdf.defVar(ncid,varinfo.var_name,...
    varinfo.data_type,dimarray);


    % Add attributes:
    attnames=fieldnames(varinfo);
    attnames(strcmp(attnames,'dim'))=[];
    for katt=1:length(attnames)
        netcdf.putAtt(ncid,varid(kvar),attnames{katt},varinfo.(attnames{katt}))
    end
end

% Write global attributes
varglobal=netcdf.getConstant('NC_GLOBAL');
globatt_all=fieldnames(PARAM.general_attributes);
for katt=1:length(globatt_all)
    netcdf.putAtt(ncid,varglobal,globatt_all{katt},PARAM.general_attributes.(globatt_all{katt}))
end


% Change to data mode to write data:
netcdf.endDef(ncid);

% Write data for each variable
for kvar=1:length(var_all)
    datavar = PARAM.data.(var_all{kvar});
    netcdf.putVar(ncid, varid(kvar), 0, length(datavar), datavar);
    
    % try
    % netcdf.putVar(ncid,varid(kvar),PARAM.data.(var_all{kvar}));
    % catch
    %     pause(1)
    % end
end

% netCDF file will automatically be closed with onCleanup function defined
% at the beginning. If onCleanup is not used, add:
%netcdf.close(ncid)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PARAM=create_param_netCDF(level)
%CREATE_PARAM_NETCDF: create general attributes, dimensions and variables
%based on the selected output data level.
% 
%   INPUTS:
%   level (string): data output level, either 'L1' or 'L2'
%
%   OUTPUT:
%   PARAM (structure array): information about variables to export (fields:
%   general_attributes, dimensions, variables). 
%
% T. Doda, 20.01.2026
%%

PARAM.general_attributes=struct("institution","Eawag-Unil",...
            "history","See history on Github",...
            "conventions","CF 1.7",...
            "comment","Data from microstructure profiles (VMP, microCTD) collected by the Eawag Aphys group and the UNIL Lakes group in Lake Zug",...
            "title","Lake Zug microstructure profiles");

% The fields of PARAM.dimensions and PARAM.variables must match the name of
% variables in the DATA structure.
% The fields dim_name and var_name are the names of dimensions and
% variables that will be exported to the netCDF file. They do not have to
% match the variable names from the DATA structure. 
if strcmp(level,'L1') % Level 1 export
    PARAM.dimensions=struct();  
    PARAM.variables=struct();
elseif strcmp(level,'L2') % Level 2 export
    PARAM.dimensions=struct('SLOW_profID',struct('dim_name','profID','dim_size',1),...
        'SLOW_dtsec',struct('dim_name','SLOW_dtsec','dim_size',netcdf.getConstant('NC_UNLIMITED')),...
        'FAST_dtsec',struct('dim_name','FAST_dtsec','dim_size',netcdf.getConstant('NC_UNLIMITED')),...
        'BINNED_depth',struct('dim_name','BINNED_depth','dim_size',netcdf.getConstant('NC_UNLIMITED')));
    
    PARAM.variables=struct('SLOW_datenum',struct('var_name','SLOW_datenum','dim','profID','unit','Number of days since January 1, 0000','long_name','Measurements starting time (UTC)','data_type','NC_DOUBLE'),...
        'SLOW_profID',struct('var_name','profID','dim','profID','unit','','long_name','Profile index','data_type','NC_DOUBLE'),...
        'SLOW_dtsec',struct('var_name','SLOW_dtsec','dim','SLOW_dtsec','unit','s','long_name','Time during profiling since the start of measurements','data_type','NC_DOUBLE'),...
        'SLOW_pressure',struct('var_name','SLOW_pressure','dim','SLOW_dtsec','unit','dbar','long_name','Corrected water pressure','data_type','NC_DOUBLE'),...
        'SLOW_raw_pressure',struct('var_name','SLOW_raw_pressure','dim','SLOW_dtsec','unit','dbar','long_name','Water pressure before any correction for atmospheric pressure','data_type','NC_DOUBLE'),...
        'FAST_dtsec',struct('var_name','FAST_dtsec','dim','FAST_dtsec','unit','s','long_name','Time during profiling since the start of measurements','data_type','NC_DOUBLE'),...
        'FAST_pressure',struct('var_name','FAST_pressure','dim','FAST_dtsec','unit','dbar','long_name','Corrected water pressure','data_type','NC_DOUBLE'),...
        'FAST_raw_pressure',struct('var_name','FAST_raw_pressure','dim','FAST_dtsec','unit','dbar','long_name','Water pressure before any correction for atmospheric pressure','data_type','NC_DOUBLE'),...
        'FAST_depth',struct('var_name','FAST_depth','dim','FAST_dtsec','unit','m','long_name','Water depth','data_type','NC_DOUBLE'),...
        'FAST_fast_T1',struct('var_name','FAST_fast_T1','dim','FAST_dtsec','unit','°C','long_name','Microstructure temperature T1 (calibrated)','data_type','NC_DOUBLE'),...
        'FAST_fast_T2',struct('var_name','FAST_fast_T2','dim','FAST_dtsec','unit','°C','long_name','Microstructure temperature T2 (calibrated)','data_type','NC_DOUBLE'),...
        'FAST_grad_T1',struct('var_name','FAST_grad_T1','dim','FAST_dtsec','unit','°C.m^(-1)','long_name','Vertical gradient of microstructure temperature T1','data_type','NC_DOUBLE'),...
        'FAST_grad_T2',struct('var_name','FAST_grad_T2','dim','FAST_dtsec','unit','°C.m^(-1)','long_name','Vertical gradient of microstructure temperature T2','data_type','NC_DOUBLE'),...
        'FAST_fast_S1',struct('var_name','FAST_fast_S1','dim','FAST_dtsec','unit','s^(-1)','long_name','Shear S1','data_type','NC_DOUBLE'),...
        'FAST_fast_S2',struct('var_name','FAST_fast_S2','dim','FAST_dtsec','unit','s^(-1)','long_name','Shear S2','data_type','NC_DOUBLE'),...
        'BINNED_depth',struct('var_name','BINNED_depth','dim','BINNED_depth','unit','m','long_name','Averaged water depth of the bins','data_type','NC_DOUBLE'),...
        'BINNED_Xi_T1',struct('var_name','BINNED_Xi_T1','dim','BINNED_depth','unit','°C^2.s^(-1)','long_name','Temperature variance dissipation rate from microstructure temperature T1','data_type','NC_DOUBLE'),...
        'BINNED_Xi_T2',struct('var_name','BINNED_Xi_T2','dim','BINNED_depth','unit','°C^2.s^(-1)','long_name','Temperature variance dissipation rate from microstructure temperature T2','data_type','NC_DOUBLE'),...
        'BINNED_eps_T1',struct('var_name','BINNED_eps_T1','dim','BINNED_depth','unit','m^2.s^(-3)','long_name','TKE dissipation rate from microstructure temperature T1','data_type','NC_DOUBLE'),...
        'BINNED_eps_T2',struct('var_name','BINNED_eps_T2','dim','BINNED_depth','unit','m^2.s^(-3)','long_name','TKE dissipation rate from microstructure temperature T2','data_type','NC_DOUBLE'),...
        'BINNED_eps_S1',struct('var_name','BINNED_eps_S1','dim','BINNED_depth','unit','m^2.s^(-3)','long_name','TKE dissipation rate from shear S1','data_type','NC_DOUBLE'),...
        'BINNED_eps_S2',struct('var_name','BINNED_eps_S2','dim','BINNED_depth','unit','m^2.s^(-3)','long_name','TKE dissipation rate from shear S2','data_type','NC_DOUBLE'),...
        'BINNED_KOsborn_T1',struct('var_name','BINNED_KOsborn_T1','dim','BINNED_depth','unit','m^2.s^(-1)','long_name','Osborn diapycnal diffusivity from microstructure temperature T1','data_type','NC_DOUBLE'),...
        'BINNED_KOsborn_T2',struct('var_name','BINNED_KOsborn_T2','dim','BINNED_depth','unit','m^2.s^(-1)','long_name','Osborn diapycnal diffusivity from microstructure temperature T2','data_type','NC_DOUBLE'),...
        'BINNED_KOsborn_S1',struct('var_name','BINNED_KOsborn_S1','dim','BINNED_depth','unit','m^2.s^(-1)','long_name','Osborn diapycnal diffusivity from shear S1','data_type','NC_DOUBLE'),...
        'BINNED_KOsborn_S2',struct('var_name','BINNED_KOsborn_S2','dim','BINNED_depth','unit','m^2.s^(-1)','long_name','Osborn diapycnal diffusivity from shear S2','data_type','NC_DOUBLE'),...
        'BINNED_KOsbornCox_T1',struct('var_name','BINNED_KOsbornCox_T1','dim','BINNED_depth','unit','m^2.s^(-1)','long_name','Osborn-Cox diapycnal diffusivity from microstructure temperature T1','data_type','NC_DOUBLE'),...
        'BINNED_KOsbornCox_T2',struct('var_name','BINNED_KOsbornCox_T2','dim','BINNED_depth','unit','m^2.s^(-1)','long_name','Osborn-Cox diapycnal diffusivity from microstructure temperature T2','data_type','NC_DOUBLE'),...
        'BINNED_LTuT1',struct('var_name','BINNED_LTuT1','dim','BINNED_depth','unit','m','long_name','Thorpe scale computed from bin-averaged Thorpe displacements in microstructure temperature T1','data_type','NC_DOUBLE'),...
        'BINNED_LTuT2',struct('var_name','BINNED_LTuT2','dim','BINNED_depth','unit','m','long_name','Thorpe scale computed from bin-averaged Thorpe displacements in microstructure temperature T2','data_type','NC_DOUBLE'),...
        'BINNED_flag_T1',struct('var_name','BINNED_flag_T1','dim','BINNED_depth','unit','-','long_name','Acceptance flag for turbulence estimates from microstructure temperature T1 (0=keep, 1=reject)','data_type','NC_DOUBLE'),...
        'BINNED_flag_T2',struct('var_name','BINNED_flag_T2','dim','BINNED_depth','unit','-','long_name','Acceptance flag for turbulence estimates from microstructure temperature T2 (0=keep, 1=reject)','data_type','NC_DOUBLE'),...
        'BINNED_flag_S1',struct('var_name','BINNED_flag_S1','dim','BINNED_depth','unit','-','long_name','Acceptance flag for turbulence estimates from shear S1 (0=keep, 1=reject)','data_type','NC_DOUBLE'),...
        'BINNED_flag_S2',struct('var_name','BINNED_flag_S2','dim','BINNED_depth','unit','-','long_name','Acceptance flag for turbulence estimates from shear S2 (0=keep, 1=reject)','data_type','NC_DOUBLE'));
else
    error('Wrong data level')
end



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PARAM=add_data(PARAM,DATA)
%ADD_DATA: add data to the variables specified in the netCDF parameters.
% 
%   INPUTS:
%   PARAM (structure array): information about variables to export (fields:
%   general_attributes, dimensions, variables). 
%   DATA (structure array): data to export to netCDF. Each field is the
%   name of a variable. If the variable is a structure array, all the variables 
%   from that structure array will be extracted separately and renamed.
%   Only the variables defined in create_param_netCDF are added.
%
%   OUTPUT:
%   PARAM (structure array): same as the input PARAM with an additional
%   field "data".
%
% T. Doda, 20.01.2026
%%

VAR=extract_variables(DATA,'structname','numeric');

fnames=fieldnames(VAR);
PARAM.data=struct();
for k=1:length(fnames)
    if ismember(fnames{k},fieldnames(PARAM.variables))
        PARAM.data.(fnames{k})=VAR.(fnames{k});
    end
end

% Remove the variables without data:
PARAM.variables = rmfield(PARAM.variables, setdiff(fieldnames(PARAM.variables), fieldnames(PARAM.data)));



end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PARAM=add_att(PARAM,ATT)
%ADD_ATT: add attributes to the netCDF parameters.
% 
%   INPUTS:
%   PARAM (structure array): information about variables to export (fields:
%   general_attributes, dimensions, variables). 
%   ATT (structure array): attributes to add to the general attributes of
%   the netCDF files. Each field is the name of an attribute. 
%
%   OUTPUT:
%   PARAM (structure array): same as the input PARAM with values of the
%   additional attributes.
%
% T. Doda, 20.01.2026
%%

% Get all the attributes (including fields from a structure, but not cell
% arrays):
VAR=extract_variables(ATT,'structname','');

fnames=fieldnames(VAR);

% Add all the attributes in VAR:
for k=1:length(fnames)
        PARAM.general_attributes.(fnames{k})=num2str(VAR.(fnames{k}));
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function VAR=extract_variables(DATA,prefix_struct,data_format)
%EXTRACT_VARIABLES: keep variables from DATA that are numerical arrays,
%including the numerical fields of stucture arrays in DATA.
% 
%   INPUTS:
%   DATA (structure array): data to export to netCDF. Each field is the
%   name of a variable. If the variable is a structure array, all the variables 
%   from that structure array will be extracted separately and renamed.
%   prefix_struct (string) [optional]: prefix to use for the name of variables stored in a
%   structure, either 'structname' (name of the structure) or '' (no
%   prefix). Default: ''.
%   data_format (string) [optional]: if "numeric", only extract numeric
%   data. Otherwise, extract everything except structures and cell arrays. Default: empty.
%
%   OUTPUT:
%   VAR (structure array): same as DATA containing only numerical arrays.
%
% T. Doda, 20.01.2026
%%
fnames=fieldnames(DATA);
VAR=struct();

if nargin<2
    prefix_struct="";
end

if nargin<3
    data_format="";
end

for k=1:length(fnames)
    if isstruct(DATA.(fnames{k}))
        struct_name=fnames{k};
        fnames_sub=fieldnames(DATA.(struct_name));
        for k_sub=1:length(fnames_sub)
            if strcmp(prefix_struct,'structname')
                varname=[struct_name,'_',fnames_sub{k_sub}];
            else
                varname=fnames_sub{k_sub};
            end
            if strcmp(data_format,'numeric')
                if isnumeric(DATA.(struct_name).(fnames_sub{k_sub}))
                    VAR.(varname)=DATA.(struct_name).(fnames_sub{k_sub});
                end
            elseif ~isstruct(DATA.(struct_name).(fnames_sub{k_sub})) && ~iscell(DATA.(struct_name).(fnames_sub{k_sub}))  
                VAR.(varname)=DATA.(struct_name).(fnames_sub{k_sub});
            end
        end
    elseif strcmp(data_format,'numeric')
        if isnumeric(DATA.(fnames{k}))
            VAR.(fnames{k})=DATA.(fnames{k});
        end
    elseif ~isstruct(DATA.(fnames{k})) && ~iscell(DATA.(fnames{k}))  
        VAR.(fnames{k})=DATA.(fnames{k});
    end
end

end




