%% cal_FP07_in_situ

%% Eawag (Tomy Doda) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   2024-12-06: modified to deal with cfg files that do not contain two
%   beta coefficients (e.g., VMP).
%%  EPFL (Sebastiano Piccolroaz) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   2020-04-22: Function modified allowing to calibrate the fast
%   thermistors without the need to save the mat-file of the data. All
%   plots have been removed. At the end of the script, the cfg file is 
%   updated with the calibrated parameters. 
%   Overall, this script allows to save computational time. The function 
%   has been tested against the original one, and the results are the same.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine calibration coefficients for a thermistor probe using in situ data.
%%
% <latex>\index{Functions!cal\_FP07\_in\_situ}</latex>
%
%%% Syntax
%   [T_0, beta, Lag] = cal_FP07_in_situ( file_name, 
%                           T_ref_string, T_string, SN, cal_info,... )
%
% * [file_name] Name of mat-file containing data used to
%       calibrate a thermistor. Note: Data needs to have been previously 
%       extracted from .P file using either quick_look or odas_p2_mat.
% * [T_ref_string] Name of vector within the mat-file (e.g. 'JAC_T') 
%       that contains the reference temperature in degrees celsius. Usually 
%       from a SBE4F thermometer or a JAC CT. 
% * [T_string] Name of thermistor to calibrate, typically 'T1' or 'T2'.
% * [SN] Serial number of thermistor.
% * [cal_info] Structure (optional) containing configuration parameters. 
%       A template is generated and returned when cal_FP07_in_situ is called 
%       with no input parameters. The parameters are described below.
% * [...] Optional configuration parameters to supplement, or override, 
%       those values included within cal_info. Inputs are accepted as
%       string/value pairs.
% * []
% * [T_0] Value of parameter T_0, used in the Steinhart-Hart equation. When
%       called with no input parameters, a structure containing default
%       input parameter values is returned for reference.
% * [beta] beta coefficients, in ascending order, of the fit to the Steinhart-Hart
%       equation. i.e. beta_1, beta_2, beta_3. 
% * [Lag] Delay, in seconds, between the thermistor and the reference
%       thermometer. Typically a negative value because the reference
%       sensor is usually behind the thermistor being calibrated.
%
%%% Description
%
% This function can be used to calibrate a FP07 thermistor probe using 
% in-situ data. The reference temperature data is usually provided by a Sea-Bird 
% SBE3F or a JAC-CT. The reference temperature data must be contained within 
% the specified mat-file. If necessary the data can be incorportated into the 
% file using a hotel file at the time of data conversion. This is typically 
% necessary when analyzing data collected with a glider. 
%
% This function processes the data as follows:
% (1) Selects the portion of the data file that will be used for the 
% calibration based on the input parameters. The range will be plotted if 
% plot_range = true. 
% (2) It then detrends the thermistor data and scales it so that it is 
% approximately aligned with the reference thermometer. A plot of the
% scaled signals will only be shown if plot_scaled = true. 
% (3) It then computes the cross-correlation coefficient between the thermistor 
% and the reference thermometer and estimates the lag between these two 
% signals. A plot will only be generated if plot_xcorr = true.
% (4) Next it does a regression based on the Steinhart-Hart equation and estimates 
% the coefficients based on the chosen order of the regression. A plot of the 
% natural logarithm of the resistance ratio against the inverse of the 
% absolute temperature is shown if plot_regress = true.
% (5) The regression coefficients are then used to convert the thermistor
% data into physical units. The depth profile of the calibrated thermistor
% data and the reference signal, and their difference, will be plotted if
% plot_result = true. 
%
% If $\texttt{cal\_FP07\_in\_situ}$ is called without input arguments, it 
% returns the default parameters used by $\texttt{cal\_FP07\_in\_situ}$. 
% You can then customize this structure to your particular processing 
% requirements. For example,
%
%    >> cal_info = cal_FP07_in_situ
%    
%    cal_info = 
% 
%             make_figures: 1
%                    order: 2
%               plot_range: 1
%             plot_regress: 1
%              plot_result: 1
%              plot_scaled: 0
%               plot_xcorr: 0
%            profile_min_P: 1
%            profile_min_W: 0.2000
%     profile_min_duration: 20
%              profile_num: 1
%             vehicle_info: []
%
% The configuration parameters (fields) within the structure $\texttt{cal\_info}$ 
% that control the behaviour of $\texttt{cal\_FP07\_in\_situ}$,are listed below. 
% They are grouped for clarity, but are all part of the single structure 
% $\texttt{cal\_info}$.
%
%%% Parameters that control the calibration
%
% * [order] Fit order to the Steinhart-Hart equation. Value can be 1,
%       2, or 3. Default = 2. For small temperature ranges, order = 1 is 
%       recommended.
%
%%% Parameters that specify a profile
% * [profile_num] Index to the requested profile from the set of detected
%       profiles.  The first profile, 1, is the default value.
% * [profile_min_P] The minimum pressure of a profile, [dbar]. Default = 1.
% * [profile_min_W] The minimum vertical speed of profiling  
%       in [dbar/s]. Default = 0.2. Use a smaller value, e.g. 0.1 dBar/s, 
%       for gliders.
% * [profile_min_duration] The minimum duration in which the minimum
%      pressure and speed must be satisfied [s]. Default = 20.
% * [vehicle_info] A structure found in all .mat data files that is used
%       to determine the direction of profiling. Default = [].
%       If empty, the information will be loaded from the datafile. If a 
%       change is desired, the profiling direction needs to be specified, 
%       i.e. Set vehicle_info.profile_dir to 'down' for vmps, 'up' for rvmps
%       or 'glide' for gliders.
%
%
%%%% Parameters that toggle data visualization
%
% * [make_figures] The parameter that determines if figures are generated.
%     make_figures = false suppresses the generation of figures to speed up
%     the data processing. Default = true.
% * [plot_range] A logical parameter that determines if a figure of the 
%     selected is generated. Default = true.
% * [plot_scaled] A logical parameter that determines if a figure of the 
%     detrended signals and lag are plotted. Default = false.
% * [plot_xcorr] A logical parameter that determines if a figure of the 
%      cross correlation is plotted. Default = false.
% * [plot_regress] A logical parameter that determines if the regression
%      data and fit are plotted. Default = true.
% * [plot_result] A logical parameter that determines if figures of the 
%       calibration results (i.e. comparison and difference) are shown. 
%       Default = true.


%___________________________
%
% Version History
%
% * 2013-12-05 (RGL) original version.
% * 2013-12-06 (RGL) added varargin to define a profile with default values.
% * 2015-04-10 (WID) revised documentation for publishing.
% * 2015-04-27 (RGL) modified to allow the specification of the fit order.
% * 2015-07-27 (WID) use texstr in place of fix_underscore.
% * 2015-07-29 (WID) return default values when called with no input
%                       parameters.
% * 2015-10-27 (RGL) Changed description section.
% * 2016-06-07 (RGL) Changed legend call, added clf to the start of figures.
% * 2016-11-10 (RGL) Changed the call to the deconvolve function so that it
%      includes the thermistor signal without pre-emphasis. odas_p2mat uses
%      the the thermistor without pre-emphasis (if it is available) for
%      conversion in to physical units.
%      Consequently, the thermistor signal without pre-emphasis MUST be
%      called with this in situ calibration function. Otherwise, there is a
%      substantial error due to the offset (order ~10 counts) that may be
%      present in the signal with pre-emphasis. This version also includes
%      a test for the existence of the signal without pre-emphasis so that
%      it does not bomb in its absence. I also beautified the legends and
%      labels.
% * 2016-11-10 (RGL) Added a low-pass filter to the thermistor data, in the
%      case of a JAC-T reference, in order to get a tighter regression. The
%      cut-off frequency is sped dependent and follows the recommendation
%      for calculating salinity. There is an improvement when the
%      thermistor is filtered.
% * 2016-11-15 (RGL) Added fc in case of Sea-Bird thermometer.
% * 2017-11-28 (RGL) Added ability to handle both type=therm abd type=t_ms.
% * 2017-11-30 (RGL) Some more display changes and warning will be
%      suppressed in case the channel without pre-emphasis does not exist.
%      Made the naming of parameters (fields) consistent with usage in
%      quick_look.
% * 2019-06-27 (JMM) Changes to function plotting and inputs. Regression
% algorithm is unaffected. Changes include:
%       - Added an initial subplot of the pressure record, highlighting the 
%           portion of the file being used
%       - Load in vehicle_info section from .mat file (datafile is required)
%       - Change default profile_min_W to 0.2m/s to match quick_look inputs
%       - Added plotting flags to be able to supress figures.
%       - Changed variable name from T2 to T_prof to be more descriptive.
%       - Added a warning message if temperature range is too small for high
%         order fit. 
% * 2019-06-03 (JMM) Updated comments and documentation. 


function [T_0,beta,Lag] = cal_FP07_in_situ_EPFL(DATA,m,T_ref_string,T_string,cfgfile,varargin)

%-----------------------------------------------------------------
% ----- Default parameters ---------------------------------------
%-----------------------------------------------------------------
default_vehicle_info         = []; % will trigger profile_dir = 'down'
default_order                = 2; % The order of the fit to the SS equation
default_make_figures         = false; % render figures for data visulization
default_plot_range           = false; % flag to plot the range used for profiles
default_plot_scaled          = false; % flag to plot detrended signals showing lag
default_plot_xcorr           = false; % flag to plot cross-correlation
default_plot_regress         = false; % flag to plot regression
default_plot_result          = false; % flag to plot result of calibration

if ~nargin
    for d = whos('default_*')'
        param = regexp(d.name, 'default_(.*)', 'tokens');
        result.(char(param{1})) = eval(d.name);
    end
    T_0 = result;
    return
end

%-----------------------------------------------------------------
% ----- Parse Inputs ---------------------------------------------
%-----------------------------------------------------------------
p = inputParser;
p.KeepUnmatched = true;
p.CaseSensitive = true;

val_numeric     = @(x) isnumeric(x) && isscalar(x);
val_string      = @(x) ischar(x);
val_logical     = @(x) islogical(x);

% Add parameters with their default values that will only be used if no
% value is specified when calling the function (values of the parameters
% are provided in varargin otherwise)
addParamValue(p, 'order',                default_order,                val_numeric); 
addParamValue(p, 'vehicle_info',         default_vehicle_info);

addParamValue(p, 'make_figures', default_make_figures, val_logical);
addParamValue(p, 'plot_range',   default_plot_range,   val_logical);
addParamValue(p, 'plot_scaled',  default_plot_scaled,  val_logical);
addParamValue(p, 'plot_xcorr',   default_plot_xcorr,   val_logical);
addParamValue(p, 'plot_regress', default_plot_regress, val_logical);
addParamValue(p, 'plot_result',  default_plot_result,  val_logical);

% -- Parse the arguments.
parse(p, varargin{:});

% -- Define simple variables to use throughout function
order                = p.Results.order;
vehicle_info         = p.Results.vehicle_info;
make_figures         = p.Results.make_figures;
plot_range           = p.Results.plot_range;
plot_scaled          = p.Results.plot_scaled;
plot_xcorr           = p.Results.plot_xcorr;
plot_regress         = p.Results.plot_regress;
plot_result          = p.Results.plot_result;
% end of input argument checking.

fig_num = 0;

fprintf('Temperature regression: order %0d\n',order)

% -----------------------------------------------------------------------
% --- Load the data -----------------------------------------------------
% -----------------------------------------------------------------------

%%%%
% We need the name of the FP07 signals with and without pre-emphasis. The
% name without pre-emphasis is usually the section in the setup.cfg-file
% that contains the processing parameters. However, sometimes this signal
% does not exist and the information has to be gathered from the section
% for the signal with pre-emphasis. For eample, with a Sea-glider.

T_with_pre_emphasis_string = [T_string '_d' T_string]; % Should be T1_dT1 or T2_dT2
T_without_pre_emphasis_string = T_string; % Should be T1 or T2

% -----------------------------------------------------------------------
% --- Get relevant temperature data (and deconvolve if necessary) -------
% -----------------------------------------------------------------------

% We must identify the name of the section in the setup.cfg-file that
% contains the processing paramters.
section_name = T_without_pre_emphasis_string;
eval(['T_without_pre_emphasis = DATA.' T_without_pre_emphasis_string ';'])

% Assign temperature vectors
eval (['T_ref = DATA.' T_ref_string ';']) % T_ref is the reference thermometer, usually SBT or JAC_T
eval(['T = DATA.' T_without_pre_emphasis_string ';']) % T is the thermistor signal without pre-emphasis

% If we have a signal with pre-emphasis, then we will use it to form the signal T.
if ~isempty(T_with_pre_emphasis_string)
    eval(['T = DATA.' T_with_pre_emphasis_string ';']) % T is the thermistor signal with pre-emphasis
    T = deconvolve(...
        T_with_pre_emphasis_string, ...
        T_without_pre_emphasis, ...
        T, ...
        DATA.fs_fast, DATA.setupfilestr);
    
    % sampling rate ratio (used to downsample)
    ratio = round(DATA.fs_fast / DATA.fs_slow);
    
    % down size to match T_ref
    T = reshape(T, ratio, []); 
    T = mean(T)';
end

%----------------------------------
% - Warning if range is too small - 
%----------------------------------
if max(T_ref(m))-min(T_ref(m))<=8 && order>1
    warning(['Temperature range is less than 8 degrees '...
        '-> Recommend using FIRST-ORDER calibration.',...
        ' Exit function (typically by CRTL+C) and re-run with additional inputs,'...
        ' - i.e. ''order'',1 - '...
        ' or modify cal_info structure.',...
        ' Then, in setup file, delete entire ''beta_2=__'' line and ',...
        ' update T_0 and beta_1 values before patching.']) %, pause
end

%-------------------------------------------------------------------------
% -- Low-Pass filter thermistor data and compute cross correlation -------
%-------------------------------------------------------------------------
% Low-pass filtering the thermistor data to make it more compatible with the JAC-T
W_mean = abs(mean(DATA.W_slow(m)));
fc = 0.73 * sqrt(W_mean / 0.62); % in Hz
[b,a] = butter(1, fc / (DATA.fs_slow/2));
T = filter(b, a, T);

% Compute cross-correlation
max_lag = round(10*DATA.fs_slow); % estimate of the max lag required to find the actual lag.
[bb, aa] = butter(2,4/(DATA.fs_slow/2)); % 4 Hz smoother to suppress high-frequency noise
[correlation, lags] = xcorr(...
    filter(bb,aa,detrend(diff(T(m)))),...
    filter(bb,aa,detrend(diff(T_ref(m)))),max_lag,'coeff');
[max_corr, m_lag] = max(abs(correlation));
junk_m = m_lag; % needed for figure
m_lag = m_lag - max_lag - 1;
Lag    = m_lag / DATA.fs_slow; % in seconds and should be negative.

%-----------------------------------------------------
% -- Do regression to get thermistor coefficients ----
%           (using Steinhart-Hart equation)
%-----------------------------------------------------

% Copy only the profile data
T_ref_prof = T_ref(m);
T_prof     = T(m);
P_prof     = DATA.P_slow(m);

% First align the T and T_ref signals using m_lag.
if m_lag >0, m_lag = 0; end % m_lag is expected to be negative because 
                            % reference sensor is physically 'behind' probes
P_prof     = P_prof(1:end+m_lag); 
T_prof     = T_prof(1:end+m_lag);
T_ref_prof = T_ref_prof(1-m_lag:end); % shift reference temperature
T_ref_regress = T_ref_prof + 273.15; % in kelvin
T_ref_regress = 1 ./ T_ref_regress;

% Now gather information about the electronics for this thermistor.
my_object = setupstr( DATA.setupfilestr );
therm_type =    (char(setupstr( my_object, section_name, 'type')));
E_B = str2double(char(setupstr( my_object, section_name, 'E_B')));
a   = str2double(char(setupstr( my_object, section_name, 'a'  )));
b   = str2double(char(setupstr( my_object, section_name, 'b'  )));
G   = str2double(char(setupstr( my_object, section_name, 'G'  )));
adc_fs   = str2double(char(setupstr( my_object, section_name, 'adc_fs'  )));
adc_bits = str2double(char(setupstr( my_object, section_name, 'adc_bits'  )));
try zero = str2double(char(setupstr( my_object, section_name, 'adc_zero'  )));catch, zero = 0; end

% Compute non-dimensional thermistor voltage
if strcmp(therm_type, 'therm')
    factor = (adc_fs / 2^adc_bits)*2 / (G*E_B);
    Z = factor*(T_prof - a)/b;
elseif strcmp(therm_type, 't_ms')
    Z = T_prof * (adc_fs/2^adc_bits) + zero;
    Z = ((Z - a)/b) *2 / (G*E_B);
end

% Compute resistance ratio for this thermistor.
RT_R0 = (1 - Z) ./ (1 + Z); 
RT_R0 = log(RT_R0);

% Generate the coefficients for this thermistor.
beta = zeros(1,order);
p = polyfit(RT_R0, T_ref_regress, order); 
pp = p; % save for later usage
p = 1 ./ p;
p = fliplr(p); % place in ascending order
T_0    = p(1);
for index = 2:order+1
    beta(index-1) = p(index);
end

% make a smooth line using coefficients (for comparison)
R = linspace(min(RT_R0), max(RT_R0), 1000);
R = R';
T_inverse_predicted = polyval(pp,R);

%------------------------------------------------------------------
% -- Use the computed co-efficients to estimate temperature -------
%------------------------------------------------------------------
% Estimate temperature values using computed calibration coeffients
T_calibrated = polyval(pp, RT_R0);
T_calibrated = 1 ./ T_calibrated;
T_calibrated = T_calibrated - 273.15;


%% EPFL change
%------------------------------------------------------------------
% -- Update of the cfg file ---------------------------------------
%------------------------------------------------------------------
A = regexp(fileread([cfgfile '.cfg']),'\n','split');
idx_beta1 = find(contains(A,'beta_1'));
idx_beta2 = find(contains(A,'beta_2'));
if isempty(idx_beta2)
    if isempty(idx_beta1)
        idx_beta1 = find(contains(A,'beta')); % Only one coefficient specified
        for k=1:length(idx_beta1) % Add beta_1 field
            A{idx_beta1(k)}='beta_1      = ';
        end
    end
    for k=1:length(idx_beta1) % Add beta_2 field
        A=[A(1:idx_beta1(k)+k-1),{'beta_2      = '},A(idx_beta1(k)+k:end)];
    end
    idx_beta2 = find(contains(A,'beta_2'));
end
idx_T = find(endsWith(A,['= ' T_string]) | contains(A,['= ' [T_string '_']]));
add_space=true;
if isempty(idx_T)
    %idx_T = find(endsWith(A,['=' T_string]) | contains(A,['=' [T_string '_']])); % No space
    idx_T = find(contains(A,['=' T_string])); % No space
    add_space=false;
end
idx_T0 = find(contains(A,'T_0'));

idx_beta1 = idx_beta1( idx_beta1 > idx_T(1) & idx_beta1 < idx_T(2)); % Beta1 coefficient of the specific T probe
idx_beta2 = idx_beta2( idx_beta2 > idx_T(1) & idx_beta2 < idx_T(2)); % Beta2 coefficient of the specific T probe
idx_T0 = idx_T0( idx_T0 > idx_T(1) & idx_T0 < idx_T(2)); % T0 coefficient of the specific T probe

if add_space
    A(idx_beta1) =  {['beta_1      = ', num2str(beta(1),'%15.6f')]};
    if order==2
        A(idx_beta2) =  {['beta_2      = ', num2str(beta(2),'%15.6f')]};
        A(idx_T0) =  {['T_0      = ', num2str(T_0,'%15.6f')]};
    else
        A(idx_T0) =  {['T_0      = ', num2str(T_0,'%15.6f')]};
        A(idx_beta2)=[]; % No beta_2
    end
else
    A(idx_beta1) =  {['beta_1=', num2str(beta(1),'%15.6f')]};
    if order==2
        A(idx_beta2) =  {['beta_2=', num2str(beta(2),'%15.6f')]};
        A(idx_T0) =  {['T_0=', num2str(T_0,'%15.6f')]};
    else
        A(idx_T0) =  {['T_0=', num2str(T_0,'%15.6f')]};
        A(idx_beta2)=[]; % No beta_2
    end

end

fid = fopen([cfgfile '.cfg'], 'w');
fprintf(fid,'%s\n',A{:});
fclose(fid);





