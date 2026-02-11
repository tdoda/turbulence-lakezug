function [BIN,SLOW,FAST] = resolve_turbulence(DATA, kprof, param,folder_out,profID,make_plot_prof,ind_plot_spectra,run_dissip,show_progress)
%RESOLVE_TURBULENCE Compute turbulence estimates for a given profile based on the function
%"resolve_microCTD_profile_all".
%
%   INPUTS:
%   DATA (structure): profiling data from a given data file (output of
%   odas_p2mat).
%   prof (int): profile index between 1 and Nprf.
%   param (structure): parameters specific to the campaign.
%   folder_out (string): path to the output folder.
%   profID (int): profile counter.
%   make_plot_prof (boolean) [optional]: = True to make profiles plots. Default value:
%   False.
%   ind_plot_spectra (int array) [optional]: indices of bins to make spectrum plot. Default value:
%   0.
%   run_dissip (boolean) [optional]: = True if turbulence analysis should
%   be done. Default value: True.
%   show_progress (boolean) [optional]: = True to show progress bar and warning messages. 
%   Default value: True.
%
%
%   OUTPUTS:
%   BIN (structure): turbulence estimates and related quantities for each bin.
%   SLOW (structure): low-frequency quantities.
%   FAST (structure): high-frequency quantities.
%
% T. Doda based on S. Piccolroaz, last version: 22.01.2026
%% Set up parameters
if nargin<6
    make_plot_prof=false;
end

if nargin<7
    ind_plot_spectra=0;
end

if nargin<8
    run_dissip=true;
end

if nargin<9
    show_progress=true;
end

info=param.info;


%defines times
time_fast0 = (0:1:length(DATA.P_fast)-1)/DATA.fs_fast; % [s]
time_slow0= (0:1:length(DATA.P_slow)-1)/DATA.fs_slow; % [s]

ind_prof_slow = DATA.ind_prof_slow;
ind_prof_fast = DATA.ind_prof_fast;
NP=DATA.Nprf;

if kprof<1 || kprof>NP
    fprintf('Wrong number of profiles')
    return
end

%gets index for the desired profile
iipf = ind_prof_fast(1,kprof):ind_prof_fast(2,kprof);
iips = ind_prof_slow(1,kprof):ind_prof_slow(2,kprof);

%date
% date = datenum(DATA.Year, DATA.Month, DATA.Day, DATA.Hour, DATA.Minute, DATA.Second );
% date = date + ind_prof_slow(1,kprof)/DATA.fs_slow/60/60/24;
% datestr(date)

%sampling frequencies
fss = DATA.fs_slow;
fsf = DATA.fs_fast;

%gets fast response sensors
timef = time_fast0(iipf);
Pf = DATA.P_fast(iipf);
if param.config.T1
    T1f = DATA.T1_fast(iipf);
    T1_dT1 = DATA.T1_dT1(iipf);
    gradT1f = DATA.gradT1(iipf);
else
    T1f = zeros(size(Pf));
    T1_dT1 = zeros(size(Pf));
    gradT1f = zeros(size(Pf));
end

if param.config.T2
    T2f = DATA.T2_fast(iipf);
    T2_dT2 = DATA.T2_dT2(iipf);
    gradT2f = DATA.gradT2(iipf);
else
    T2f = zeros(size(Pf));
    T2_dT2 = zeros(size(Pf));
    gradT2f = zeros(size(Pf));
end

if param.config.S1
    sh1 = DATA.sh1(iipf);
else
    sh1 = zeros(size(Pf));
end

if param.config.S2
    sh2 = DATA.sh2(iipf);
else
    sh2 = zeros(size(Pf));
end

Ax = DATA.Ax(iipf);
Ay = DATA.Ay(iipf);
Wf = DATA.W_fast(iipf);
AA = [Ax,Ay];


%gets slow response sensors
times = time_slow0(iips);
Ws = DATA.W_slow(iips);
Ps = DATA.P_slow(iips);
T_JAC = DATA.(param.CTD_T)(iips);
C_JAC = DATA.(param.CTD_C)(iips);
C_corr=DATA.Cond_corr(iips);

if strcmp(param.instrument,'microCTD') % MicroCTD
    Incl_X = DATA.Incl_X(iips);
    Incl_Y = DATA.Incl_Y(iips);
elseif strcmp(param.instrument,'VMP') % VMP: no inclination data
    Incl_X = NaN(size(T_JAC));
    Incl_Y = NaN(size(T_JAC));
else
    error("Wrong instrument name")
end

%salinity and density
sgt=DATA.rhoTS(iips);
Ss=DATA.Sal(iips);
depth=DATA.depth(iips);

%accurate T in fast sensors grid
T_JAC_fast = interp1(Ps,T_JAC,Pf);
i1 = find(isfinite(T_JAC_fast),1,'first');
if i1>1
    T_JAC_fast(1:i1-1) = T_JAC_fast(i1); % Repeat first value
end
i2 = find(isfinite(T_JAC_fast),1,'last');
if i2<length(T_JAC_fast)
    T_JAC_fast(i2+1:end) = T_JAC_fast(i2); % Repeat last value
end

% Chl and Turb
if strcmp(param.instrument,'microCTD') % MicroCTD, fast sensors
    Chl_fast = DATA.(param.CTD_Chl)(iipf);
    Turb_fast = DATA.(param.CTD_Turb)(iipf);
    % Convert to slow rate:
    Chl_slow=interp1(Pf,Chl_fast,Ps);
    Turb_slow=interp1(Pf,Turb_fast,Ps);
elseif strcmp(param.instrument,'VMP') % VMP, slow sensors
    % Chl = zeros(size(Ax));
    % Turb = zeros(size(Ax));
    Chl_slow = DATA.(param.CTD_Chl)(iips);
    Turb_slow = DATA.(param.CTD_Turb)(iips);
    % Convert to fast rate:
    Chl_fast=interp1(Ps,Chl_slow,Pf);
    Turb_fast=interp1(Ps,Turb_slow,Pf);
else
    error("Wrong instrument name")
end

%sets maximum depth
% if ~ isfield(info,'pmax')
%     info.pmax = round(max(Ps)-info.dpD/2)+1;
% end


%% Create FAST structure
FAST.datenum = datenum(DATA.Year, DATA.Month, DATA.Day, DATA.Hour, DATA.Minute, DATA.Second);
FAST.filename = DATA.filename;
if ismember(info.prof_dir,'up')
    FAST.upward = 1;
else
    FAST.upward = 0;
end
%FAST.dtsec=timef;
FAST.dtsec=DATA.t_fast(iipf); % s, time from odas_p2mat()
FAST.pressure = Pf;
FAST.raw_pressure=DATA.P_fast_raw(iipf); 
FAST.velocity=Wf;
FAST.fast_T1 = T1f;
FAST.grad_T1 = gradT1f;
FAST.fast_T2 = T2f;
FAST.grad_T2 = gradT2f;
FAST.chlorophyll = Chl_fast;
FAST.turbidity = Turb_fast;
FAST.fast_S1 = sh1;
FAST.fast_S2 = sh2;
FAST.A_x = Ax;
FAST.A_y = Ay;
FAST.depth = interp1(Ps,depth,Pf,'linear','extrap');
FAST.sigmatT1 = DATA.rhoT1S(iipf);
FAST.sigmatT2 = DATA.rhoT2S(iipf);

%% Create SLOW structure
SLOW.datenum = FAST.datenum;
SLOW.filename = DATA.filename;
if ismember(info.prof_dir,'up')
    SLOW.upward = 1;
else
    SLOW.upward = 0;
end
%SLOW.dtsec=times; % s
SLOW.dtsec=DATA.t_slow(iips); % s, time from odas_p2mat()
SLOW.pressure = Ps;
SLOW.raw_pressure=DATA.P_slow_raw(iips);
SLOW.velocity=Ws;
SLOW.depth = depth;
SLOW.temperature = T_JAC;
SLOW.conductivity = C_corr; %uS/cm
SLOW.conductivity20 = DATA.Cond_20(iips); %uS/cm
SLOW.salinity = 1000*Ss; %mg/l
SLOW.density = sgt;
SLOW.chlorophyll = Chl_slow;
SLOW.turbidity = Turb_slow;
SLOW.Incl_x = Incl_X;
SLOW.Incl_y = Incl_Y;

%% Calculates displacements for thorpe length
disp('>>>>> Calculation of Thorpe displacements')
if ismember(info.prof_dir,'down')
    [sort_rho, isd] = sort(sgt); % Sort slow density profile, from low to high density
    displ = Ps - Ps(isd); % displacements [m] = distance between the depth i and the depth where the sorted density is measured
    % Depth of sorted density value = current depth i - displacement i

    [sort_uT1, ist1] = sort(T1f,'descend'); % Sort fast T1 profile, from high to low temperature
    displuT1 = Pf - Pf(ist1);

    [sort_uT2, ist2] = sort(T2f,'descend'); % Sort fast T2 profile, from high to low temperature
    displuT2 = Pf - Pf(ist2);

    [sort_T_JAC, ~] = sort(T_JAC,'descend'); % Sort JAC_T profile (not for Thorpe displacements), from high to low temperature
else
    [sort_rho, isd] = sort(sgt,'descend'); % Sort slow density profile, from high to low density
    displ = Ps - Ps(isd);

    [sort_uT1, ist1] = sort(T1f,'ascend'); % Sort fast T1 profile, from low to high temperature
    displuT1 = Pf - Pf(ist1);

    [sort_uT2, ist2] = sort(T2f,'ascend');% Sort fast T2 profile, from low to high temperature
    displuT2 = Pf - Pf(ist2);

    [sort_T_JAC, ~] = sort(T_JAC,'ascend'); % Sort slow JAC_T profile (not for Thorpe displacements), from low to high temperature
end

% centered Thorpe scale
Ps_centered=Ps-displ*0.5; % Centred pressure values = pressure of the middle of the Thorpe segments (~center of overturns)
Lc=zeros(1,length(Ps)); count_c=zeros(1,length(Ps));
for i_Ps=1:length(Ps)
    [~,i_c] = min(abs(Ps-Ps_centered(i_Ps))); % Closest pressure value to the centered pressure values
    Lc(i_c) = Lc(i_c) + abs(displ(i_Ps));
    count_c(i_c) = count_c(i_c) + 1;
end
Lc = Lc./count_c; Lc(count_c==0)=0; % Thorpe scale profile computed as averaged Thorpe displacement for depths with several overturns
SLOW.LTc = Lc';
SLOW.LTavg=sqrt(mean(displ.^2,"omitnan")); % Added usual Thorpe scale definition, one value per profile (T. Doda)
% Note: LTavg is similar to the averaged non zero LTc


%% Binned temperature, salinity and density
disp('>>>>> Create bins')

%defines the pressure vector where to calculate
pres = info.pmin:info.dp:info.pmax; % Pressure of the bins center (space between bins: overlapping dp)
if ismember(info.prof_dir,'up') % For upward profiles, start the profile after the acceleration
    idx75=find(Ps<0.75*info.pmax,1,'first'); % Index of the depth corresponding to 25% of the entire profile (75 % from the surface). Note: Ps is decreasing (upward profile).
    [~,idx]=min(Ws(1:idx75)); % Get the minimum speed over these first 25% (not affected by wave-induced oscillations above)
    n_pres = find(pres<=Ps(idx),1,'last'); % Maximal depth where the profile starts (for upward)
else
    n_pres=length(pres); % Use pmax as the maximal depth for downward profile
end


BIN.datenum = FAST.datenum;
BIN.date=cellstr(datestr(BIN.datenum,'yyyy-mm-dd'));
BIN.time=cellstr(datestr(BIN.datenum,'HH:MM:SS'));
BIN.filename = DATA.filename;
if ismember(info.prof_dir,'up')
    BIN.upward = 1;
else
    BIN.upward = 0;
end
BIN.filename=cellstr(BIN.filename);
if BIN.upward==1
    BIN.direction={'upward'};
else
    BIN.direction={'downward'};
end

BIN.pressure = pres; % Pressure of each bin [dbar]
% BIN.depth = pres_av(Ps,depth,pres,info.dp,2.7); % Average depth for each bin [m]
% BIN.temperature = pres_av(Ps,T_JAC,pres,info.dp,2.7); % Average temperature for each bin [°C]
% BIN.conductivity = 1000*pres_av(Ps,C_JAC,pres,info.dp,2.7); % Average condcutivity for each bin [uS/cm]
% BIN.salinity = 1000*pres_av(Ps,Ss,pres,info.dp,2.7); % Average salinity for each bin [mg/l]
% BIN.density = pres_av(Ps,sgt,pres,info.dp,2.7); % Average density for each bin [kg/m3]

% Use dpD instead of dp for bin thickness:
BIN.depth = pres_av(Ps,depth,pres,info.dpD,2.7); % Average depth for each bin [m]
BIN.temperature = pres_av(Ps,T_JAC,pres,info.dpD,2.7); % Average temperature for each bin [°C]
BIN.conductivity = 1000*pres_av(Ps,C_JAC,pres,info.dpD,2.7); % Average condcutivity for each bin [uS/cm]
BIN.salinity = 1000*pres_av(Ps,Ss,pres,info.dpD,2.7); % Average salinity for each bin [mg/l]
BIN.density = pres_av(Ps,sgt,pres,info.dpD,2.7); % Average density for each bin [kg/m3]

Thorpe=sqrt(pres_av(Ps,displ.^2,pres,info.dpD,0)); % Thorpe length for each bin (average absolute displacement) from slow density profile [m]
BIN.LT = Thorpe;

for i=1:length(pres) % Loop across the bins
    % It is better to use FP07 because when the lake is homogeneous,
    % the resolution of JAC-T is too low, resulting in grT=1E-16
    BIN.avggradT1(i) = mean_grad(Pf,sort_uT1,pres(i),max(Thorpe(i),info.dpD));   % SEB: evaluated on the resorted T profile, over the largest length scale between bin size and Thorpe length
    BIN.avggradT2(i) = mean_grad(Pf,sort_uT2,pres(i),max(Thorpe(i),info.dpD));   % SEB: evaluated on the resorted T profile, over the largest length scale between bin size and Thorpe length
    if param.config.T1 && param.config.T2
        BIN.grT12(i)=nanmean([BIN.avggradT1(i),  BIN.avggradT2(i)]); % Mean gradient between the two FP07
    elseif ~param.config.T2 % Only T1
        BIN.grT12(i)=BIN.avggradT1(i);
    elseif ~param.config.T1 % Only T2
        BIN.grT12(i)=BIN.avggradT2(i);
    else % No FP07
        BIN.grT12(i)=NaN(size(BIN.avggradT1(i)));
    end
    BIN.avggradT(i) = mean_grad(Ps,sort_T_JAC,pres(i),max(Thorpe(i),info.dpD));   % SEB: evaluated on the resorted T profile, over the largest length scale between bin size and Thorpe length
    BIN.N2(i) = -9.81*mean_grad(Ps,sort_rho,pres(i),max(Thorpe(i),info.dpD))/1000;
end

% BIN.chlorophyll = pres_av(Pf,Chl,pres,info.dp,2.7);
% BIN.turbidity = pres_av(Pf,Turb,pres,info.dp,2.7);

% Use dpD instead of dp for bin thickness:
BIN.chlorophyll = pres_av(Pf,Chl_fast,pres,info.dpD,2.7);
BIN.turbidity = pres_av(Pf,Turb_fast,pres,info.dpD,2.7);

BIN.LTuT1 = sqrt(pres_av(Pf,displuT1.^2,pres,info.dpD,0)); % Thorpe length scale from T1
BIN.LTuT2 = sqrt(pres_av(Pf,displuT2.^2,pres,info.dpD,0)); % Thorpe length scale from T2
BIN.acc = -mean_grad_time(Pf,Wf,pres,info.dpD,timef);  % instrument acceleration [dbar.s^(-2)]

pmaxplot=pres(n_pres); %pmaxplot is only used for the plot


%% Filtering
disp('>>>>> Filtering')
%filters shear and highpasses microstructure
% use the odas function "despike" to remove short-duration spikes
[sh1, ~, ~, ~ ] =  despike(sh1, info.despike_sh(1), info.despike_sh(2), fsf, round(info.despike_sh(3)*fsf));
if param.config.S2
    [sh2, ~, ~, ~ ] =  despike(sh2, info.despike_sh(1), info.despike_sh(2), fsf, round(info.despike_sh(3)*fsf));
end



mW = (max(Pf)-min(Pf))/(max(timef)-min(timef)); % Mean speed
if info.k_HP_cut>0
    f_HP_cut = info.k_HP_cut*mW; % cutoff frequency [s-1]
    [bh,ah] = butter(1, f_HP_cut/(fsf/2), 'high'); % Defines a 1st order Butterworth highpass filter with cutoff frequency normalized by half of the sampling frequency (Matlab function)
    if param.config.S1
        sh1_hp = filter(bh, ah, sh1); % Filters the data to enhance frequencies above f_HP_cut and remove low-frequency noise (Matlab function)
        sh1_hp = flipud(sh1_hp);
        sh1_hp = filter(bh, ah, sh1_hp);
        sh1_hp = flipud(sh1_hp);
    else
        sh1_hp=sh1;
    end

    if param.config.S2
        sh2_hp = filter(bh, ah, sh2);
        sh2_hp = flipud(sh2_hp);
        sh2_hp = filter(bh, ah, sh2_hp);
        sh2_hp = flipud(sh2_hp);
    else
        sh2_hp=sh2;
    end

else
    sh1_hp = sh1;
    sh2_hp = sh2;

end

FAST.fast_S1_hp=sh1_hp;
FAST.fast_S2_hp=sh2_hp;

% No filtering of the temperature data
T1f_hp=T1f;
T2f_hp=T2f;

% Despike the piezo-accelerometer signals
piezo_accel_num = size(AA,2);
if  ~isempty(AA) && info.despike_A(1) ~= inf
    for probe = 1:piezo_accel_num % For each accelerometer direction
        [AA(:,probe), ~, ~, ~]  = ...
            despike(AA(:,probe),  info.despike_A(1), ...
            info.despike_A(2), fsf, round(info.despike_A(3)*fsf));
    end
end

if ~param.config.S2 % Only one shear: select one of the two vibration signals

    %-- identify the piezo-accelerometer signal to be used for noise filtering
    % Compute autospectra of the detrended signals of sh1, Ax and Ay (i.e., residuals with
    % respect to linear fit)
    PSDsh = csd_odas(detrend(sh1_hp),detrend(sh1_hp),1024,fsf,[],512,'linear');
    PSDA1 = csd_odas(detrend(AA(:,1)),detrend(AA(:,1)),1024,fsf,[],512,'linear');
    PSDA2 = csd_odas(detrend(AA(:,2)),detrend(AA(:,2)),1024,fsf,[],512,'linear');
    % Compute cross-spectra sh1-Ax and sh1-Ay:
    CSDshA1 = csd_odas(detrend(sh1_hp),detrend(AA(:,1)),1024,fsf,[],512,'linear');
    [CSDshA2,fA] = csd_odas(detrend(sh1_hp),detrend(AA(:,2)),1024,fsf,[],512,'linear');
    % Compute magnitude-squared coherence ranging from 0 to 1 and measuring how strongly 
    % accelerometer and shear signals are linearly related at each frequency
    COHshA1 = abs(CSDshA1.^2./(PSDA1.*PSDsh));
    COHshA2 = abs(CSDshA2.^2./(PSDA2.*PSDsh));
    AAxy=AA; % store both signals
    if nanmean(COHshA1)>nanmean(COHshA2)
        AA = AA(:,1);
    else
        AA = AA(:,2);
    end

end

FAST.AA_filt = AA;

%% Plot

if make_plot_prof
    plot_profile(DATA,SLOW,FAST,param,pmaxplot)

    if ~exist([folder_out,'Figures'],'dir')
        mkdir([folder_out,'Figures'])
    end
    saveas(gcf,[folder_out,'Figures/profile',num2str(kprof,'%02d'),'.fig']);
    exportgraphics(gcf,[folder_out,'Figures/profile',num2str(kprof,'%02d'),'.png'],'Resolution',400);
end

%% Defines output variables in BIN structure

BIN.ThermCond= nan(1,length(pres));
BIN.Diff= nan(1,length(pres));
BIN.DynVisco= nan(1,length(pres));
BIN.KinVisco= nan(1,length(pres));

BIN.speed = nan(1,length(pres));
BIN.WW= nan(1,length(pres));
BIN.maxgradT = nan(1,length(pres));
BIN.eps_S1 = nan(1,length(pres));
BIN.kB_S1 = nan(1,length(pres));
BIN.MAD_S1 = nan(1,length(pres));
BIN.MADc = nan(1,length(pres));
BIN.flag_S1 = nan(1,length(pres));
BIN.kL_S1 = nan(1,length(pres));
BIN.kU_S1 = nan(1,length(pres));

%if param.config.S2
BIN.eps_S2 = nan(1,length(pres));
BIN.kB_S2 = nan(1,length(pres));
BIN.MAD_S2 = nan(1,length(pres));
BIN.flag_S2 = nan(1,length(pres));
BIN.kL_S2 = nan(1,length(pres));
BIN.kU_S2 = nan(1,length(pres));
%end

BIN.Xi_ST1 = nan(1,length(pres));
BIN.Xi_T1 = nan(1,length(pres));
BIN.kB_T1 = nan(1,length(pres));
%BIN.sXif1 = nan(1,length(pres));
%BIN.sKB1 = nan(1,length(pres));
BIN.Xiv1 = nan(1,length(pres));
BIN.kU_T1 = nan(1,length(pres));
BIN.eps_T1 = nan(1,length(pres));
BIN.epsT1max = nan(1,length(pres));
BIN.MAD_ST1 = nan(1,length(pres));
BIN.MAD_T1 = nan(1,length(pres));
%BIN.LKH1 = nan(1,length(pres));
BIN.LR_T1 = nan(1,length(pres));
BIN.krange_T1 = nan(1,length(pres));
BIN.flag_T1 = nan(1,length(pres));
BIN.kpeak_T1 = nan(1,length(pres));
BIN.kL_T1 = nan(1,length(pres));

BIN.Xi_ST2 = nan(1,length(pres));
BIN.Xi_T2 = nan(1,length(pres));
BIN.kB_T2 = nan(1,length(pres));
BIN.sXif2 = nan(1,length(pres));
BIN.sKB2 = nan(1,length(pres));
BIN.Xiv2 = nan(1,length(pres));
BIN.kU_T2 = nan(1,length(pres));
BIN.eps_T2 = nan(1,length(pres));
BIN.epsT2max = nan(1,length(pres));
BIN.MAD_ST2 = nan(1,length(pres));
BIN.MAD_T2 = nan(1,length(pres));
BIN.LKH2 = nan(1,length(pres));
BIN.LR_T2 = nan(1,length(pres));
BIN.krange_T2 = nan(1,length(pres));
BIN.flag_T2 = nan(1,length(pres));
BIN.kpeak_T2 = nan(1,length(pres));
BIN.kL_T2 = nan(1,length(pres));
BIN.flag_vibration = nan(1,length(pres));
BIN.flag_inclination = nan(1,length(pres));
BIN.flag_acceleration = nan(1,length(pres));
BIN.flag_ST1 = nan(1,length(pres));
BIN.flag_ST2 = nan(1,length(pres));

disp('>>>>> Compute turbulence quantities from shear and FP07')
for i = 1:n_pres %length(pres)
    if show_progress
        fprintf('>>>>>>> Progress: %d/%d (%0.2f %%)\n',i,n_pres,i/n_pres*100)
    end
    jp = find(Pf>=pres(i)-info.dpD/2 & Pf<=pres(i)+info.dpD/2); % cell indices that are part of bin i in fast channels
    jps = find(Ps>=pres(i)-info.dpD/2 & Ps<=pres(i)+info.dpD/2);   % cell indices that are part of bin i in slow channels

    % SEB: added evaluation of "maximum" gradT for pyroelectric analysis
    if ~isempty(jps) & length(jps)>=3 % At least 3 slow cells in bin i
        tmp = movavg(T_JAC(jps),'linear',3); % Smoothed temperature data
        BIN.maxgradT(i) = max(-gradient(tmp)./gradient(Ps(jps))); % Maximal smoothed gradient
    end

    % SEB: added evaluation of D and viscosity (based on http://web.mit.edu/seawater/)
    avT = nanmean(T_JAC(jps)); % Average temperature in bin i
    avS = nanmean(Ss(jps)); % Average salinity
    avsigmat = nanmean(sgt(jps)); % Average density
    cond = SW_Conductivity(avT,'C',avS,'ppt'); % Thermal conductivity
    BIN.ThermCond(i)=cond;
    cp = SW_SpcHeat(avT,'C',avS,'ppt',0.101325,'MPa');   % The last argument is sat vapor pressure. It can be kept constant as suggested in SW_Diffusivity.m
    D = cond/(avsigmat*cp); % Thermal diffusivity
    BIN.Diff(i)=D;

    dyn_visco = SW_Viscosity(avT,'C',avS,'ppt'); % dynamic viscosity mu
    BIN.DynVisco(i)=dyn_visco;
    kin_visco = dyn_visco/avsigmat; % kinematic viscosity nu
    BIN.KinVisco(i)=kin_visco;

    % Speed (added by T. Doda)
    BIN.speed(i)=nanmean(DATA.speed_slow(iips(jps)));

    if isfield(info,'num_fft')
        Nfft = floor(length(jp)/((info.num_fft+1)/2));
    else
        Nfft = floor(length(jp)/2);
    end

    if isfield(info,'overlap_pct') % If overlap percentage for spectral analysis is specified
        overlap = round(info.overlap_pct/100*Nfft);
    else
        overlap = round(Nfft/2); % Default overlap for spectral analysis: 50 %
    end

    if length(jp)<Nfft | length(jp)<128 % Not enough values in the bin
        continue
    end

    % SEB: calculate the profiling speed to be used below. I prefer this way
    % than using Wf or Ws, since these are smoothed and for short bins
    % may be critical. However, in case just define WW=mean(abs(Wf(jp)))
    % Note that in the main script I removed the smoothing, however I
    % leave this definition to be safe.
    WW = mean( abs( (Pf(jp(end)) - Pf(jp(1)))/(timef(jp(end)) - timef(jp(1))) ) );
    BIN.WW(i)=WW;
    %% acceleration, vibration and inclination flags
    % flagging for too much vibrations, flag if more than 5% of
    % vibration per bin is >200 or if standard deviation is >100
    % (i.e. 2 sigma C.L. = 95%)
    Axi = Ax(jp);
    Ayi = Ay(jp);
    Axidev = nanstd(Axi);
    Ayidev = nanstd(Ayi);
    Aximean = abs(nanmean(Axi));
    Ayimean = abs(nanmean(Ayi));
    len_Axi = nansum(abs(Axi)>200); % Number of values with high Ax
    len_Ayi = nansum(abs(Ayi)>200); % Number of values with high Ay
    if len_Axi/length(jp) > 0.05
        BIN.flag_vibration(i) = 1;
    elseif (Axidev + Aximean) > 100
        BIN.flag_vibration(i) = 1;
    elseif len_Ayi/length(jp) > 0.05
        BIN.flag_vibration(i) = 1;
    elseif (Ayidev + Ayimean) > 100
        BIN.flag_vibration(i) = 1;
    else
        BIN.flag_vibration(i) = 0;
    end

    % flagging for inclination, same as vibr for incl. > 5°
    Incli = Incl_X(jps);
    Inclidev = nanstd(Incli);
    Inclimean = abs(nanmean(Incli));
    len_Incli = sum(abs(Incli)>5);
    if len_Incli/length(jps) > 0.05
        BIN.flag_inclination(i) = 1;
    elseif (Inclidev + Inclimean) > 2.5
        BIN.flag_inclination(i) = 1;
    else
        BIN.flag_inclination(i) = 0;
    end

    %flagging for acceleration of instrument
    %flag if acceleration is > 0.01 m/s/s (preliminary value, specific
    %to lake Garda downward profiles, conservative estimate)
    if strcmp(info.system,'Gar') & ismember(info.prof_dir,'down')
        if abs(BIN.acc(i)) > 0.01
            BIN.flag_acceleration(i) = 1;
        else
            BIN.flag_acceleration(i) = 0;
        end
    else
        BIN.flag_acceleration(i) = NaN;
    end

    
    %% Shear spectral calculations
    close all % Close all the figures
    if ismember(i,ind_plot_spectra)
        make_plot_spectra=true;
    else
        make_plot_spectra=false;
    end
    if  run_dissip
        if strcmp(info.Nasmyth_spec,'EPFL')
            if param.config.S1
                if param.config.S2
                    try
                        % [BIN.eps_S2(i), BIN.MAD_S2(i), ~,BIN.flag_S2(i),  BIN.kL_S2(i),BIN.kU_S2(i)] = ...
                        %     TKE_dis_spec(Pf(jp),[sh1_hp(jp) sh2_hp(jp)],AA(jp,:),0.1,14,info.fAA,kin_visco,WW, Nfft, overlap, info.noise_corr,'sh2',make_plot_spectra,Pf(jp),T1f(jp),folder_out,filename,profID);
                        % [BIN.eps_S1(i), BIN.MAD_S1(i), BIN.MADc(i),BIN.flag_S1(i), BIN.kL_S1(i),BIN.kU_S1(i)] = ...
                        %     TKE_dis_spec(Pf(jp),[sh1_hp(jp) sh2_hp(jp)],AA(jp,:),0.1,14,info.fAA,kin_visco,WW, Nfft, overlap, info.noise_corr,'sh1',make_plot_spectra,Pf(jp),T1f(jp),folder_out,filename,profID);
                        [BIN.eps_S1(i), BIN.MAD_S1(i), BIN.MADc(i),flag_good_S1, BIN.kL_S1(i),BIN.kU_S1(i),BIN.SPECTRUMS1(i)] = ...
                            TKE_dis_spec(Pf(jp),[sh1_hp(jp) sh2_hp(jp)],AA(jp,:),info.minKS,info.maxKS,info.fAA,kin_visco,WW, Nfft, overlap, info.noise_corr,'sh1',make_plot_spectra,Pf(jp),T1f(jp));
                        [BIN.eps_S2(i), BIN.MAD_S2(i), ~,flag_good_S2,  BIN.kL_S2(i),BIN.kU_S2(i),BIN.SPECTRUMS2(i)] = ...
                            TKE_dis_spec(Pf(jp),[sh1_hp(jp) sh2_hp(jp)],AA(jp,:),info.minKS,info.maxKS,info.fAA,kin_visco,WW, Nfft, overlap, info.noise_corr,'sh2',make_plot_spectra,Pf(jp),T1f(jp));
                        BIN.kB_S1(i)=1/(2*pi())*(BIN.eps_S1(i)/(kin_visco*D^2))^(1/4);
                        BIN.kB_S2(i)=1/(2*pi())*(BIN.eps_S2(i)/(kin_visco*D^2))^(1/4);
                        
                        % Flag indices defined as =1 for bad data (opposite of TKE_dis_spec):
                        if ~isnan(flag_good_S1)
                            BIN.flag_S1(i)=~flag_good_S1;
                        end
                        if ~isnan(flag_good_S2)
                            BIN.flag_S2(i)=~flag_good_S2;
                        end
                        
                        % Batchelor wavenumber determined from shear probe used
                        % to calculate Xi_ST:
                        if (BIN.flag_S1(i)==0 && BIN.flag_S2(i)==1)
                            meanKBSH=BIN.kB_S1(i);
                        elseif (BIN.flag_S1(i)==1 && BIN.flag_S2(i)==0)
                            meanKBSH=BIN.kB_S2(i);
                        else % if both accepted or both rejected (anyway I want to have a value)
                            meanKBSH=mean([BIN.kB_S1(i), BIN.kB_S2(i)]);
                        end
                    catch
                        warning('Shear (S1,S2) spectral calculations did not work in this bin')
                    end
                else
                    try
                        % [BIN.eps_S1(i), BIN.MAD_S1(i), BIN.MADc(i),BIN.flag_S1(i),  BIN.kL_S1(i),BIN.kU_S1(i)] = ...
                        %     TKE_dis_spec(Pf(jp),sh1_hp(jp),AA(jp),0.1,14,info.fAA,kin_visco,WW, Nfft, overlap, info.noise_corr,'sh_1',make_plot_spectra,Pf(jp),T1f(jp),folder_out,filename,profID);
                        [BIN.eps_S1(i), BIN.MAD_S1(i), BIN.MADc(i),flag_good_S1,  BIN.kL_S1(i),BIN.kU_S1(i),BIN.SPECTRUMS1(i)] = ...
                            TKE_dis_spec(Pf(jp),sh1_hp(jp),AA(jp,:),info.minKS,info.maxKS,info.fAA,kin_visco,WW, Nfft, overlap, info.noise_corr,'sh_1',make_plot_spectra,Pf(jp),T1f(jp));
                        
                        % Flag indices defined as =1 for bad data (opposite of TKE_dis_spec):
                        if ~isnan(flag_good_S1)
                            BIN.flag_S1(i)=~flag_good_S1;
                        end
                        
                        
                        BIN.kB_S1(i)=1/(2*pi())*(BIN.eps_S1(i)/(kin_visco*D^2))^(1/4);
                        meanKBSH=BIN.kB_S1(i); % Batchelor wavenumber determined from shear probe used to calculate Xi_ST
                    catch
                        warning('Shear (S1) spectral calculations did not work in this bin')
                    end
                end
            end
        else
            error('Nasmyth_spec is not EPFL')
        end
    end
    %% FP07 spectral calculations
    if ~exist('meanKBSH','var') % Batchelor wavenumber has not been calculated from shear probe --> set it to zero
        meanKBSH=0;
    end

    if  run_dissip
        if param.config.T1
            % [BIN.Xiv1(i),BIN.Xi_ST1(i),BIN.Xi_T1(i),BIN.kB_T1(i),BIN.eps_T1(i),BIN.MAD_ST1(i),BIN.MAD_T1(i),~,BIN.LR_T1(i),BIN.kL_T1(i),BIN.kU_T1(i),BIN.krange_T1(i), BIN.kpeak_T1(i),BIN.flag_T1(i)] =...
            %     gradT_dis_spec(Pf(jp),gradT1f(jp),info.minKT,info.fAA,meanKBSH,WW, Nfft, overlap,info.Tspec,info.q,info.time_res,info.time_corr,info.npoles,info.int_range,D,kin_visco,T1_dT1,'T1_dT1',DATA.setupfilestr,make_plot_spectra,Pf(jp),T1f(jp),folder_out,filename,profID);
            try
                [BIN.Xiv1(i),BIN.Xi_ST1(i),BIN.Xi_T1(i),BIN.kB_T1(i),BIN.eps_T1(i),BIN.MAD_ST1(i),BIN.MAD_T1(i),BIN.MADc(i),BIN.LR_T1(i),BIN.kL_T1(i),BIN.kU_T1(i),BIN.krange_T1(i), BIN.kpeak_T1(i),flag_good_T1,BIN.SPECTRUMT1(i)] =...
                    gradT_dis_spec(Pf(jp),gradT1f(jp),info.minKT,info.fAA,meanKBSH,WW, Nfft, overlap,info.Tspec,info.q,info.time_res,info.time_corr,info.npoles,info.int_range,D,kin_visco,T1_dT1,'T1_dT1',DATA.setupfilestr,make_plot_spectra,Pf(jp),T1f(jp),info.ksfact,info.Snfact);
                
                % Flag indices defined as =1 for bad data (opposite of gradT_dis_spec):
                if ~isnan(flag_good_T1)
                    BIN.flag_T1(i)=~flag_good_T1;
                end
                
                BIN.eps_T1(i) = kin_visco*D^2*(2*pi()*BIN.kB_T1(i))^4;
                BIN.epsT1max(i) = kin_visco*D^2*(2*pi()*info.fAA/WW*info.kmax_factor)^4;
            catch
                if show_progress
                    warning('FP07-T1 spectral calculations did not work in this bin')
                end
            end
        end

        if param.config.T2
            try
                % [BIN.Xiv2(i),BIN.Xi_ST2(i),BIN.Xi_T2(i),BIN.kB_T2(i),BIN.eps_T2(i),BIN.MAD_ST2(i),BIN.MAD_T2(i),~,BIN.LR_T2(i),BIN.kL_T2(i),BIN.kU_T2(i),BIN.krange_T2(i), BIN.kpeak_T2(i),BIN.flag_T2(i)] =...
                %     gradT_dis_spec(Pf(jp),gradT2f(jp),info.minKT,info.fAA,meanKBSH,WW, Nfft, overlap,info.Tspec,info.q,info.time_res,info.time_corr,info.npoles,info.int_range,D,kin_visco,T2_dT2,'T2_dT2',DATA.setupfilestr,make_plot_spectra,Pf(jp),T2f(jp),folder_out,filename,profID);
                %
                [BIN.Xiv2(i),BIN.Xi_ST2(i),BIN.Xi_T2(i),BIN.kB_T2(i),BIN.eps_T2(i),BIN.MAD_ST2(i),BIN.MAD_T2(i),MADc_val,BIN.LR_T2(i),BIN.kL_T2(i),BIN.kU_T2(i),BIN.krange_T2(i), BIN.kpeak_T2(i),flag_good_T2,BIN.SPECTRUMT2(i)] =...
                    gradT_dis_spec(Pf(jp),gradT2f(jp),info.minKT,info.fAA,meanKBSH,WW, Nfft, overlap,info.Tspec,info.q,info.time_res,info.time_corr,info.npoles,info.int_range,D,kin_visco,T2_dT2,'T2_dT2',DATA.setupfilestr,make_plot_spectra,Pf(jp),T2f(jp),info.ksfact,info.Snfact);
                % Flag indices defined as =1 for bad data (opposite of gradT_dis_spec):
                if ~isnan(flag_good_T2)
                    BIN.flag_T2(i)=~flag_good_T2;
                end
                
                if isnan(BIN.MADc(i))
                    BIN.MADc(i)=MADc_val; % Use the threshold value from T2
                end
                BIN.eps_T2(i) = kin_visco*D^2*(2*pi()*BIN.kB_T2(i))^4;
                BIN.epsT2max(i) = kin_visco*D^2*(2*pi()*info.fAA/WW*info.kmax_factor)^4;
            catch
                if show_progress
                    warning('FP07-T2 spectral calculations did not work in this bin')
                end
            end
        end

        

        % Save all open figures
        hfig = get(0, 'Children');
        for kf=1:length(hfig)
            path_spectra=[folder_out,'Figures/spectra_profile_',num2str(kprof,'%02d')];
            if ~exist(path_spectra,"dir")
                mkdir(path_spectra) % Create folder
            end
            % saveas(gcf,[path_spectra,'/bin',num2str(i),'_fig',num2str(kf),'.fig']);
            exportgraphics(gcf,[path_spectra,'/bin',num2str(i),'_press',num2str(BIN.pressure(i)),'_fig',num2str(kf),'.png'],'Resolution',400);


        end


    end
end



%% Flagging of shear and FP07 data (spectral analysis)
if run_dissip
    BIN.flag_ST1=0*BIN.flag_T1;
    if param.config.S2
        idx=find( BIN.MAD_ST1>=2*BIN.MADc | (BIN.flag_S1==1 &  BIN.flag_S2==1) );
    else
        idx=find( BIN.MAD_ST1>=2*BIN.MADc | BIN.flag_S1==1 );
    end
    BIN.flag_ST1(idx)=1;
    idx=find(isnan(BIN.Xi_ST1) & BIN.flag_ST1>=0);
    BIN.flag_ST1(idx)=NaN;

    BIN.flag_ST2=0*BIN.flag_T2;
    if param.config.S2
        idx=find( BIN.MAD_ST2>=2*BIN.MADc | (BIN.flag_S1==1 &  BIN.flag_S2==1) );
    else
        idx=find( BIN.MAD_ST2>=2*BIN.MADc | BIN.flag_S1==1 );
    end
    BIN.flag_ST2(idx)=1;
    idx=find(isnan(BIN.Xi_ST2) & BIN.flag_ST2>=0);
    BIN.flag_ST2(idx)=NaN;

    % Create txt file if flag(s) exist(s)
    % Number of bins with bad data:
    flagvib = length(BIN.flag_vibration(BIN.flag_vibration==1));
    flaginc = length(BIN.flag_inclination(BIN.flag_inclination==1));
    flagacc = length(BIN.flag_acceleration(BIN.flag_acceleration==1));
    flagsh1 = length(BIN.flag_S1(BIN.flag_S1==1));
    flagT1 = length(BIN.flag_T1(BIN.flag_T1==1));
    flagT2 = length(BIN.flag_T2(BIN.flag_T2==1));
    flagST1 = length(BIN.flag_ST1(BIN.flag_T1==1));
    flagST2 = length(BIN.flag_ST2(BIN.flag_T2==1));
    flagsh2 = 0; %% Define it to be 0 if only sh1 is present, so no conflict below
    if param.config.S2
        flagsh2 = length(BIN.flag_S2(BIN.flag_S2==1));
    end
    if flagvib == 0 & flaginc == 0 & flagacc == 0 & flagT1 == 0 & flagT2 == 0 & flagsh1 == 0 & flagsh2 == 0 & flagST1 == 0 & flagST2 == 0
        %no flags, do nothing
    else
        flagname = sprintf('flag_profile%.0f.txt', kprof);
        fido = fopen(strcat(folder_out,flagname),'wt');
        if flagvib > 0
            fprintf(fido,'flagged due to excessive vibrations\n');
        end
        if flaginc > 0
            fprintf(fido,'flagged due to excessive inclination\n');
        end
        if flagacc > 0
            fprintf(fido,'flagged due to excessive acceleration of instrumet\n');
        end
        if flagsh1 > 0
            fprintf(fido,'flagged for S1 fit\n');
        end
        if flagsh2 > 0
            fprintf(fido,'flagged for S2 fit\n');
        end
        if flagT1 > 0
            fprintf(fido,'flagged for T1 fit\n');
        end
        if flagT2 > 0
            fprintf(fido,'flagged for T2 fit\n');
        end
        if flagST1 > 0
            fprintf(fido,'flagged for ST1\n');
        end
        if flagST2 > 0
            fprintf(fido,'flagged for ST2');
        end
        fclose(fido);
    end

    % Display flagging % data removed:
    disp('************************************')
    fprintf('Percentage of data to keep: T1 = %0.2f%%, T2 = %0.2f%%, S1 = %0.2f%%, S2 = %0.2f%%\n', ...
        sum(BIN.flag_T1==0)/sum(~isnan(BIN.depth))*100,...
        sum(BIN.flag_T2==0)/sum(~isnan(BIN.depth))*100,...
        sum(BIN.flag_S1==0)/sum(~isnan(BIN.depth))*100,...
        sum(BIN.flag_S2==0)/sum(~isnan(BIN.depth))*100);
     disp('************************************')
    %% Additional variables: K-OsbornCox and Ozmidov length scale
    BIN.KOsborn_S1 = 0.2*BIN.eps_S1.*(BIN.N2).^-1; % Gamma = 0.2
    BIN.KOsborn_T1=0.2*BIN.eps_T1.*(BIN.N2).^-1;
    BIN.KOsborn_T2=0.2*BIN.eps_T2.*(BIN.N2).^-1;
    BIN.KOsbornCox_T1 = 0.5*BIN.Xi_T1.*(BIN.avggradT1).^-2;
    BIN.KOsbornCox_T2 = 0.5*BIN.Xi_T2.*(BIN.avggradT2).^-2;
    BIN.KOsbornCox_ST1 = 0.5*BIN.Xi_ST1.*(BIN.avggradT1).^-2;
    BIN.KOsbornCox_ST2 = 0.5*BIN.Xi_ST2.*(BIN.avggradT2).^-2;
    BIN.LO_S1 = (BIN.eps_S1./BIN.N2.^(3/2)).^(0.5);
    if param.config.S2
        BIN.LO_S2 = (BIN.eps_S2./BIN.N2.^(3/2)).^(0.5);
        BIN.KOsborn_S2 = 0.2*BIN.eps_S2.*(BIN.N2).^-1; % Gamma = 0.2
    end


    %% Plots bin profile
    if make_plot_prof
        plot_bin_profile(SLOW,BIN,param)

        if ~exist([folder_out,'Figures'],'dir')
            mkdir([folder_out,'Figures'])
        end
        saveas(gcf,[folder_out,'Figures/results_',num2str(kprof,'%02d'),info.prof_dir,'.fig'])
        exportgraphics(gcf,[folder_out,'Figures/results_',num2str(kprof,'%02d'),info.prof_dir,'.png'],'Resolution',400)
    end
end

%% Add metadata
% uncomment longitude/latitude entries once they are programmed
%        BIN.latitude = lat;
%        BIN.longitude = lon;
%        FAST.latitude = lat;
%        FAST.longitude = lon;
%        SLOW.latitude = lat;
%        SLOW.longitude = lon;
BIN.bin_size = info.dpD;
BIN.bin_overlap = info.dp;
BIN.profID=profID;
SLOW.profID=profID;
FAST.profID=profID;
FAST.date = BIN.date;
SLOW.date = BIN.date;
FAST.time = BIN.time;
SLOW.time = BIN.time;
FAST.direction = BIN.direction;
SLOW.direction = BIN.direction;


end
