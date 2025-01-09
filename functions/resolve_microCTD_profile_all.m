function [BIN,SLOW,FAST] = resolve_microCTD_profile(DATA, inp, info,dataf_name,PLOT,folder_out,profID)
    despike_sh  = [ 8  0.5 0.04];
    despike_A = [8 0.5000 0.0400];
    if ~ isfield(info,'system')
        info.system = 'Oce';
    end
    if ~ isfield(info,'minvel_detect')
        info.minvel_detect = 0.1;
    end
    if ~ isfield(info,'mindur_detect')
        info.mindur_detect = 30;
    end
    if ~ isfield(info,'pmin')
        info.pmin = 1;
    end
    if ~ isfield(info,'dp')
        info.pint = 1;
    end
    if ~ isfield(info,'dpD')
        info.pintD = 2;
    end
    if ~ isfield(info,'prof_dir')
        info.prof_dir = 'down';
    end
    if ~ isfield(info,'k_HP_cut')
        % SEB: low-frequency motions of the free-falling profiler (fluctuations with the scales corresponding 
        % to half of the instrument length) can be a source of low-frequency noise in the microstructure 
        % shear signa DOI: 10.1016/j.pocean.2006.07.003. MicroCTD length = 1 m --> k_HP_cut = 0.5 cpm
        info.k_HP_cut = 0.5;
    end
    if ~ isfield(info,'minKT')
        info.minKT = 1;
    end
    if ~ isfield(info,'fAA')
        info.fAA = 110;
    end
    if ~ isfield(info,'Tmethod')
        info.Tmethod = 'B';
    end
    if ~ isfield(info,'Tspec')
        info.Tspec = 'B';
    end
    if ~ isfield(info,'noisep_T1')
        info.noisep_T1 = [-10.24,-0.89,info.fAA];
    else
        info.noisep_T1 = [info.noisep_T1, info.fAA];
    end
    if ~ isfield(info,'noisep_T2')
        info.noisep_T2 = [-10.08,-0.97,info.fAA];
    else
        info.noisep_T2 = [info.noisep_T2, info.fAA];
    end
    if ~isfield(info,'time_res')
        info.time_res = nan;
    end
    if ~ isfield(info,'peak_rem_T1')
        info.peak_rem_T1 = [0,0];
    end
    if ~ isfield(info,'peak_rem_T2')
        info.peak_rem_T2 = [0,0];
    end
    if ~ isfield(info,'kmax_fac')
        info.kmax_factor = 1/1.66;
    end

    %adds some more noise
    info.noisep_T1 = info.noisep_T1 + [0.25,0,0];
    info.noisep_T2 = info.noisep_T2 + [0.25,0,0];
        
    %constants
    %visco = 1e-6;
    %SEB: changed D (see below)
    %D = 1.44e-7; 

    
    %defines times
    time_fast0 = [0:1:length(DATA.P_fast)-1]/DATA.fs_fast;
    time_slow0= [0:1:length(DATA.P_slow)-1]/DATA.fs_slow;

    %gets profiles
    iPs0 = get_profile(DATA.P_slow,DATA.W_slow,0,info.minvel_detect,info.prof_dir,info.mindur_detect,DATA.fs_slow);
    iPf0 = get_profile(DATA.P_fast,DATA.W_fast,0,info.minvel_detect,info.prof_dir,info.mindur_detect,DATA.fs_fast);
    NP = size(iPf0);
    NP = NP(2);
    if inp<0 || inp>NP
        fprintf('Wrong number of profiles')
        return
    end
    
    %gets index for the desired profile
    if NP>0
        iipf = iPf0(1,inp):iPf0(2,inp);
        iips = iPs0(1,inp):iPs0(2,inp);
    else
        iipf = [1:1:length(DATA.P_fast)];
        iips = [1:1:length(DATA.P_slow)-1];
    end
    
    %date
    date = datenum(DATA.Year, DATA.Month, DATA.Day, DATA.Hour, DATA.Minute, DATA.Second );
    date = date + iPs0(1,inp)/DATA.fs_slow/60/60/24;
    datestr(date)
    
    %sampling frequencies
    fss = DATA.fs_slow;
    fsf = DATA.fs_fast;  

    %gets fast response sensors
    timef = time_fast0(iipf);
    Pf = DATA.P_fast(iipf);
    T1f = DATA.T1_fast(iipf);
    T2f = DATA.T2_fast(iipf);
    T1_dT1 = DATA.T1_dT1(iipf);
    T2_dT2 = DATA.T2_dT2(iipf);    
    gradT1f = DATA.gradT1(iipf);
    gradT2f = DATA.gradT2(iipf);
    sh1 = DATA.sh1(iipf);
    if isfield(DATA,'sh2')
        sh2 = DATA.sh2(iipf);
    end
    Ax = DATA.Ax(iipf);
    Ay = DATA.Ay(iipf);
    Wf = DATA.W_fast(iipf);
    AA = [Ax,Ay];
    if any(strcmpi(fieldnames(DATA),'Chlorophyll')) % MicroCTD
        Chl = DATA.Chlorophyll(iipf);
        Turb = DATA.Turbidity(iipf);
    elseif any(strcmpi(fieldnames(DATA),'NTU')) % VMP
        Chl = zeros(size(Ax));
        Turb = zeros(size(Ax));
    end
    
    %gets slow response sensors
    times = time_slow0(iips);
    Ws = DATA.W_slow(iips);
    Ps = DATA.P_slow(iips);
    if any(strcmpi(fieldnames(DATA),'JAC_T')) % MicroCTD
        T_JAC = DATA.JAC_T(iips);
        C_JAC = DATA.JAC_C(iips);
        Incl_X = DATA.Incl_X(iips);
        Incl_Y = DATA.Incl_Y(iips);
    elseif any(strcmpi(fieldnames(DATA),'SBT')) % VMP
        T_JAC = DATA.SBT(iips);
        C_JAC = DATA.SBC(iips);
        Incl_X = zeros(size(T_JAC));
        Incl_Y = zeros(size(T_JAC));
    end

    %accurate T in fast sensors grid
    T_JAC_fast = interp1(Ps,T_JAC,Pf);
    i1 = find(isfinite(T_JAC_fast),1,'first');
    if i1>1
        T_JAC_fast(1:i1-1) = T_JAC_fast(i1);
    end
    i2 = find(isfinite(T_JAC_fast),1,'last');
    if i2<length(T_JAC_fast)
        T_JAC_fast(i2+1:end) = T_JAC_fast(i2);
    end

    %sets maximum depth
    if ~ isfield(info,'pmax')
        info.pmax = round(max(Ps)-info.dpD/2)+1;
    end

    %calculates salinity and density
    %[rhos,Ss,depths]=rho_salinity_Geneva(T_JAC,C_JAC,Ps);
    [Ss,Cmatch] = salinity_JAC(Ps,T_JAC,C_JAC);

    sgt = sigmat(T_JAC, Ss)+1000;
    mrho = cumsum(sgt)./[1:length(sgt)]';
    depth=10000*Ps./(mrho*9.81);
    
    if info.system == 'Lem' | info.system == 'Rot'
        if info.system == 'Lem'
            %density for Leman
            [sgt,Ss,depth]=rho_salinity_Geneva(T_JAC,Cmatch*1000,Ps);
        elseif info.system == 'Rot'
            %density for Rotsee
            [sgt,Ss,depth]=rho_salinity_Rotsee(T_JAC,Cmatch*1000,Ps);
        end
        [sgtT1,SsT1,depthT1]=rho_salinity_Geneva(T1f,interp1(Ps,Cmatch*1000,Pf,'linear','extrap'),Pf);
        [sgtT2,SsT2,depthT2]=rho_salinity_Geneva(T2f,interp1(Ps,Cmatch*1000,Pf,'linear','extrap'),Pf);
    elseif info.system == 'Zue'
        [sgt,Ss,depth]=rho_salinity_Zurich(T_JAC,Cmatch*1000,Ps);
    elseif info.system == 'Gar'
        [sgt,Ss,depth]=rho_salinity_Garda(T_JAC,Cmatch*1000,Ps);
    end
    
    %Fast for output
    FAST.datenum = date;
    FAST.filename = dataf_name;
    if ismember(info.prof_dir,'up')
        FAST.upward = 1;
    else
        FAST.upward = 0;
    end
    FAST.pressure = Pf;
    FAST.fast_T1 = T1f;
    FAST.fast_T2 = T2f;
    FAST.grad_T1 = gradT1f;
    FAST.grad_T2 = gradT2f;    
    FAST.chlorophyll = Chl;
    FAST.turbidity = Turb;
    FAST.fast_S1 = sh1;
    FAST.A_x = Ax;
    FAST.A_y = Ay;
    if isfield(DATA,'sh2')
        FAST.fast_S2 = sh2;
    else
        FAST.depth = interp1(Ps,depth,Pf,'linear','extrap'); % Good compromise!
        FAST.sigmatT1 = sgtT1;
        FAST.sigmatT2 = sgtT2;
    end
        
    %slow data for output
    SLOW.datenum = date;
    SLOW.filename = dataf_name;
    if ismember(info.prof_dir,'up')
        SLOW.upward = 1;
    else
        SLOW.upward = 0;
    end
    SLOW.pressure = Ps;
    SLOW.depth = depth;
    SLOW.temperature = T_JAC;
    SLOW.conductivity = Cmatch*1000; %uS/cm
    SLOW.salinity = 1000*Ss; %mg/l
    SLOW.density = sgt;
    SLOW.Incl_x = Incl_X;
    SLOW.Incl_y = Incl_Y;
    
    %calculates displacements for thorpe length
    if ismember(info.prof_dir,'down')
        [sort_rho, isd] = sort(sgt);
        displ = Ps - Ps(isd);
        
        [sort_uT1, ist1] = sort(T1f,'descend');
        displuT1 = Pf - Pf(ist1);
        
        [sort_uT2, ist2] = sort(T2f,'descend');
        displuT2 = Pf - Pf(ist2);
        
        [sort_T_JAC, ~] = sort(T_JAC,'descend');        
    else
        [sort_rho, isd] = sort(sgt,'descend');
        displ = Ps - Ps(isd);
        
        [sort_uT1, ist1] = sort(T1f,'ascend');
        displuT1 = Pf - Pf(ist1);
        
        [sort_uT2, ist2] = sort(T2f,'ascend');
        displuT2 = Pf - Pf(ist2);

        [sort_T_JAC, ~] = sort(T_JAC,'ascend');                
    end
    % centered Thorpe scale
    Ps_centered=Ps-displ*0.5;
    Lc=zeros(1,length(Ps)); count_c=zeros(1,length(Ps));
    for i_Ps=1:length(Ps)
        [~,i_c] = min(abs(Ps-Ps_centered(i_Ps)));
        Lc(i_c) = Lc(i_c) + abs(displ(i_Ps));
        count_c(i_c) = count_c(i_c) + 1;
    end
    Lc = Lc./count_c; Lc(count_c==0)=0;
    SLOW.LTc = Lc'; 
    

    %binned temperature, salinity and density
    %defines the pressure vector where to calculate
    pplot=[info.pmin info.pmax];
    pres = [info.pmin:info.dp:info.pmax];
    if ismember(info.prof_dir,'up') % For upward profiles, start the profile after the acceleration
        idx75=find(Ps<0.75*info.pmax,1,'first'); % Note: Ps is decreasing (upward profile). 
        [~,idx]=min(Ws(1:idx75)); % Focus just on the first 25% (due to possible wave-induced oscillations)
        n_pres = find(pres<=Ps(idx),1,'last');
    else
        n_pres=length(pres);
    end
    pmaxplot=pres(n_pres);
 
    BIN.datenum = date;
    BIN.date=cellstr(datestr(BIN.datenum,'yyyy-mm-dd'));
    BIN.time=cellstr(datestr(BIN.datenum,'HH:MM:SS'));
    BIN.filename = dataf_name;
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
    

     FAST.filename = BIN.filename;
     SLOW.filename = BIN.filename;

    
    
    

    BIN.pressure = pres;
    BIN.depth = pres_av(Ps,depth,pres,info.dp,2.7);
    BIN.temperature = pres_av(Ps,T_JAC,pres,info.dp,2.7);
    BIN.conductivity = 1000*pres_av(Ps,C_JAC,pres,info.dp,2.7); % uS/cm
    BIN.salinity = 1000*pres_av(Ps,Ss,pres,info.dp,2.7); %mg/l
    BIN.density = pres_av(Ps,sgt,pres,info.dp,2.7);

    Thorpe=sqrt(pres_av(Ps,displ.^2,pres,info.dpD,0));
    for i=1:length(pres)
        % It is better to use FP07 because when the lake is homogeneous,
        % the resolution of JAC-T is too low, resulting in grT=1E-16
        BIN.avggradT1(i) = mean_grad(Pf,sort_uT1,pres(i),max(Thorpe(i),info.dpD));   % SEB: evaluated on the resorted T profile
        BIN.avggradT2(i) = mean_grad(Pf,sort_uT2,pres(i),max(Thorpe(i),info.dpD));   % SEB: evaluated on the resorted T profile
        BIN.grT12(i)=nanmean([BIN.avggradT1(i),  BIN.avggradT2(i)]);
        BIN.avggradT(i) = mean_grad(Ps,sort_T_JAC,pres(i),max(Thorpe(i),info.dpD));   % SEB: evaluated on the resorted T profile
        BIN.N2(i) = -9.81*mean_grad(Ps,sort_rho,pres(i),max(Thorpe(i),info.dpD))/1000;
    end
    BIN.LT = Thorpe; 
    BIN.chlorophyll = pres_av(Pf,Chl,pres,info.dp,2.7);
    BIN.turbidity = pres_av(Pf,Turb,pres,info.dp,2.7);
    BIN.LTuT1 = sqrt(pres_av(Pf,displuT1.^2,pres,info.dpD,0));
    BIN.LTuT2 = sqrt(pres_av(Pf,displuT2.^2,pres,info.dpD,0));
  
    BIN.acc = -mean_grad_time(Pf,Wf,pres,info.dpD,timef);  % instrument acceleration
    

    %plots raw 
    figure(1)
    clf
    subplot(4,4,[1 4])
    plot(DATA.t_slow/60,DATA.P_slow)
    ylabel('p (db)')
    hold on
    plot(DATA.t_slow(iips)/60,Ps); ylim(pplot)
    xlabel('time [min]')
    subplot(4,4,[5,9,13])
    plot(sh1,Pf); hold on; 
    if isfield(DATA,'sh2')
        plot(sh2+10,Pf); hold on; 
    end
    plot(get(gca,'xlim'),[pmaxplot pmaxplot],'--k')
    set(gca,'ydir','reverse');
    ylabel('p(db)'); xlabel('sh (1/s)'); ylim(pplot);
%     subplot(4,4,[6,10,14])
%     plot(AA(:,1),Pf); hold on
%     plot(AA(:,2),Pf); set(gca,'ydir','reverse');
%     xlabel('Acc (1/s)')
    subplot(4,4,[7,11,15]);
    if isfield(DATA,'sh2')
        plot(gradT1f,Pf)
        hold on
        plot(gradT1f+10,Pf);
    else
        plot(T1f,Pf)
        hold on
        plot(T2f,Pf); 
    end
    plot(get(gca,'xlim'),[pmaxplot pmaxplot],'--k')
    set(gca,'ydir','reverse');
    xlabel('T fast (°C)'); ylim(pplot)
    set(gca,'yticklabel',[])
    ax1=subplot(4,4,[8,12,16]);     
    plot(Wf,Pf,'-k'); hold on;    
    plot(get(gca,'xlim'),[pmaxplot pmaxplot],'--k')
    set(gca,'ydir','reverse');
    xlabel('W (m/s)');yticklabels([]);ylim(pplot);
    ax2=axes('Position',ax1.Position,'XAxisLocation','top',...
        'YAxisLocation','right','color','none',...
        'xColor','r','yColor','k');   
    set(ax1,'box','off')    
    line(Incl_X,DATA.P_slow(iips),'color','r'); set(gca,'ydir','reverse');
    xlabel('Inclination [°]')
    yticklabels([]); xlim([-5 5]);

    %filters shear and highpasses microstructure
    [sh1, ~, ~, ~ ] =  despike(sh1, despike_sh(1), despike_sh(2), fsf, round(despike_sh(3)*fsf));
    if isfield(DATA,'sh2')
        [sh2, ~, ~, ~ ] =  despike(sh2, despike_sh(1), despike_sh(2), fsf, round(despike_sh(3)*fsf));
    end
    
    mW = (max(Pf)-min(Pf))/(max(timef)-min(timef));
    if info.k_HP_cut>0
        f_HP_cut = info.k_HP_cut*mW;
        [bh,ah] = butter(1, f_HP_cut/(fsf/2), 'high');
        sh1_hp = filter(bh, ah, sh1);
        sh1_hp = flipud(sh1_hp);
        sh1_hp = filter(bh, ah, sh1_hp);
        sh1_hp = flipud(sh1_hp);
        if isfield(DATA,'sh2')
            sh2_hp = filter(bh, ah, sh2);
            sh2_hp = flipud(sh2_hp);
            sh2_hp = filter(bh, ah, sh2_hp);
            sh2_hp = flipud(sh2_hp);
        end
        
        T2f_hp=T2f;
        T1f_hp=T1f;
        
    else
        sh1_hp = sh1;
        if isfield(DATA,'sh2')
            sh2_hp = sh2;
        end
        T1f_hp = T1f;
        T2f_hp = T2f;
    end
    
    
    %-- despike the piezo-accelerometer signals
    piezo_accel_num = size(AA,2);
    if  ~isempty(AA) && despike_A(1) ~= inf
        for probe = 1:piezo_accel_num
            [AA(:,probe), ~, ~, ~]  = ...
                despike(AA(:,probe),  despike_A(1), ...
                despike_A(2), fsf, round(despike_A(3)*fsf));
        end
    end
    
    if ~isfield(DATA,'sh2')
            %-- identify the piezo-accelerometer signal to be used for noise filtering
        PSDsh = csd_odas(detrend(sh1_hp),detrend(sh1_hp),1024,fsf,[],512,'linear');
        PSDA1 = csd_odas(detrend(AA(:,1)),detrend(AA(:,1)),1024,fsf,[],512,'linear');
        PSDA2 = csd_odas(detrend(AA(:,2)),detrend(AA(:,2)),1024,fsf,[],512,'linear');
        CSDshA1 = csd_odas(detrend(sh1_hp),detrend(AA(:,1)),1024,fsf,[],512,'linear');
        [CSDshA2,fA] = csd_odas(detrend(sh1_hp),detrend(AA(:,2)),1024,fsf,[],512,'linear');
        COHshA1 = abs(CSDshA1.^2./(PSDA1.*PSDsh));
        COHshA2 = abs(CSDshA2.^2./(PSDA2.*PSDsh));
        AAxy=AA; % store both signals
        if nanmean(COHshA1)>nanmean(COHshA2)
            AA = AA(:,1);
        else
            AA = AA(:,2);
        end
    end
    % add despiked and HP signals to the plot
    subplot(4,4,[5,9,13])
    plot(sh1,Pf)
        
    if isfield(DATA,'sh2')
        plot(sh2+10,Pf)
    end

    subplot(4,4,[6,10,14])
    plot(AA,Pf); hold on
    set(gca,'ydir','reverse');
    plot([-200 -200],[min(Pf) max(Pf)],'--k'); plot([200 200],[min(Pf) max(Pf)],'--k');    
    xlim([-500 500]);  xlabel('Acc (1/s)'); ylim(pplot);
    plot(get(gca,'xlim'),[pmaxplot pmaxplot],'--k')
    saveas(gcf,[folder_out,'/profile',num2str(inp,'%02d'),'.png'])
    
    
    %defines output variables
    BIN.speed = nan(1,length(pres));
    BIN.maxgradT = nan(1,length(pres));
    BIN.eps_S1 = nan(1,length(pres));
    BIN.kB_S1 = nan(1,length(pres));
    BIN.MAD_S1 = nan(1,length(pres));
    BIN.MADc = nan(1,length(pres));
    BIN.flag_S1 = nan(1,length(pres));
	BIN.kL_S1 = nan(1,length(pres));
    BIN.kU_S1 = nan(1,length(pres));
    
    if isfield(DATA,'sh2')
        BIN.eps_S2 = nan(1,length(pres));
        BIN.kB_S2 = nan(1,length(pres));
        BIN.MAD_S2 = nan(1,length(pres));
        BIN.flag_S2 = nan(1,length(pres));
        BIN.kL_S2 = nan(1,length(pres));
        BIN.kU_S2 = nan(1,length(pres));
    end    
    
    BIN.Xi_ST1 = nan(1,length(pres));
    BIN.Xi_T1 = nan(1,length(pres));
    BIN.kB_T1 = nan(1,length(pres));
    BIN.sXif1 = nan(1,length(pres));
    BIN.sKB1 = nan(1,length(pres));
    BIN.Xiv1 = nan(1,length(pres));
    BIN.kU_T1 = nan(1,length(pres));
    BIN.eps_T1 = nan(1,length(pres));
    BIN.epsT1max = nan(1,length(pres));
    BIN.MAD_ST1 = nan(1,length(pres));
    BIN.MAD_T1 = nan(1,length(pres));
    BIN.LKH1 = nan(1,length(pres));
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
      
    for i = 1:n_pres; %length(pres)
        jp = find(Pf>=pres(i)-info.dpD/2 & Pf<=pres(i)+info.dpD/2);
        jps = find(Ps>=pres(i)-info.dpD/2 & Ps<=pres(i)+info.dpD/2);   % slow channel
        
        % SEB: added evaluation of "maximum" gradT for pyroelectric analysis
        if ~isempty(jps) & length(jps)>=3
            tmp = movavg(T_JAC(jps),'linear',3);
            BIN.maxgradT(i) = max(-gradient(tmp)./gradient(Ps(jps)));
        end
        
        % SEB: added evaluation of D and viscosity (based on http://web.mit.edu/seawater/)
        avT = nanmean(T_JAC(jps));
        avS = nanmean(Ss(jps));
        avsigmat = nanmean(sgt(jps));
        cond = SW_Conductivity(avT,'C',avS,'ppt');
        cp = SW_SpcHeat(avT,'C',avS,'ppt',0.101325,'MPa');   % The last argument is sat vapor pressure. It can kept constant as suggested in SW_Diffusivity.m
        D = cond/(avsigmat*cp);
    
        visco = SW_Viscosity(avT,'C',avS,'ppt'); % dynamic viscosity mu
        visco = visco/avsigmat; % kinematic viscosity nu
        
        if isfield(info,'num_fft')
           Nfft = floor(length(jp)/((info.num_fft+1)/2));
        else
           Nfft = floor(length(jp)/2);
        end
        
        if isfield(info,'overlap')
           overlap = info.overlap;
        else
           overlap = round(Nfft/2);
        end

        if length(jp)<Nfft | length(jp)<128
            continue
        end
        
        % SEB: calculate the profiling speed to be used below. I prefer this way
        % than using Wf or Ws, since these are smoothed and for short bins
        % may be critical. However, in case just define WW=mean(abs(Wf(jp)))
        % Note that in the main script I removed the smmoothing, however I
        % leave this definition to be safe.
        WW = mean( abs( (Pf(jp(end)) - Pf(jp(1)))/(timef(jp(end)) - timef(jp(1))) ) );
        
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
        len_Axi = nansum(abs(Axi)>200);
        len_Ayi = nansum(abs(Ayi)>200);
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
        
        % flaggin for inclination, same as vibr for incl. > 5°
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
        
        %% ODAS epsilon calculation
       try

            if strcmp(info.Nasmyth_spec,'EPFL')
                if isfield(DATA,'sh2')
                    [BIN.eps_S1(i), BIN.MAD_S1(i), BIN.MADc(i),BIN.flag_S1(i), BIN.kL_S1(i),BIN.kU_S1(i)] = ...
                        TKE_dis_spec(Pf(jp),[sh1_hp(jp) sh2_hp(jp)],AA(jp,:),0.1,14,info.fAA,visco,WW, Nfft, overlap, info.noise_corr,'sh1',PLOT,Pf(jp),T1f(jp),folder_out,dataf_name,profID);
                    [BIN.eps_S2(i), BIN.MAD_S2(i), ~,BIN.flag_S2(i),  BIN.kL_S2(i),BIN.kU_S2(i)] = ...
                        TKE_dis_spec(Pf(jp),[sh1_hp(jp) sh2_hp(jp)],AA(jp,:),0.1,14,info.fAA,visco,WW, Nfft, overlap, info.noise_corr,'sh2',PLOT,Pf(jp),T1f(jp),folder_out,dataf_name,profID);
                else
                    [BIN.eps_S1(i), BIN.MAD_S1(i), BIN.MADc(i),BIN.flag_S1(i),  BIN.kL_S1(i),BIN.kU_S1(i)] = ...
                        TKE_dis_spec(Pf(jp),sh1_hp(jp),AA(jp),0.1,14,info.fAA,visco,WW, Nfft, overlap, info.noise_corr,'sh_1',PLOT,Pf(jp),T1f(jp),folder_out,dataf_name,profID);
                end
            else
                error('Nasmyth_spec is not EPFL')
            end
        end

            BIN.kB_S1(i)=1/(2*pi())*(BIN.eps_S1(i)/(visco*D^2))^(1/4);
            if isfield(DATA,'sh2')
                BIN.kB_S2(i)=1/(2*pi())*(BIN.eps_S2(i)/(visco*D^2))^(1/4);
                if (BIN.flag_S1(i)==0 && BIN.flag_S2(i)==1)
                    meanKBSH=BIN.kB_S1(i);
                elseif (BIN.flag_S1(i)==1 && BIN.flag_S2(i)==0)
                    meanKBSH=BIN.kB_S2(i);
                else % if both accepted or both rejected (anyway I want to have a value)
                    meanKBSH=mean([BIN.kB_S1(i), BIN.kB_S2(i)]);
                end
            else
                meanKBSH=BIN.kB_S1(i);
            end
            

                    %% FP07 calculations
                    
                    
            try
                [BIN.Xiv1(i),BIN.Xi_ST1(i),BIN.Xi_T1(i),BIN.kB_T1(i),BIN.eps_T1(i),BIN.MAD_ST1(i),BIN.MAD_T1(i),~,BIN.LR_T1(i),BIN.kL_T1(i),BIN.kU_T1(i),BIN.krange_T1(i), BIN.kpeak_T1(i),BIN.flag_T1(i)] =...
                    gradT_dis_spec(Pf(jp),gradT1f(jp),info.minKT,info.fAA,meanKBSH,WW, Nfft, overlap,info.Tspec,info.q,info.time_res,info.time_corr,info.npoles,info.int_range,D,visco,T1_dT1,'T1_dT1',DATA.setupfilestr,PLOT,Pf(jp),T1f(jp),folder_out,dataf_name,profID);
                
                BIN.eps_T1(i) = visco*D^2*(2*pi()*BIN.kB_T1(i))^4;
                BIN.epsT1max(i) = visco*D^2*(2*pi()*info.fAA/WW*info.kmax_factor)^4;

                [BIN.Xiv2(i),BIN.Xi_ST2(i),BIN.Xi_T2(i),BIN.kB_T2(i),BIN.eps_T2(i),BIN.MAD_ST2(i),BIN.MAD_T2(i),~,BIN.LR_T2(i),BIN.kL_T2(i),BIN.kU_T2(i),BIN.krange_T2(i), BIN.kpeak_T2(i),BIN.flag_T2(i)] =...
                    gradT_dis_spec(Pf(jp),gradT2f(jp),info.minKT,info.fAA,meanKBSH,WW, Nfft, overlap,info.Tspec,info.q,info.time_res,info.time_corr,info.npoles,info.int_range,D,visco,T2_dT2,'T2_dT2',DATA.setupfilestr,PLOT,Pf(jp),T2f(jp),folder_out,dataf_name,profID);

                BIN.eps_T2(i) = visco*D^2*(2*pi()*BIN.kB_T2(i))^4;
                BIN.epsT2max(i) = visco*D^2*(2*pi()*info.fAA/WW*info.kmax_factor)^4;
                
           end
     end
    

    BIN.flag_ST1=0*BIN.flag_T1;
    if isfield(DATA,'sh2')
        idx=find( BIN.MAD_ST1>=2*BIN.MADc | (BIN.flag_S1==1 &  BIN.flag_S2==1) );
    else
        idx=find( BIN.MAD_ST1>=2*BIN.MADc | BIN.flag_S1==1 );
    end
    BIN.flag_ST1(idx)=1;
    idx=find(isnan(BIN.Xi_ST1) & BIN.flag_ST1>=0);
    BIN.flag_ST1(idx)=NaN;
    
    BIN.flag_ST2=0*BIN.flag_T2;
    if isfield(DATA,'sh2')
        idx=find( BIN.MAD_ST2>=2*BIN.MADc | (BIN.flag_S1==1 &  BIN.flag_S2==1) );
    else
        idx=find( BIN.MAD_ST2>=2*BIN.MADc | BIN.flag_S1==1 );
    end
    BIN.flag_ST2(idx)=1;    
    idx=find(isnan(BIN.Xi_ST2) & BIN.flag_ST2>=0);
    BIN.flag_ST2(idx)=NaN;

    % Create txt file if flag(s) exist(s)
    flagvib = length(BIN.flag_vibration(BIN.flag_vibration==1));
    flaginc = length(BIN.flag_inclination(BIN.flag_inclination==1));
    flagacc = length(BIN.flag_acceleration(BIN.flag_acceleration==1));
    flagsh1 = length(BIN.flag_S1(BIN.flag_S1==1));
    flagT1 = length(BIN.flag_T1(BIN.flag_T1==1));
    flagT2 = length(BIN.flag_T2(BIN.flag_T2==1));
    flagST1 = length(BIN.flag_ST1(BIN.flag_T1==1));
    flagST2 = length(BIN.flag_ST2(BIN.flag_T2==1));
    flagsh2 = 0; %% Define it to be 0 if only sh1 is present, so no conflict below
    if isfield(DATA,'sh2')
        flagsh2 = length(BIN.flag_S2(BIN.flag_S2==1));
    end
    if flagvib == 0 & flaginc == 0 & flagacc == 0 & flagT1 == 0 & flagT2 == 0 & flagsh1 == 0 & flagsh2 == 0 & flagST1 == 0 & flagST2 == 0
        %no flags, do nothing
    else
        flagname = sprintf('flag_profile%.0f.txt', inp);
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
            fprintf(fido,'flagged for S1 fit\n')
        end
        if flagsh2 > 0
            fprintf(fido,'flagged for S2 fit\n')
        end
        if flagT1 > 0
            fprintf(fido,'flagged for T1 fit\n')
        end
        if flagT2 > 0
            fprintf(fido,'flagged for T2 fit\n')
        end
        if flagST1 > 0
            fprintf(fido,'flagged for ST1\n')
        end
        if flagST2 > 0
            fprintf(fido,'flagged for ST2')
        end
        fclose(fido);
    end
    
    %% plots
    BIN.KOsborn_S1 = 0.2*BIN.eps_S1.*(BIN.N2).^-1; % Gamma = 0.2
    BIN.KOsborn_T1=0.2*BIN.eps_T1.*(BIN.N2).^-1;
    BIN.KOsborn_T2=0.2*BIN.eps_T2.*(BIN.N2).^-1;
    BIN.KOsbornCox_T1 = 0.5*BIN.Xi_T1.*(BIN.avggradT1).^-2;
    BIN.KOsbornCox_T2 = 0.5*BIN.Xi_T2.*(BIN.avggradT2).^-2;
    BIN.KOsbornCox_ST1 = 0.5*BIN.Xi_ST1.*(BIN.avggradT1).^-2;
    BIN.KOsbornCox_ST2 = 0.5*BIN.Xi_ST2.*(BIN.avggradT2).^-2;
    BIN.LO_S1 = (BIN.eps_S1./BIN.N2.^(3/2)).^(0.5);
    if isfield(DATA,'sh2')
        BIN.LO_S2 = (BIN.eps_S2./BIN.N2.^(3/2)).^(0.5);
        BIN.KOsborn_S2 = 0.2*BIN.eps_S2.*(BIN.N2).^-1; % Gamma = 0.2
    end
    
    pmin = info.pmin - info.dpD/2;
    pmax = pres( find(isfinite(BIN.temperature),1,'last')) + info.dpD/2;

    %% plots profile
    figure(2)
    clf
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperSize', [29 20]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 29 20]);

    ax1=subplot(2,4,1);
    plot(T_JAC, Ps,'-k','linewidth',1)
    xlabel('T (°C)')
    ylabel('p (db)')
    set(gca,'YDir','reverse')
    ylim([pmin,pmax])
    ax2=axes('Position',ax1.Position);
    set(ax1,'box','off')
    plot(Ss, Ps,'r', 'parent' , ax2,'linewidth',1)
    set(ax2,'XAxisLocation','top',...
            'YAxisLocation','right',...
            'Color','none',...
            'XColor','r','YColor','k');
    yticklabels([])
    ylim([pmin,pmax])
    xlabel('Salinity')
    set(gca,'yticklabel',[])
    set(gca,'YDir','reverse')

    ax3=subplot(2,4,2);
    plot(BIN.speed,BIN.pressure,'.-k','linewidth',1,'markersize',4)
    xlabel('W (db/s)')
    ylim([pmin,pmax])
    set(gca,'yticklabel',[])
    set(gca,'YDir','reverse')
    grid('on')

    ax4=axes('Position',ax3.Position);   
    xlabel('Chl-Turb','color','g')
    set(ax3,'box','off')
    plot(BIN.chlorophyll, BIN.pressure,'g', 'parent' , ax4,'linewidth',1)
    hold on
    plot(BIN.turbidity, BIN.pressure,'r', 'parent' , ax4,'linewidth',1)
    set(ax4,'XAxisLocation','top',...
        'YAxisLocation','right',...
        'Color','none',...
        'XColor','g','YColor','k');
    yticklabels([])
    ylim([pmin,pmax])
    set(gca,'yticklabel',[])
    set(gca,'YDir','reverse')


    subplot(2,4,3)
    plot(BIN.eps_T1,BIN.pressure,'.-','linewidth',1,'markersize',4)
    hold on
    plot(BIN.eps_T2,BIN.pressure,'.-','linewidth',1,'markersize',4)
    plot(BIN.eps_S1,BIN.pressure,'.-','linewidth',1,'markersize',4)
    if isfield(DATA,'sh2')
        plot(BIN.eps_S2,BIN.pressure,'.-','linewidth',1,'markersize',4)
    end
    plot(0.5*(BIN.epsT1max+BIN.epsT2max),BIN.pressure,'-','linewidth',1,'markersize',4,'color',[0.5,0.5,0.5])
    plot(BIN.eps_T1(BIN.flag_T1==0),BIN.pressure(BIN.flag_T1==0),'ob','markersize',3) 
    plot(BIN.eps_T2(BIN.flag_T2==0),BIN.pressure(BIN.flag_T2==0),'or','markersize',3) 
    if isfield(DATA,'sh2')
        leg=legend('T01','T02','sh1','sh2','maxT','fontsize',7,'location','best');
    else
        leg=legend('T01','T02','shear','maxT','fontsize',7,'location','best');
    end
    leg.ItemTokenSize = [7,7]; 
    xlabel('\epsilon (m^2 s^{-3})')
    set(gca,'yticklabel',[])
    set(gca,'YDir','reverse')
    set(gca, 'xscale','log')
    xlim([1e-12,1e-4])
    ylim([pmin,pmax])
    set(gca,'xtick',10.^[-12:2:-4])
    grid('on')
    title(datestr(date))

    subplot(2,4,4)
    l1=plot(BIN.Xi_T1,BIN.pressure,'.-','linewidth',1,'markersize',4);
    hold on
    l2=plot(BIN.Xi_T2,BIN.pressure,'.-','linewidth',1,'markersize',4);
    plot(BIN.Xi_ST1,BIN.pressure,'--','linewidth',1,'color',l1.Color)
    plot(BIN.Xi_ST2,BIN.pressure,'--','linewidth',1,'color',l2.Color)
    leg=legend('T01_f','T02_f','fontsize',7,'location','best');
    leg.ItemTokenSize = [7,7]; 
    xlabel('\chi (K^2 s^{-1})')
    set(gca,'yticklabel',[])
    set(gca,'YDir','reverse')
    set(gca,'xscale','log')
    xlim([1e-12,1e-3])
    ylim([pmin,pmax])
    set(gca,'xtick',10.^[-11:2:-3])
    grid('on')

    subplot(2,4,5)
    plot(BIN.KOsbornCox_T1,BIN.pressure,'.-','linewidth',1,'markersize',4)
    hold on
    plot(BIN.KOsbornCox_T2,BIN.pressure,'.-','linewidth',1,'markersize',4)
    plot(BIN.KOsborn_S1,BIN.pressure,'.-','linewidth',1,'markersize',4)
    if isfield(DATA,'sh2')
        plot(BIN.KOsborn_S2,BIN.pressure,'.-','linewidth',1,'markersize',4)
    end
    xlabel('K (m^2 s^{-1})')
    ylabel('p (db)')
    set(gca,'YDir','reverse')
    set(gca,'xscale','log')
    xlim([1e-9,1e0])
    line([D,D],ylim(),'color',[0.5,0.5,0.5])
    if isfield(DATA,'sh2')
        leg=legend('T01 (O&C)','T02 (O&C)','sh1 (O)','sh2 (O)','molecular','fontsize',7,'location','best');
    else
        leg=legend('T01 (O&C)','T02 (O&C)','she (O)','molecular','fontsize',7,'location','best');
    end
    leg.ItemTokenSize = [7,7]; 
    ylim([pmin,pmax])
     set(gca,'xtick',10.^[-7:2:0])
    grid('on')

     subplot(2,4,6)
% $$$     l1=plot(BIN.kB_T1,BIN.pressure,'.-','linewidth',1,'markersize',4);
% $$$     hold on
% $$$     l2=plot(BIN.kB_T2,BIN.pressure,'.-','linewidth',1,'markersize',4);
% $$$     plot(BIN.KBSH,BIN.pressure,'.-','linewidth',1,'markersize',4)
% $$$     plot(BIN.kU_T1,BIN.pressure,'--','linewidth',1,'color',l1.Color())
% $$$     plot(BIN.kU_T2,BIN.pressure,'--','linewidth',1,'color',l2.Color())
% $$$     xlabel('K_B (m)')
% $$$     set(gca,'YDir','reverse')
% $$$     %set(gca,'xscale','log')
% $$$     grid('on')
% $$$     ylim([pmin,pmax])
    
    plot(BIN.LTuT1,BIN.pressure,'.-','linewidth',1,'markersize',4)
    hold on
    plot(BIN.LTuT2,BIN.pressure,'.-','linewidth',1,'markersize',4)
        
    if isfield(DATA,'sh2')
        plot(BIN.LO_S1,BIN.pressure,'.-','linewidth',1,'markersize',4)
        plot(BIN.LO_S2,BIN.pressure,'.-','linewidth',1,'markersize',4)
        leg=legend('L_T^{T1}','L_T^{T2}','L_O^{sh1}','L_O^{sh2}','fontsize',7,'location','best');
    else
        plot(BIN.LO_S1,BIN.pressure,'.-','linewidth',1,'markersize',4)
        leg=legend('L_T^{T1}','L_T^{T2}','L_O','fontsize',7,'location','best');
    end
    leg.ItemTokenSize = [7,7]; 
    xlabel('L_T, L_O (m)')
    set(gca,'YDir','reverse')
    set(gca,'xscale','log')
    grid('on')
    ylim([pmin,pmax])
    
    subplot(2,4,7)
    plot(BIN.MAD_T1, BIN.pressure,'.-','linewidth',1,'markersize',4)
    hold on
    plot(BIN.MAD_T2, BIN.pressure,'.-','linewidth',1,'markersize',4)
    plot(BIN.MAD_S1, BIN.pressure,'.-','linewidth',1,'markersize',4)
    if isfield(DATA,'sh2')
        plot(BIN.MAD_S2, BIN.pressure,'.-','linewidth',1,'markersize',4)
    end
    line(2*[min(BIN.MADc),min(BIN.MADc)],ylim, 'color','k')
    line([min(BIN.MADc),min(BIN.MADc)],ylim, 'color','k')
    xlabel('MAD')
    xlim([0,2])
    set(gca,'yticklabel',[])
    set(gca,'YDir','reverse')
    grid('on')
    ylim([pmin,pmax])
    
    subplot(2,4,8)
    plot(BIN.LR_T1, BIN.pressure,'.-','linewidth',1,'markersize',4)
    hold on
    plot(BIN.LR_T2, BIN.pressure,'.-','linewidth',1,'markersize',4)
    xlabel('Likelihood ratio')
    set(gca,'yticklabel',[])
    set(gca,'YDir','reverse')
    grid('on')
    ylim([pmin,pmax])
    
    saveas(gcf,[folder_out,'/results_',num2str(inp,'%02d'),info.prof_dir,'.png'])
    

    %hold on
    %semilogx(epsilon2,-pres)
    %semilogx(epsilonN,-pres)
