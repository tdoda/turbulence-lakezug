function [BIN,SLOW,FAST,AAxy] = resolve_microCTD_profile(DATA, inp, info,dataf_name,PLOT,folder_out)
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
    if inp<0 | inp>NP
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
    sh = DATA.sh1(iipf);
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
    
    if info.system == 'Lem'
        %density for Leman
        [sgt,Ss,depth]=rho_salinity_Geneva(T_JAC,Cmatch*1000,Ps);
        [sgtT1,SsT1,depthT1]=rho_salinity_Geneva(T1f,interp1(Ps,Cmatch*1000,Pf,'linear','extrap'),Pf);
        [sgtT2,SsT2,depthT2]=rho_salinity_Geneva(T2f,interp1(Ps,Cmatch*1000,Pf,'linear','extrap'),Pf);
    elseif info.system == 'Zue'
        [sgt,Ss,depth]=rho_salinity_Zurich(T_JAC,Cmatch*1000,Ps);
    end  
    
    %Fast for output
    FAST.date = date;
    FAST.filename = dataf_name;
    if ismember(info.prof_dir,'up')
        FAST.upward = 1;
    else
        FAST.upward = 0;
    end
    FAST.pres = Pf;
    FAST.depth = interp1(Ps,depth,Pf,'linear','extrap'); % Good compromise!
    FAST.sigmatT1 = sgtT1;
    FAST.sigmatT2 = sgtT2;
    FAST.T1 = T1f;
    FAST.T2 = T2f;
    FAST.gradT1 = gradT1f;
    FAST.gradT2 = gradT2f;    
    FAST.Chl = Chl;
    FAST.Turb = Turb;
    
    %slow data for output
    SLOW.date = date;
    SLOW.filename = dataf_name;
    if ismember(info.prof_dir,'up')
        SLOW.upward = 1;
    else
        SLOW.upward = 0;
    end
    SLOW.pres = Ps;
    SLOW.depth = depth;
    SLOW.T = T_JAC;
    SLOW.C = Cmatch;
    SLOW.S = Ss;
    SLOW.sigmat = sgt;
    SLOW.Incl_X = Incl_X;
    SLOW.Incl_Y = Incl_Y;
    
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

    BIN.date = date;
    BIN.filename = dataf_name;
    if ismember(info.prof_dir,'up')
        BIN.upward = 1;
    else
        BIN.upward = 0;
    end
    BIN.pres = pres;
    BIN.depth = pres_av(Ps,depth,pres,info.dp,2.7);
    BIN.T = pres_av(Ps,T_JAC,pres,info.dp,2.7);
    BIN.C = pres_av(Ps,C_JAC,pres,info.dp,2.7);
    BIN.S = pres_av(Ps,Ss,pres,info.dp,2.7);
    BIN.sigmat = pres_av(Ps,sgt,pres,info.dp,2.7);
    Thorpe=sqrt(pres_av(Ps,displ.^2,pres,info.dpD,0));
    %     BIN.grT = mean_grad(Ps,T_JAC,pres,info.dpD);
    for i=1:length(pres)
        % It is better to use FP07 because when the lake is homogeneous,
        % the resolution of JAC-T is too low, resulting in grT=1E-16
        BIN.grT1(i) = mean_grad(Pf,sort_uT1,pres(i),max(Thorpe(i),info.dpD));   % SEB: evaluated on the resorted T profile
        BIN.grT2(i) = mean_grad(Pf,sort_uT2,pres(i),max(Thorpe(i),info.dpD));   % SEB: evaluated on the resorted T profile
        BIN.grT12(i)=nanmean([BIN.grT1(i),  BIN.grT2(i)]);
        
        BIN.grT(i) = mean_grad(Ps,sort_T_JAC,pres(i),max(Thorpe(i),info.dpD));   % SEB: evaluated on the resorted T profile
        BIN.N2(i) = -9.81*mean_grad(Ps,sort_rho,pres(i),max(Thorpe(i),info.dpD))/1000;
    end
%     BIN.grT = mean_grad(Ps,sort_T_JAC,pres,info.dpD);   % SEB: evaluated on the resorted T profile    
%     BIN.N2 = -9.81*mean_grad(Ps,sort_rho,pres,info.dpD)/1000;
    BIN.Chl = pres_av(Pf,Chl,pres,info.dp,2.7);
    BIN.Turb = pres_av(Pf,Turb,pres,info.dp,2.7);
    BIN.LT = Thorpe;
%    BIN.LTc = sqrt(pres_av(Ps,Lc.^2,pres,info.dpD,0));
    BIN.LTuT1 = sqrt(pres_av(Pf,displuT1.^2,pres,info.dpD,0));
    BIN.LTuT2 = sqrt(pres_av(Pf,displuT2.^2,pres,info.dpD,0));
    BIN.acc = -mean_grad(Pf,Wf,pres,info.dpD);  % instrument acceleration

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
    plot(sh,Pf); hold on; 
    plot(get(gca,'xlim'),[pmaxplot pmaxplot],'--k')
    set(gca,'ydir','reverse');
    ylabel('p(db)'); xlabel('sh (1/s)'); ylim(pplot);
%     subplot(4,4,[6,10,14])
%     plot(AA(:,1),Pf); hold on
%     plot(AA(:,2),Pf); set(gca,'ydir','reverse');
%     xlabel('Acc (1/s)')
    subplot(4,4,[7,11,15]);
    plot(T1f,Pf)
    hold on
    plot(T2f,Pf); 
    plot(get(gca,'xlim'),[pmaxplot pmaxplot],'--k')
    set(gca,'ydir','reverse');
    xlabel('T fast (°C)'); ylim(pplot)
    set(gca,'yticklabel',[])
    ax1=subplot(4,4,[8,12,16]);     
    plot(Wf,Pf,'-k'); hold on;    
    plot(get(gca,'xlim'),[pmaxplot pmaxplot],'--k')
    set(gca,'ydir','reverse');
    xlabel('W (m/s)');yticklabels([]);ylim(pplot); xlim(sort(sign(mean(Wf))*[0.4 0.7]));
    ax2=axes('Position',ax1.Position,'XAxisLocation','top',...
        'YAxisLocation','right','color','none',...
        'xColor','r','yColor','k');   
    set(ax1,'box','off')    
    line(Incl_X,DATA.P_slow(iips),'color','r'); set(gca,'ydir','reverse');
    xlabel('Inclination [°]')
    yticklabels([]); xlim([-5 5]);

    %filters shear and highpasses microstructure
    [sh, ~, ~, ~ ] =  despike(sh, despike_sh(1), despike_sh(2), fsf, round(despike_sh(3)*fsf));
    
    mW = (max(Pf)-min(Pf))/(max(timef)-min(timef));
    if info.k_HP_cut>0
        f_HP_cut = info.k_HP_cut*mW;
        [bh,ah] = butter(1, f_HP_cut/(fsf/2), 'high');
        
        sh_hp = filter(bh, ah, sh);
        sh_hp = flipud(sh_hp);
        sh_hp = filter(bh, ah, sh_hp);
        sh_hp = flipud(sh_hp);
        
        T2f_hp=T2f;
        T1f_hp=T1f;
        
    else
        sh_hp = sh;
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
    
    %-- identify the piezo-accelerometer signal to be used for noise filtering
    PSDsh = csd_odas(detrend(sh_hp),detrend(sh_hp),1024,fsf,[],512,'linear');
    PSDA1 = csd_odas(detrend(AA(:,1)),detrend(AA(:,1)),1024,fsf,[],512,'linear');
    PSDA2 = csd_odas(detrend(AA(:,2)),detrend(AA(:,2)),1024,fsf,[],512,'linear');
    CSDshA1 = csd_odas(detrend(sh_hp),detrend(AA(:,1)),1024,fsf,[],512,'linear');
    [CSDshA2,fA] = csd_odas(detrend(sh_hp),detrend(AA(:,2)),1024,fsf,[],512,'linear');
    COHshA1 = abs(CSDshA1.^2./(PSDA1.*PSDsh));
    COHshA2 = abs(CSDshA2.^2./(PSDA2.*PSDsh));
    AAxy=AA; % store both signals
    if nanmean(COHshA1)>nanmean(COHshA2)
        AA = AA(:,1);
    else
        AA = AA(:,2);
    end
    
    % add despiked and HP signals to the plot
    subplot(4,4,[5,9,13])
    plot(sh,Pf)
%     plot(sh_hp,Pf)
    %     legend('sh','sh_{desp}','sh_{HP}')
    subplot(4,4,[6,10,14])
    plot(AA,Pf); hold on
    set(gca,'ydir','reverse');
    plot([-200 -200],[min(Pf) max(Pf)],'--k'); plot([200 200],[min(Pf) max(Pf)],'--k');    
    xlim([-500 500]);  xlabel('Acc (1/s)'); ylim(pplot);
    plot(get(gca,'xlim'),[pmaxplot pmaxplot],'--k')
    saveas(gcf,[folder_out,'/profile',num2str(inp,'%02d'),'.png'])
    
    
    %defines output variables
    BIN.maxgrT = nan(1,length(pres));
        
    BIN.epsSH = nan(1,length(pres));
    BIN.KBSH = nan(1,length(pres));
    BIN.W = nan(1,length(pres));
    BIN.MADsh = nan(1,length(pres));
    BIN.MADcsh = nan(1,length(pres));
    BIN.fit_flag_sh = nan(1,length(pres));
	BIN.K1sh = nan(1,length(pres));
    BIN.K3sh = nan(1,length(pres));
    
    BIN.Xic1 = nan(1,length(pres));
    BIN.Xif1 = nan(1,length(pres));
    BIN.KB1 = nan(1,length(pres));
    BIN.sXif1 = nan(1,length(pres));
    BIN.sKB1 = nan(1,length(pres));
    BIN.Xiv1 = nan(1,length(pres));
    BIN.maxK1 = nan(1,length(pres));
    BIN.epsT1 = nan(1,length(pres));
    BIN.epsT1max = nan(1,length(pres));
    BIN.MAD1 = nan(1,length(pres));
    BIN.MADf1 = nan(1,length(pres));
    BIN.LKH1 = nan(1,length(pres));
    BIN.LKHratio1 = nan(1,length(pres));
    BIN.MAD1 = nan(1,length(pres));
    BIN.MADc1 = nan(1,length(pres));
    BIN.rangeK1 = nan(1,length(pres));
    BIN.fit_flag_T1 = nan(1,length(pres));
    BIN.Kp1 = nan(1,length(pres));
    BIN.minK1 = nan(1,length(pres));
    
    BIN.Xic2 = nan(1,length(pres));
    BIN.Xif2 = nan(1,length(pres));
    BIN.KB2 = nan(1,length(pres));
    BIN.sXif2 = nan(1,length(pres));
    BIN.sKB2 = nan(1,length(pres));
    BIN.Xiv2 = nan(1,length(pres));
    BIN.maxK2 = nan(1,length(pres));
    BIN.epsT2 = nan(1,length(pres));
    BIN.epsT2max = nan(1,length(pres));
    BIN.MAD2 = nan(1,length(pres));
    BIN.MADf2 = nan(1,length(pres));
    BIN.LKH2 = nan(1,length(pres));
    BIN.LKHratio2 = nan(1,length(pres));
    BIN.MAD2 = nan(1,length(pres));
    BIN.MADc2 = nan(1,length(pres));
    BIN.rangeK2 = nan(1,length(pres));
    BIN.fit_flag_T2 = nan(1,length(pres));
    BIN.Kp2 = nan(1,length(pres));
    BIN.minK2 = nan(1,length(pres));
    BIN.flag_vibr = nan(1,length(pres));
    BIN.flag_incl = nan(1,length(pres));
    BIN.flag_acc = nan(1,length(pres));
    
%     for i = 1:n_pres; %length(pres)
%         jp = find(Pf>=pres(i)-info.dpD/2 & Pf<=pres(i)+info.dpD/2);
%         jps = find(Ps>=pres(i)-info.dpD/2 & Ps<=pres(i)+info.dpD/2);   % slow channel
%         
%         % SEB: added evaluation of "maximum" gradT for pyroelectric analysis
%         if ~isempty(jps) & length(jps)>=3
%             tmp = movavg(T_JAC(jps),'linear',3);
%             BIN.maxgrT(i) = max(-gradient(tmp)./gradient(Ps(jps)));
%         end
%         % SEB: added evaluation of D and viscosity (based on http://web.mit.edu/seawater/)
%         avT = nanmean(T_JAC(jps));
%         avS = nanmean(Ss(jps));
%         avsigmat = nanmean(sgt(jps));
%         cond = SW_Conductivity(avT,'C',avS,'ppt');
%         cp = SW_SpcHeat(avT,'C',avS,'ppt',0.101325,'MPa');   % The last argument is sat vapor pressure. It can kept constant as suggested in SW_Diffusivity.m
%         D = cond/(avsigmat*cp);
%     
%         visco = SW_Viscosity(avT,'C',avS,'ppt'); % dynamic viscosity mu
%         visco = visco/avsigmat; % kinematic viscosity nu
%         
%         if isfield(info,'num_fft')
%            Nfft = floor(length(jp)/((info.num_fft+1)/2));
%         else
%            Nfft = floor(length(jp)/2);
%         end
%         
%         if isfield(info,'overlap')
%            overlap = info.overlap;
%         else
%            overlap = round(Nfft/2);
%         end
% 
%         if length(jp)<Nfft | length(jp)<128
%             continue
%         end
%         
%         % SEB: calculate the profiling speed to be used below. I prefer this way
%         % than using Wf or Ws, since these are smoothed and for short bins
%         % may be critical. However, in case just define WW=mean(abs(Wf(jp)))
%         % Note that in the main script I removed the smmoothing, however I
%         % leave this definition to be safe.
%         WW = mean( abs( (Pf(jp(end)) - Pf(jp(1)))/(timef(jp(end)) - timef(jp(1))) ) );
%         
%         % flagging for too much vibrations, flag if more than 5% of
%         % vibration per bin is >200 or if standard deviation is >100
%         % (i.e. 2 sigma C.L. = 95%)
%         AAi = AA(jp);
%         AAidev = nanstd(AAi);
%         AAimean = abs(nanmean(AAi));
%         len_AAi = nansum(abs(AAi)>200);
%         if len_AAi/length(jp) > 0.05
%             BIN.flag_vibr(i) = 1;
%         elseif (AAidev + AAimean) > 100
%             BIN.flag_vibr(i) = 1;
%         else
%             BIN.flag_vibr(i) = 0;
%         end
%         
%         % flaggin for inclination, same as vibr for incl. > 2°
%         Incli = Incl_X(jps);
%         Inclidev = nanstd(Incli);
%         Inclimean = abs(nanmean(Incli));
%         len_Incli = sum(abs(Incli)>2);
%         if len_Incli/length(jps) > 0.05
%             BIN.flag_incl(i) = 1;
%         elseif (Inclidev + Inclimean) > 2
%             BIN.flag_incl(i) = 1;
%         else
%             BIN.flag_incl(i) = 0;
%         end
%         
%         %flagging for acceleration of instrument$
%         %flag if variation is >20% at 95% C.L.
%         Wfi = Wf(jp);
%         Wfidev = nanstd(Wfi);
%         Wfimean = abs(nanmean(Wfi));
%         len_WWi = nansum(abs(Wfi/Wfimean)>1.2)+nansum(abs(Wfi/Wfimean)<0.8);
%         if len_WWi/length(jp) > 0.05
%             BIN.flag_acc(i) = 1;
%         elseif (Wfidev/Wfimean) > 0.1
%             BIN.flag_acc(i) = 1;
%         else
%             BIN.flag_acc(i) = 0;
%         end
%         
%         
%         %% ODAS epsilon calculation
%         try
% %             [BIN.epsSH(i), BIN.W(i), BIN.fit_flag_sh(i), BIN.MADsh(i), BIN.MADcsh(i)] = dis_spec_ODAS(sh_hp(jp), AAxy(jp,:), fsf, fss, WW, T_JAC_fast(jp),timef(jp), Pf(jp), visco,Nfft,  PLOT, folder_out);
%         if strcmp(info.Nasmyth_spec,'ODAS')
%             [BIN.epsSH(i), BIN.W(i), BIN.fit_flag_sh(i), BIN.MADsh(i), BIN.MADcsh(i)] = ...
%                 dis_spec_ODAS(sh_hp(jp), AA(jp), fsf, fss, WW, T_JAC_fast(jp),timef(jp), Pf(jp), visco,Nfft, PLOT, folder_out, info.noise_corr);
%         elseif strcmp(info.Nasmyth_spec,'EPFL')
%             [BIN.epsSH(i), BIN.W(i), BIN.fit_flag_sh(i), BIN.MADsh(i), BIN.MADcsh(i), BIN.K1sh(i),BIN.K3sh(i)] = ...
%                 dis_spec_EPFL(Pf(jp),sh_hp(jp),AA(jp),0.1,14,info.fAA,visco,WW, Nfft, overlap, PLOT,folder_out, info.noise_corr,'sh_1');
%         end
%         
%             BIN.KBSH(i)=1/(2*pi())*(BIN.epsSH(i)/(visco*D^2))^(1/4); 
%         end
%         %% FP07 calculations
% %         try
% %             [BIN.Xiv1(i),BIN.Xic1(i),BIN.Xif1(i),BIN.KB1(i),BIN.fit_flag_T1(i),BIN.sXif1(i), BIN.sKB1(i), BIN.MAD1(i), BIN.MADf1(i),BIN.MADc1(i),BIN.LKH1(i), BIN.LKHratio1(i), BIN.maxK1(i)] =Xi_spec(Pf(jp),T1f_hp(jp),info.minKT,info.fAA,BIN.KBSH(i),WW,info.noisep_T1,  Nfft, overlap,info.Tspec,info.Tmethod,info.time_res, info.peak_rem_T1,PLOT,timef(jp));
% %             [BIN.Xiv1(i),BIN.Xic1(i),BIN.Xif1(i),BIN.KB1(i),BIN.fit_flag_T1(i),BIN.sXif1(i), BIN.sKB1(i), BIN.MAD1(i), BIN.MADf1(i),BIN.MADc1(i),BIN.LKH1(i), BIN.LKHratio1(i), BIN.maxK1(i)] =Xi_spec_gradT(Pf(jp),gradT1f(jp),info.minKT,info.fAA,BIN.KBSH(i),WW,info.noisep_T1,  Nfft, overlap,info.Tspec,info.Tmethod,info.time_res, info.peak_rem_T1,PLOT);
%             [BIN.Xiv1(i),BIN.Xic1(i),BIN.Xif1(i),BIN.KB1(i),BIN.fit_flag_T1(i),BIN.sXif1(i), BIN.sKB1(i), BIN.MAD1(i), BIN.MADf1(i),BIN.MADc1(i),BIN.LKH1(i), BIN.LKHratio1(i), BIN.maxK1(i), BIN.rangeK1(i), BIN.Kp1(i), BIN.minK1(i)] =...
%                 Xi_spec_gradT_ODAS(Pf(jp),gradT1f(jp),info.minKT,info.fAA,BIN.KBSH(i),WW,info.noisep_T1, Nfft, overlap,info.Tspec,info.q,info.Tmethod,info.time_res,info.time_corr,info.npoles,info.peak_rem_T1,PLOT,timef(jp),D,visco,T1_dT1,'T1_dT1',DATA.setupfilestr,info.int_range,folder_out);
%           
%             BIN.epsT1(i) = visco*D^2*(2*pi()*BIN.KB1(i))^4;
% %             BIN.epsT1max(i) = visco*D^2*(2*pi()*BIN.maxK1(i)*info.kmax_factor)^4;
%             BIN.epsT1max(i) = visco*D^2*(2*pi()*info.fAA/WW*info.kmax_factor)^4;
% %         end
%         
% %         try
% %             [BIN.Xiv2(i),BIN.Xic2(i),BIN.Xif2(i),BIN.KB2(i),BIN.fit_flag_T2(i), BIN.sXif2(i), BIN.sKB2(i), BIN.MAD2(i), BIN.MADf2(i),BIN.MADc2(i),BIN.LKH1(i), BIN.LKHratio2(i), BIN.maxK2(i)] =Xi_spec(Pf(jp),T2f_hp(jp),info.minKT,info.fAA,BIN.KBSH(i),WW,info.noisep_T2,  Nfft, overlap,info.Tspec,info.Tmethod,info.time_res, info.peak_rem_T2,PLOT,timef(jp));
% %             [BIN.Xiv2(i),BIN.Xic2(i),BIN.Xif2(i),BIN.KB2(i),BIN.fit_flag_T2(i), BIN.sXif2(i), BIN.sKB2(i), BIN.MAD2(i), BIN.MADf2(i),BIN.MADc2(i),BIN.LKH1(i), BIN.LKHratio2(i), BIN.maxK2(i)] =Xi_spec_gradT(Pf(jp),gradT2f(jp),info.minKT,info.fAA,BIN.KBSH(i),WW,info.noisep_T2,  Nfft, overlap,info.Tspec,info.Tmethod,info.time_res, info.peak_rem_T2,PLOT);
%             [BIN.Xiv2(i),BIN.Xic2(i),BIN.Xif2(i),BIN.KB2(i),BIN.fit_flag_T2(i), BIN.sXif2(i), BIN.sKB2(i), BIN.MAD2(i), BIN.MADf2(i),BIN.MADc2(i),BIN.LKH1(i), BIN.LKHratio2(i), BIN.maxK2(i), BIN.rangeK2(i), BIN.Kp2(i), BIN.minK2(i)] =...
%                 Xi_spec_gradT_ODAS(Pf(jp),gradT2f(jp),info.minKT,info.fAA,BIN.KBSH(i),WW,info.noisep_T2, Nfft, overlap,info.Tspec,info.q,info.Tmethod,info.time_res,info.time_corr,info.npoles,info.peak_rem_T2,PLOT,timef(jp),D,visco,T2_dT2,'T2_dT2',DATA.setupfilestr,info.int_range,folder_out);
%          
%             BIN.epsT2(i) = visco*D^2*(2*pi()*BIN.KB2(i))^4;
% %             BIN.epsT2max(i) = visco*D^2*(2*pi()*BIN.maxK2(i)*info.kmax_factor)^4;
%             BIN.epsT2max(i) = visco*D^2*(2*pi()*info.fAA/WW*info.kmax_factor)^4;
% %         end
%     end
%     
%     % Create txt file if flag(s) exist(s)
%     flagvib = nansum(BIN.flag_vibr);
%     flaginc = nansum(BIN.flag_incl);
%     flagacc = nansum(BIN.flag_acc);
%     if flagvib == 0 & flaginc == 0 & flagacc == 0
%         %no flags, do nothing
%     else
%         flagname = sprintf('flag_profile%.0f.txt', inp);
%         fido = fopen(strcat(folder_out,flagname),'wt');
%         if flagvib > 0
%             fprintf(fido,'flagged due to bad vibrations\n');
%         else
%         end
%         if flaginc > 0
%             fprintf(fido,'flagged due to bad inclination\n');
%         else
%         end
%         if flagacc > 0
%             fprintf(fido,'flagged due to bad speed of instrumet\n');
%         else
%         end
%         fclose(fido);
%     end
%     
%     %% plots
%     BIN.Krho = 0.2*BIN.epsSH.*(BIN.N2).^-1;
% %     BIN.KTf1 = 0.5*BIN.Xif1.*(BIN.grT).^-2;
% %     BIN.KTf2 = 0.5*BIN.Xif2.*(BIN.grT).^-2;
% %     BIN.KTc1 = 0.5*BIN.Xic1.*(BIN.grT).^-2;
% %     BIN.KTc2 = 0.5*BIN.Xic2.*(BIN.grT).^-2;
%     BIN.KTf1 = 0.5*BIN.Xif1.*(BIN.grT12).^-2;
%     BIN.KTf2 = 0.5*BIN.Xif2.*(BIN.grT12).^-2;
%     BIN.KTc1 = 0.5*BIN.Xic1.*(BIN.grT12).^-2;
%     BIN.KTc2 = 0.5*BIN.Xic2.*(BIN.grT12).^-2;
%     BIN.LO = (BIN.epsSH./BIN.N2.^(3/2)).^(0.5);
% 
%     pmin = info.pmin - info.dpD/2;
%     pmax = pres( find(isfinite(BIN.T),1,'last')) + info.dpD/2;
%     
%     %% plots profile
%     figure(2)
%     clf
%     set(gcf, 'PaperUnits', 'centimeters');
%     set(gcf, 'PaperSize', [29 20]);
%     set(gcf, 'PaperPositionMode', 'manual');
%     set(gcf, 'PaperPosition', [0 0 29 20]);
% 
%     ax1=subplot(2,4,1);
%     plot(T_JAC, Ps,'-k','linewidth',1)
%     xlabel('T (°C)')
%     ylabel('p (db)')
%     set(gca,'YDir','reverse')
%     ylim([pmin,pmax])
%     ax2=axes('Position',ax1.Position);
%     set(ax1,'box','off')
%     plot(Ss, Ps,'r', 'parent' , ax2,'linewidth',1)
%     set(ax2,'XAxisLocation','top',...
%             'YAxisLocation','right',...
%             'Color','none',...
%             'XColor','r','YColor','k');
%     yticklabels([])
%     ylim([pmin,pmax])
%     xlabel('Salinity')
%     set(gca,'yticklabel',[])
%     set(gca,'YDir','reverse')
%         
%     ax3=subplot(2,4,2);
%     plot(BIN.W,BIN.pres,'.-k','linewidth',1,'markersize',4)
%     xlabel('W (db/s)')
%     ylim([pmin,pmax])
%     set(gca,'yticklabel',[])
%     set(gca,'YDir','reverse')
%     grid('on')
% 
%     ax4=axes('Position',ax3.Position);   
%     xlabel('Chl-Turb','color','g')
%     set(ax3,'box','off')
%     plot(BIN.Chl, BIN.pres,'g', 'parent' , ax4,'linewidth',1)
%     hold on
%     plot(BIN.Turb, BIN.pres,'r', 'parent' , ax4,'linewidth',1)
%     set(ax4,'XAxisLocation','top',...
%         'YAxisLocation','right',...
%         'Color','none',...
%         'XColor','g','YColor','k');
%     yticklabels([])
%     ylim([pmin,pmax])
%     set(gca,'yticklabel',[])
%     set(gca,'YDir','reverse')
%     
%     
%     subplot(2,4,3)
%     plot(BIN.epsT1,BIN.pres,'.-','linewidth',1,'markersize',4)
%     hold on
%     plot(BIN.epsT2,BIN.pres,'.-','linewidth',1,'markersize',4)
%     plot(BIN.epsSH,BIN.pres,'.-','linewidth',1,'markersize',4)
%     plot(0.5*(BIN.epsT1max+BIN.epsT2max),BIN.pres,'-','linewidth',1,'markersize',4,'color',[0.5,0.5,0.5])
%     plot(BIN.epsT1(BIN.fit_flag_T1==1),BIN.pres(BIN.fit_flag_T1==1),'ob','markersize',3) 
%     plot(BIN.epsT2(BIN.fit_flag_T2==1),BIN.pres(BIN.fit_flag_T2==1),'or','markersize',3) 
%     leg=legend('T01','T02','shear','maxT','fontsize',7,'location','best');
%     leg.ItemTokenSize = [7,7]; 
%     xlabel('\epsilon (m^2 s^{-3})')
%     set(gca,'yticklabel',[])
%     set(gca,'YDir','reverse')
%     set(gca, 'xscale','log')
%     xlim([1e-12,1e-4])
%     ylim([pmin,pmax])
%     set(gca,'xtick',10.^[-12:2:-4])
%     grid('on')
%     title(datestr(date))
%     
%     subplot(2,4,4)
%     l1=plot(BIN.Xif1,BIN.pres,'.-','linewidth',1,'markersize',4);
%     hold on
%     l2=plot(BIN.Xif2,BIN.pres,'.-','linewidth',1,'markersize',4);
%     plot(BIN.Xic1,BIN.pres,'--','linewidth',1,'color',l1.Color())
%     plot(BIN.Xic2,BIN.pres,'--','linewidth',1,'color',l2.Color())
%     leg=legend('T01_f','T02_f','fontsize',7,'location','best');
%     leg.ItemTokenSize = [7,7]; 
%     xlabel('\chi (K^2 s^{-1})')
%     set(gca,'yticklabel',[])
%     set(gca,'YDir','reverse')
%     set(gca,'xscale','log')
%     xlim([1e-12,1e-3])
%     ylim([pmin,pmax])
%     set(gca,'xtick',10.^[-11:2:-3])
%     grid('on')
%     
%     subplot(2,4,5)
%     plot(BIN.KTf1,BIN.pres,'.-','linewidth',1,'markersize',4)
%     hold on
%     plot(BIN.KTf2,BIN.pres,'.-','linewidth',1,'markersize',4)
%     plot(BIN.Krho,BIN.pres,'.-','linewidth',1,'markersize',4)
%     xlabel('K (m^2 s^{-1})')
%     ylabel('p (db)')
%     set(gca,'YDir','reverse')
%     set(gca,'xscale','log')
%     xlim([1e-9,1e0])
%     line([D,D],ylim(),'color',[0.5,0.5,0.5])
%     leg=legend('T01 (O&C)','T02 (O&C)','she (O)','molecular','fontsize',7,'location','best');
%     leg.ItemTokenSize = [7,7]; 
%     ylim([pmin,pmax])
%      set(gca,'xtick',10.^[-7:2:0])
%     grid('on')
%     
%      subplot(2,4,6)
% % $$$     l1=plot(BIN.KB1,BIN.pres,'.-','linewidth',1,'markersize',4);
% % $$$     hold on
% % $$$     l2=plot(BIN.KB2,BIN.pres,'.-','linewidth',1,'markersize',4);
% % $$$     plot(BIN.KBSH,BIN.pres,'.-','linewidth',1,'markersize',4)
% % $$$     plot(BIN.maxK1,BIN.pres,'--','linewidth',1,'color',l1.Color())
% % $$$     plot(BIN.maxK2,BIN.pres,'--','linewidth',1,'color',l2.Color())
% % $$$     xlabel('K_B (m)')
% % $$$     set(gca,'YDir','reverse')
% % $$$     %set(gca,'xscale','log')
% % $$$     grid('on')
% % $$$     ylim([pmin,pmax])
%     
%     plot(BIN.LTuT1,BIN.pres,'.-','linewidth',1,'markersize',4)
%     hold on
%     plot(BIN.LTuT2,BIN.pres,'.-','linewidth',1,'markersize',4)
%     plot(BIN.LO,BIN.pres,'.-','linewidth',1,'markersize',4)
%     leg=legend('L_T^{T1}','L_T^{T2}','L_O','fontsize',7,'location','best');
%     leg.ItemTokenSize = [7,7]; 
%     xlabel('L_T, L_O (m)')
%     set(gca,'YDir','reverse')
%     set(gca,'xscale','log')
%     grid('on')
%     ylim([pmin,pmax])
%     
%     subplot(2,4,7)
%     plot(BIN.MADf1, BIN.pres,'.-','linewidth',1,'markersize',4)
%     hold on
%     plot(BIN.MADf2, BIN.pres,'.-','linewidth',1,'markersize',4)
%     plot(BIN.MADsh, BIN.pres,'.-','linewidth',1,'markersize',4)
%     line(2*[min(BIN.MADc1),min(BIN.MADc1)],ylim, 'color','k')
%     line([min(BIN.MADc1),min(BIN.MADc1)],ylim, 'color','k')
%     xlabel('MAD')
%     xlim([0,2])
%     set(gca,'yticklabel',[])
%     set(gca,'YDir','reverse')
%     grid('on')
%     ylim([pmin,pmax])
%     
%     subplot(2,4,8)
%     plot(BIN.LKHratio1, BIN.pres,'.-','linewidth',1,'markersize',4)
%     hold on
%     plot(BIN.LKHratio2, BIN.pres,'.-','linewidth',1,'markersize',4)
%     xlabel('Likelihood ratio')
%     set(gca,'yticklabel',[])
%     set(gca,'YDir','reverse')
%     grid('on')
%     ylim([pmin,pmax])
%     
%     saveas(gcf,[folder_out,'/results_',num2str(inp,'%02d'),info.prof_dir,'.png'])
    

    %hold on
    %semilogx(epsilon2,-pres)
    %semilogx(epsilonN,-pres)