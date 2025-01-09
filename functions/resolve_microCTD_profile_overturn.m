function [BIN,SLOW,FAST] = resolve_microCTD_profile(DATA, inp, info,dataf_name,PLOT,folder_out)
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
    if ~ isfield(info,'mindp')
        info.mindp = 1;
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
    iPf0 = get_profile(DATA.P_fast,DATA.W_fast,0,info.minvel_detect,info.prof_dir,info.mindur_detect,DATA.fs_fast);
    iPs0 = get_profile(DATA.P_slow,DATA.W_slow,0,info.minvel_detect,info.prof_dir,info.mindur_detect,DATA.fs_slow);
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
    
    %calculates displacements for thorpe length
    if ismember(info.prof_dir,'down')
        [sort_rho, isd] = sort(sgt);
        displ = Ps - Ps(isd);
        
        [sort_uT1, ist1] = sort(T1f,'descend');
        displuT1 = Pf - Pf(ist1);
        
        [sort_uT2, ist2] = sort(T2f,'descend');
        displuT2 = Pf - Pf(ist2);
    else
        [sort_rho, isd] = sort(sgt,'descend');
        displ = Ps - Ps(isd);
        
        [sort_uT1, ist1] = sort(T1f,'ascend');
        displuT1 = Pf - Pf(ist1);
        
        [sort_uT2, ist2] = sort(T2f,'ascend');
        displuT2 = Pf - Pf(ist2);      
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
    
    if(1)
    % Segmentation Imberger and Ivey 1991
    info_seg.window=ceil(DATA.fs_fast/mean(Ws)*0.1);   % n counts in 10 cm relative to the mean profiling speed
    info_seg.poly=4;
    info_seg.spike_minsep=info_seg.window*2; info_seg.spike_minlength=info_seg.window*0.1;
    info_seg.segment_minsize=ceil(DATA.fs_fast/mean(Ws)); info_seg.segment_maxsize=ceil(DATA.fs_fast/mean(Ws)*5); 
    
    scale = 2*DATA.fs_fast/info_seg.window;    
    nscan=length(iipf);
    scanend = nscan-info_seg.window;
    dd_dt=zeros(nscan,3);
    
    gradT1f_hp = filter(bh, ah, gradT1f); gradT1f_hp = flipud(gradT1f_hp); gradT1f_hp = filter(bh, ah, gradT1f_hp); gradT1f_hp = flipud(gradT1f_hp);
    gradT2f_hp = filter(bh, ah, gradT2f); gradT2f_hp = flipud(gradT2f_hp); gradT2f_hp = filter(bh, ah, gradT2f_hp); gradT2f_hp = flipud(gradT2f_hp);
            
    NR = info_seg.window;
    NRprime = NR - info_seg.poly;
    for iscan=info_seg.window+1:scanend   % faster than calling it with a for loop and a fradX matrix
        dd_dt(iscan,1)=scale*calc_AR_distance(gradT1f_hp,iscan, NR, NRprime, info_seg.poly);
    end    
    for iscan=info_seg.window+1:scanend
        dd_dt(iscan,2)=scale*calc_AR_distance(gradT2f_hp,iscan, NR, NRprime, info_seg.poly);
    end
    for iscan=info_seg.window+1:scanend
        dd_dt(iscan,3)=scale*calc_AR_distance(sh_hp,iscan, NR, NRprime, info_seg.poly);
    end
    
    gradX=[gradT1f_hp, gradT2f_hp,sh_hp];
    figure
    subplot(1,5,1)
    plot(T_JAC,Ps); set(gca,'ydir','reverse'); xlabel('T [°C]'); ylabel('Depth [mbar]'); hold on
    subplot(1,5,2)
    plot(Lc, Ps); set(gca,'ydir','reverse'); xlabel('Thorpe displ. cent. [m]'); ylabel('Depth [mbar]'); hold on   
    for igrad=3:-1:1
        subplot(1,5,igrad+2)
        semilogx(dd_dt(:,igrad),Pf); hold on; set(gca,'ydir','reverse'); xlabel('dd_dt [°C]','Interpreter', 'none'); ylabel('Depth [mbar]')
        info_seg.spike_threshold=prctile(dd_dt(:,igrad),99); 
        xlim([1 1e5])
        [over_start,over_stop,n_segments]=find_segments(dd_dt(:,igrad),info_seg);
        for i_segments=1:n_segments
            idx = over_start(i_segments):over_stop(i_segments);
            pres(i_segments) = mean([Pf(over_start(i_segments)),Pf(over_stop(i_segments))]);
            Deltapres(i_segments) =  Pf(over_stop(i_segments)) - Pf(over_start(i_segments));
            x=1; y=Pf(over_start(i_segments));
            w= info_seg.spike_threshold; h=abs(Pf(over_stop(i_segments))-Pf(over_start(i_segments)));
            patch([x w w x],[y+h, y+h, y, y],[(igrad-1)/2 1/igrad (igrad-1)/2],'FaceAlpha',.3)
        end
    end
    saveas(gcf,[folder_out,'/segments',num2str(inp,'%02d'),'.png'])
    else
    % identify overturn events
    cumdispl=cumsum(displ);
    idx0=find(abs(cumdispl)<1e-8);
    if idx0(1)~=1
        idx0=[0; idx0];
    end
    breaks=find(diff(idx0)>1); 
    over_start=idx0(breaks)+1;
    over_end=idx0(breaks+1);
    % combine overturns if too short
    Deltapres_tmp=Ps(over_end) - Ps(over_start);
    idx_short=find(Deltapres_tmp< info.mindp); idx_long=find(Deltapres_tmp>= info.mindp);
    breaks=find(diff(idx_short)>1);
    breaks_=[0; breaks(1:end-1)];
    over_start_short=over_start(idx_short(breaks_+1));
    over_end_short=over_end(idx_short(breaks));
    % split combined overturns in segments of info.mindp m
    Deltapres_tmp=Ps(over_end_short) - Ps(over_start_short);
    over_start_split=[]; over_end_split=[];
    for i_split=1:length(Deltapres_tmp)
        n_split=floor(Deltapres_tmp(i_split)/(5*info.mindp));
        if n_split==1
            over_start_split=[over_start_split; over_start_short(i_split)];
            over_end_split=[over_end_split; over_end_short(i_split)];
        elseif n_split>1
            n_tot=over_end_short(i_split)-over_start_short(i_split)+1;
            n_seg=floor(n_tot/n_split);
            for ii_split=1:n_split
                over_start_split=[over_start_split; over_start_short(i_split)+(ii_split-1)*n_seg];
                if ii_split==n_split
                    over_end_split=[over_end_split; over_end_short(i_split)];
                else
                    over_end_split=[over_end_split; over_start_short(i_split)+ii_split*n_seg-1];
                end
            end
        end
    end
        
    over_start=sort([over_start_split; over_start(idx_long)]);
    over_end=sort([over_end_split; over_end(idx_long)]);
    Deltapres_tmp=Ps(over_end) - Ps(over_start);
    idx_short=find(Deltapres_tmp<info.mindp*0.95); % some confidence because falling rate is not const
    over_start(idx_short)=[];
    over_end(idx_short)=[];
    
    LT=zeros(1,length(sgt));
    n_overturns=length(over_start);
    for i_overturns=1:n_overturns
        idx = over_start(i_overturns):over_end(i_overturns);
        LT(idx)=sqrt(mean(displ(idx).^2));
        pres(i_overturns) = mean([Ps(over_start(i_overturns)),Ps(over_end(i_overturns))]);
        Deltapres(i_overturns) =  Ps(over_end(i_overturns)) - Ps(over_start(i_overturns));
    end
      figure
    subplot(1,4,1)
    plot(T_JAC,Ps); set(gca,'ydir','reverse'); xlabel('T [°C]'); ylabel('Depth [mbar]'); hold on
    subplot(1,4,3)
    plot(Lc, Ps); set(gca,'ydir','reverse'); xlabel('Thorpe displ. cent. [m]'); ylabel('Depth [mbar]'); hold on   
     subplot(1,4,2)
    plot(sgt, Ps); set(gca,'ydir','reverse'); xlabel('Density [kg/m3]'); ylabel('Depth [mbar]'); hold on   
        subplot(1,4,4);  set(gca,'ydir','reverse');
           for i_overturns=1:n_overturns
            idx = over_start(i_overturns):over_end(i_overturns);
            pres(i_overturns) = mean([Ps(over_start(i_overturns)),Ps(over_end(i_overturns))]);
            Deltapres(i_overturns) =  Ps(over_end(i_overturns)) - Ps(over_start(i_overturns));
            x=0; y=Ps(over_start(i_overturns));
            w= 1; h=abs(Ps(over_end(i_overturns))-Ps(over_start(i_overturns)));
            patch([x w w x],[y+h, y+h, y, y],[0 1 0],'FaceAlpha',.3)
        end
    saveas(gcf,[folder_out,'/segments',num2str(inp,'%02d'),'.png'])
    
    
    end

    %defines the pressure vector where to calculate
    if ismember(info.prof_dir,'up') % For upward profiles, start the profile after the acceleration
        idx75=find(Ps<0.75*info.pmax,1,'first'); % Note: Ps is decreasing (upward profile). 
        [~,idx]=min(Ws(1:idx75)); % Focus just on the first 25% (due to possible wave-induced oscillations)
        n_pres = find(pres<=Ps(idx),1,'last');
    else
        n_pres=length(pres);
    end
    pmaxplot=pres(n_pres);

   %plots raw 
    pplot=[info.pmin info.pmax];
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
   subplot(4,4,[6,10,14])
    plot(AA,Pf); hold on
    set(gca,'ydir','reverse');
    plot([-200 -200],[min(Pf) max(Pf)],'--k'); plot([200 200],[min(Pf) max(Pf)],'--k');    
    xlim([-500 500]);  xlabel('Acc (1/s)'); ylim(pplot);
    plot(get(gca,'xlim'),[pmaxplot pmaxplot],'--k')
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
    xlabel('W (m/s)');yticklabels([]);ylim(pplot);
    ax2=axes('Position',ax1.Position,'XAxisLocation','top',...
        'YAxisLocation','right','color','none',...
        'xColor','r','yColor','k');   
    set(ax1,'box','off')    
    line(Incl_X,DATA.P_slow(iips),'color','r'); set(gca,'ydir','reverse');
    xlabel('Inclination [°]')
    yticklabels([]); xlim([-5 5]);
    saveas(gcf,[folder_out,'/profile',num2str(inp,'%02d'),'.png'])

    %binned temperature, salinity and density
    BIN.date = date;
    BIN.filename = dataf_name;
    if ismember(info.prof_dir,'up')
        BIN.upward = 1;
    else
        BIN.upward = 0;
    end
    BIN.pres = pres;
    BIN.depth = pres_av_overturns(Ps,depth,pres,Deltapres,2.7);
    BIN.T = pres_av_overturns(Ps,T_JAC,pres,Deltapres,2.7);
    BIN.C = pres_av_overturns(Ps,C_JAC,pres,Deltapres,2.7);
    BIN.S = pres_av_overturns(Ps,Ss,pres,Deltapres,2.7);
    BIN.sigmat = pres_av_overturns(Ps,sgt,pres,Deltapres,2.7);
    BIN.grT = mean_grad_overturns(Ps,T_JAC,pres,Deltapres);
    BIN.N2 = -9.81*mean_grad_overturns(Ps,sort_rho,pres,Deltapres)/1000;
    BIN.Chl = pres_av_overturns(Pf,Chl,pres,Deltapres,2.7);
    BIN.Turb = pres_av_overturns(Pf,Turb,pres,Deltapres,2.7);
    BIN.LT = sqrt(pres_av_overturns(Ps,displ.^2,pres,Deltapres,0));
    BIN.LTc = sqrt(pres_av(Ps,Lc.^2,pres,info.dpD,0));
    BIN.LTuT1 = sqrt(pres_av_overturns(Pf,displuT1.^2,pres,Deltapres,0));
    BIN.LTuT2 = sqrt(pres_av_overturns(Pf,displuT2.^2,pres,Deltapres,0));
    BIN.acc = -mean_grad_overturns(Pf,Wf,pres,Deltapres);  % instrument acceleration  
    
    %defines output variables
    BIN.epsSH = nan(1,length(pres));
    BIN.KBSH = nan(1,length(pres));
    BIN.W = nan(1,length(pres));
    BIN.MADsh = nan(1,length(pres));
    BIN.MADcsh = nan(1,length(pres));
    BIN.fit_flag_sh = nan(1,length(pres));
    
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
    
    for i = 1:n_pres; %length(pres)
        jp = find(Pf>=pres(i)-Deltapres(i)/2 & Pf<=pres(i)+Deltapres(i)/2);
        jps = find(Ps>=pres(i)-Deltapres(i)/2 & Ps<=pres(i)+Deltapres(i)/2);   % slow channel
        
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
%         Nfft=min(Nfft,floor(length(jp)/ceil(length(jp)/2048)));
        
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
        WW = mean( abs( (Pf(jp(end)) - Pf(jp(1)))/(timef(jp(end)) - timef(jp(1))) ) );
        
        %% ODAS epsilon calculation
        try
%             [BIN.epsSH(i), BIN.W(i), BIN.fit_flag_sh(i), BIN.MADsh(i), BIN.MADcsh(i)] = dis_spec_ODAS(sh_hp(jp), AAxy(jp,:), fsf, fss, WW, T_JAC_fast(jp),timef(jp), Pf(jp), visco,Nfft,  PLOT, folder_out);
        if strcmp(info.Nasmyth_spec,'ODAS')
            [BIN.epsSH(i), BIN.W(i), BIN.fit_flag_sh(i), BIN.MADsh(i), BIN.MADcsh(i)] = ...
                dis_spec_ODAS(sh_hp(jp), AA(jp), fsf, fss, WW, T_JAC_fast(jp),timef(jp), Pf(jp), visco,Nfft, PLOT, folder_out, info.noise_corr);
        elseif strcmp(info.Nasmyth_spec,'EPFL')
            [BIN.epsSH(i), BIN.W(i), BIN.fit_flag_sh(i), BIN.MADsh(i), BIN.MADcsh(i)] = ...
                dis_spec_EPFL(Pf(jp),sh_hp(jp),AA(jp),0.1,14,30,visco,WW, Nfft, overlap, PLOT,folder_out, info.noise_corr,'sh_1');
        end
        
            BIN.KBSH(i)=1/(2*pi())*(BIN.epsSH(i)/(visco*D^2))^(1/4); 
        end
        %% FP07 calculations
        try
%             [BIN.Xiv1(i),BIN.Xic1(i),BIN.Xif1(i),BIN.KB1(i),BIN.fit_flag_T1(i),BIN.sXif1(i), BIN.sKB1(i), BIN.MAD1(i), BIN.MADf1(i),BIN.MADc1(i),BIN.LKH1(i), BIN.LKHratio1(i), BIN.maxK1(i)] =Xi_spec(Pf(jp),T1f_hp(jp),info.minKT,info.fAA,BIN.KBSH(i),WW,info.noisep_T1,  Nfft, overlap,info.Tspec,info.Tmethod,info.time_res, info.peak_rem_T1,PLOT,timef(jp));
%             [BIN.Xiv1(i),BIN.Xic1(i),BIN.Xif1(i),BIN.KB1(i),BIN.fit_flag_T1(i),BIN.sXif1(i), BIN.sKB1(i), BIN.MAD1(i), BIN.MADf1(i),BIN.MADc1(i),BIN.LKH1(i), BIN.LKHratio1(i), BIN.maxK1(i)] =Xi_spec_gradT(Pf(jp),gradT1f(jp),info.minKT,info.fAA,BIN.KBSH(i),WW,info.noisep_T1,  Nfft, overlap,info.Tspec,info.Tmethod,info.time_res, info.peak_rem_T1,PLOT);
            [BIN.Xiv1(i),BIN.Xic1(i),BIN.Xif1(i),BIN.KB1(i),BIN.fit_flag_T1(i),BIN.sXif1(i), BIN.sKB1(i), BIN.MAD1(i), BIN.MADf1(i),BIN.MADc1(i),BIN.LKH1(i), BIN.LKHratio1(i), BIN.maxK1(i), BIN.rangeK1(i)] =...
                Xi_spec_gradT_ODAS(Pf(jp),gradT1f(jp),info.minKT,info.fAA,BIN.KBSH(i),WW,info.noisep_T1, Nfft, overlap,info.Tspec,info.q,info.Tmethod,info.time_res,info.time_corr,info.npoles,info.peak_rem_T1,PLOT,timef(jp),D,visco,T1_dT1,'T1_dT1',DATA.setupfilestr,info.int_range,folder_out);
          
            BIN.epsT1(i) = visco*D^2*(2*pi()*BIN.KB1(i))^4;
            BIN.epsT1max(i) = visco*D^2*(2*pi()*BIN.maxK1(i)*info.kmax_factor)^4;
        end
        
        try
%             [BIN.Xiv2(i),BIN.Xic2(i),BIN.Xif2(i),BIN.KB2(i),BIN.fit_flag_T2(i), BIN.sXif2(i), BIN.sKB2(i), BIN.MAD2(i), BIN.MADf2(i),BIN.MADc2(i),BIN.LKH2(i), BIN.LKHratio2(i), BIN.maxK2(i)] =Xi_spec(Pf(jp),T2f_hp(jp),info.minKT,info.fAA,BIN.KBSH(i),WW,info.noisep_T2,  Nfft, overlap,info.Tspec,info.Tmethod,info.time_res, info.peak_rem_T2,PLOT,timef(jp));
%             [BIN.Xiv2(i),BIN.Xic2(i),BIN.Xif2(i),BIN.KB2(i),BIN.fit_flag_T2(i), BIN.sXif2(i), BIN.sKB2(i), BIN.MAD2(i), BIN.MADf2(i),BIN.MADc2(i),BIN.LKH2(i), BIN.LKHratio2(i), BIN.maxK2(i)] =Xi_spec_gradT(Pf(jp),gradT2f(jp),info.minKT,info.fAA,BIN.KBSH(i),WW,info.noisep_T2,  Nfft, overlap,info.Tspec,info.Tmethod,info.time_res, info.peak_rem_T2,PLOT);
            [BIN.Xiv2(i),BIN.Xic2(i),BIN.Xif2(i),BIN.KB2(i),BIN.fit_flag_T2(i), BIN.sXif2(i), BIN.sKB2(i), BIN.MAD2(i), BIN.MADf2(i),BIN.MADc2(i),BIN.LKH2(i), BIN.LKHratio2(i), BIN.maxK2(i), BIN.rangeK2(i)] =...
                Xi_spec_gradT_ODAS(Pf(jp),gradT2f(jp),info.minKT,info.fAA,BIN.KBSH(i),WW,info.noisep_T2, Nfft, overlap,info.Tspec,info.q,info.Tmethod,info.time_res,info.time_corr,info.npoles,info.peak_rem_T2,PLOT,timef(jp),D,visco,T2_dT2,'T2_dT2',DATA.setupfilestr,info.int_range,folder_out);
         
            BIN.epsT2(i) = visco*D^2*(2*pi()*BIN.KB2(i))^4;
            BIN.epsT2max(i) = visco*D^2*(2*pi()*BIN.maxK2(i)*info.kmax_factor)^4;
        end
    end
    
    % add extra rows in BIN to ensure that the length is the same for all profiles
    nmax_pres=100;
    if nmax_pres<n_pres
        error('nmax_pres<n_pres. Stop!')
    end
    fields_name = fieldnames(BIN);
    for ifields=4:length(fields_name)
        BIN.(fields_name{ifields})(n_pres+1:nmax_pres)=NaN;
    end        
    
    %% plots
    BIN.Krho = 0.2*BIN.epsSH.*(BIN.N2).^-1;
    BIN.KTf1 = 0.5*BIN.Xif1.*(BIN.grT).^-2;
    BIN.KTf2 = 0.5*BIN.Xif2.*(BIN.grT).^-2;
    BIN.KTc1 = 0.5*BIN.Xic1.*(BIN.grT).^-2;
    BIN.KTc2 = 0.5*BIN.Xic2.*(BIN.grT).^-2;
    BIN.LO = (BIN.epsSH./BIN.N2.^(3/2)).^(0.5);

    pmin = info.pmin - Deltapres(1)/2;
    pmax = info.pmax; % pres( find(isfinite(BIN.T),1,'last')) + Deltapres(end)/2;
    
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
    plot(BIN.W,BIN.pres,'.-k','linewidth',1,'markersize',4)
    xlabel('W (db/s)')
    ylim([pmin,pmax])
    set(gca,'yticklabel',[])
    set(gca,'YDir','reverse')
    grid('on')

    ax4=axes('Position',ax3.Position);   
    xlabel('Chl-Turb','color','g')
    set(ax3,'box','off')
    plot(BIN.Chl, BIN.pres,'g', 'parent' , ax4,'linewidth',1)
    hold on
    plot(BIN.Turb, BIN.pres,'r', 'parent' , ax4,'linewidth',1)
    set(ax4,'XAxisLocation','top',...
        'YAxisLocation','right',...
        'Color','none',...
        'XColor','g','YColor','k');
    yticklabels([])
    ylim([pmin,pmax])
    set(gca,'yticklabel',[])
    set(gca,'YDir','reverse')
    
    
    subplot(2,4,3)
    plot(BIN.epsT1,BIN.pres,'.-','linewidth',1,'markersize',4)
    hold on
    plot(BIN.epsT2,BIN.pres,'.-','linewidth',1,'markersize',4)
    plot(BIN.epsSH,BIN.pres,'.-','linewidth',1,'markersize',4)
    plot(0.5*(BIN.epsT1max+BIN.epsT2max),BIN.pres,'-','linewidth',1,'markersize',4,'color',[0.5,0.5,0.5])
    plot(BIN.epsT1(BIN.fit_flag_T1==1),BIN.pres(BIN.fit_flag_T1==1),'ob','markersize',3) 
    plot(BIN.epsT2(BIN.fit_flag_T2==1),BIN.pres(BIN.fit_flag_T2==1),'or','markersize',3)
    

%     for i_overturns=1:n_overturns
%            
%     x=1E-14; y=Ps(over_start(i_overturns)); 
%     w=(BIN.epsT1(i_overturns)); h=abs(Ps(over_end(i_overturns))-Ps(over_start(i_overturns)));
% %     rectangle('position',[x,y,w,h],'EdgeColor',[0,0.4470,0.7410],'FaceColor',[0,0.4470,0.7410]); hold on
%     patch([x w w x],[y+h, y+h, y, y],[0,0.4470,0.7410],'FaceAlpha',.3)
%     
%     
%     x=1E-14; y=Ps(over_start(i_overturns));
%     w=(BIN.epsT2(i_overturns)); h=abs(Ps(over_end(i_overturns))-Ps(over_start(i_overturns)));
% %     rectangle('position',[x,y,w,h],'EdgeColor',[0.8500    0.3250    0.0980],'FaceColor',[0.8500    0.3250    0.0980]); hold on
%         patch([x w w x],[y+h, y+h, y, y],[0.8500    0.3250    0.0980],'FaceAlpha',.3)
% 
%     x=1E-14; y=Ps(over_start(i_overturns));
%     w=(BIN.epsSH(i_overturns)); h=abs(Ps(over_end(i_overturns))-Ps(over_start(i_overturns)));
% %     rectangle('position',[x,y,w,h],'EdgeColor',[0.9290    0.6940    0.1250],'FaceColor',[0.9290    0.6940    0.1250]); hold on
%             patch([x w w x],[y+h, y+h, y, y],[0.9290    0.6940    0.1250],'FaceAlpha',.3)
% 
%     end

    leg=legend('T01','T02','shear','maxT','fontsize',7,'location','best');
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
    l1=plot(BIN.Xif1,BIN.pres,'.-','linewidth',1,'markersize',4);
    hold on
    l2=plot(BIN.Xif2,BIN.pres,'.-','linewidth',1,'markersize',4);
    plot(BIN.Xic1,BIN.pres,'--','linewidth',1,'color',l1.Color())
    plot(BIN.Xic2,BIN.pres,'--','linewidth',1,'color',l2.Color())
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
    plot(BIN.KTf1,BIN.pres,'.-','linewidth',1,'markersize',4)
    hold on
    plot(BIN.KTf2,BIN.pres,'.-','linewidth',1,'markersize',4)
    plot(BIN.Krho,BIN.pres,'.-','linewidth',1,'markersize',4)
    xlabel('K (m^2 s^{-1})')
    ylabel('p (db)')
    set(gca,'YDir','reverse')
    set(gca,'xscale','log')
    xlim([1e-9,1e0])
    line([D,D],ylim(),'color',[0.5,0.5,0.5])
    leg=legend('T01 (O&C)','T02 (O&C)','she (O)','molecular','fontsize',7,'location','best');
    leg.ItemTokenSize = [7,7]; 
    ylim([pmin,pmax])
     set(gca,'xtick',10.^[-7:2:0])
    grid('on')
    
     subplot(2,4,6)
% $$$     l1=plot(BIN.KB1,BIN.pres,'.-','linewidth',1,'markersize',4);
% $$$     hold on
% $$$     l2=plot(BIN.KB2,BIN.pres,'.-','linewidth',1,'markersize',4);
% $$$     plot(BIN.KBSH,BIN.pres,'.-','linewidth',1,'markersize',4)
% $$$     plot(BIN.maxK1,BIN.pres,'--','linewidth',1,'color',l1.Color())
% $$$     plot(BIN.maxK2,BIN.pres,'--','linewidth',1,'color',l2.Color())
% $$$     xlabel('K_B (m)')
% $$$     set(gca,'YDir','reverse')
% $$$     %set(gca,'xscale','log')
% $$$     grid('on')
% $$$     ylim([pmin,pmax])
    
    plot(BIN.LTuT1,BIN.pres,'.-','linewidth',1,'markersize',4)
    hold on
    plot(BIN.LTuT2,BIN.pres,'.-','linewidth',1,'markersize',4)
    plot(BIN.LO,BIN.pres,'.-','linewidth',1,'markersize',4)
    leg=legend('L_T^{T1}','L_T^{T2}','L_O','fontsize',7,'location','best');
    leg.ItemTokenSize = [7,7]; 
    xlabel('L_T, L_O (m)')
    set(gca,'YDir','reverse')
    set(gca,'xscale','log')
    grid('on')
    ylim([pmin,pmax])
    
    subplot(2,4,7)
    plot(BIN.MADf1, BIN.pres,'.-','linewidth',1,'markersize',4)
    hold on
    plot(BIN.MADf2, BIN.pres,'.-','linewidth',1,'markersize',4)
    plot(BIN.MADsh, BIN.pres,'.-','linewidth',1,'markersize',4)
    line(2*[min(BIN.MADc1),min(BIN.MADc1)],ylim, 'color','k')
    line([min(BIN.MADc1),min(BIN.MADc1)],ylim, 'color','k')
    xlabel('MAD')
    xlim([0,2])
    set(gca,'yticklabel',[])
    set(gca,'YDir','reverse')
    grid('on')
    ylim([pmin,pmax])
    
    subplot(2,4,8)
    plot(BIN.LKHratio1, BIN.pres,'.-','linewidth',1,'markersize',4)
    hold on
    plot(BIN.LKHratio2, BIN.pres,'.-','linewidth',1,'markersize',4)
    xlabel('Likelihood ratio')
    set(gca,'yticklabel',[])
    set(gca,'YDir','reverse')
    grid('on')
    ylim([pmin,pmax])
    
    saveas(gcf,[folder_out,'/results_',num2str(inp,'%02d'),info.prof_dir,'.png'])
    

    %hold on
    %semilogx(epsilon2,-pres)
    %semilogx(epsilonN,-pres)