function [BIN,SLOW] = resolve_VMP_profile_shear(DATA, inp, info,dataf_name,PLOT)
    DATA
    despike_sh  = [ 8  0.5 0.04];
    despike_A = [8 0.5000 0.0400];
    if ~ isfield(info,'system')
        info.system = 'Oce';
    end
    if ~ isfield(info,'minvel_detect')
        info.minvel_detect = 0.05;
    end
    if ~ isfield(info,'mindur_detect')
        info.minvel_detect = 30;
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
        info.k_HP_cut = 2;
    end
    if ~ isfield(info,'fAA')
        info.fAA = 110;
    end
    if ~ isfield(info,'minKT')
        info.minKT = 1;
    end
    if ~ isfield(info,'Tmethod')
        info.Tmethod = 'B';
    end
    if ~ isfield(info,'Tspec')
        info.Tspec = 'B';
    end
    if ~ isfield(info,'noisep_T1')
        %info.noisep_T1 = [-10.24,-0.89,info.fAA];
        info.noisep_T1 = [-10.5,-0.6,info.fAA];
    else
        info.noisep_T1 = [info.noisep_T1, info.fAA];
    end
    if ~ isfield(info,'noisep_T2')
        info.noisep_T2 = [-10.5,-0.60,info.fAA];
        %info.noisep_T2 = [-10.08,-0.97,info.fAA];
    else
        info.noisep_T2 = [info.noisep_T2, info.fAA];
    end
    if ~isfield(info,'time_res')
        info.time_res = nan;
    end
    if ~isfield(info,'chl_cal')
        info.chl_cal = [0.063,10];
    end
    if ~isfield(info,'turb_cal')
        info.turb_cal = [0.057,5];
    end
    if ~ isfield(info,'kmax_fac')
        info.kmax_factor = 1/1.66;
    end
    

    %adds some more noise
    %info.noisep_T1 = info.noisep_T1 + [0.25,0,0];
    %info.noisep_T2 = info.noisep_T2 + [0.25,0,0];
    
    %constants
    %visco = 1e-6;
    D = 1.44e-7;
    
    %defines times
    time_fast0 = [0:1:length(DATA.P_fast)-1]/DATA.fs_fast;
    time_slow0= [0:1:length(DATA.P_slow)-1]/DATA.fs_slow;
    
    %gets profiles
    %gets profiles
    iPf0 = get_profile(DATA.P_fast,DATA.W_fast,0,info.minvel_detect,info.prof_dir,info.mindur_detect,DATA.fs_fast);
    iPs0 = get_profile(DATA.P_slow,DATA.W_slow,0,info.minvel_detect,info.prof_dir,info.mindur_detect,DATA.fs_slow);
    NP = size(iPf0);
    NP = NP(2);
    if inp<0 > inp>NP
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

    %sampling frequencies
    fss = DATA.fs_slow;
    fsf = DATA.fs_fast;  
    
    %date
    date = datenum(DATA.Year, DATA.Month, DATA.Day, DATA.Hour, DATA.Minute, DATA.Second );

    %gets fast response sensors
    timef = time_fast0(iipf);
    Pf = DATA.P_fast(iipf);
    T1f = DATA.T1_fast(iipf);
    T2f = DATA.T2_fast(iipf);
    sh1 = DATA.sh1(iipf);
    sh2 = DATA.sh2(iipf);
    Ax = DATA.Ax(iipf);
    Ay = DATA.Ay(iipf);
    AA = [Ax,Ay];
    Wf = DATA.W_fast(iipf);
    
    %pseudo shear (from sensors)
    %psh_x = Ax./Wf;
    %psh_y = Ay./Wf;
    %pshear = Ay./Wf;
    
    %from velocity
    %acc = nan(size(Wf));
    %acc(2:end-1) = (Wf(3:end)-Wf(1:end-2))*(fsf/2);
    %acc(1) = (Wf(2)-Wf(1))*fsf;
    %acc(end) = (Wf(end)-Wf(end-1))*fsf;
    %pshear = acc./Wf;



    %gets slow response sensors
    times = time_slow0(iips);
    Ps = DATA.P_slow(iips);
    T_SB = DATA.SBT(iips);
    C_SB = DATA.SBC(iips);
    
    
    
    %pressure correction
    if isfield(info,'pres_cor')
       Ps = Ps*info.pres_cor(1) + info.pres_cor(2); 
       Pf = Pf*info.pres_cor(1) + info.pres_cor(2);
    end
    
    %accurate T in fast sensors grid
    T_SB_fast = interp1(Ps,T_SB,Pf);
    i1 = find(isfinite(T_SB_fast),1,'first');
    if i1>1
        T_SB_fast(1:i1-1) = T_SB_fast(i1);
    end
    i2 = find(isfinite(T_SB_fast),1,'last');
    if i2<length(T_SB_fast)
        T_SB_fast(i2+1:end) = T_SB_fast(i2);
    end


    %sets maximum depth
    if ~ isfield(info,'pmax')
        info.pmax = round(max(Ps)-info.dpD/2)+1;
    end

    %calculates salinity and density
    %[rhos,Ss,depths]=rho_salinity_Geneva(T_SB,C_SB,Ps);
    S_SB = salinity(Ps,T_SB,C_SB);
    sgt = sigmat(T_SB, S_SB)+1000;
    mrho = cumsum(sgt)./[1:length(sgt)]';
    depth=10000*Ps./(mrho*9.81);
    
    if info.system == 'Lem'
        %density for Leman
        [sgt,S_SB,depth]=rho_salinity_Geneva(T_SB,C_SB*1000,Ps);
    elseif info.system == 'Zue'
        [sgt,S_SB,depth]=rho_salinity_Zurich(T_SB,C_SB*1000,Ps);
    end
    
    %slow data for output
    SLOW.date = date;
    SLOW.filename = dataf_name;
    SLOW.depth = depth;
    SLOW.pres = Ps;
    SLOW.T = T_SB;
    SLOW.C = C_SB;
    SLOW.S = S_SB;
    SLOW.sigmat = sgt;

    
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
    

    %binned temperature, salinity and density
    %defines the presure vector where to calculate
    pres = [info.pmin:info.dp:info.pmax];
    BIN.pres = pres;
    BIN.date = date;
    BIN.filename = dataf_name;
    BIN.depth = pres_av(Ps,depth,pres,info.dp,2.7);
    BIN.T=pres_av(Ps,T_SB,pres,info.dp,2.7);
    BIN.C = pres_av(Ps,C_SB,pres,info.dp,2.7);
    BIN.S = pres_av(Ps,S_SB,pres,info.dp,2.7);
    BIN.sigmat = pres_av(Ps,sgt,pres,info.dp,2.7);
    BIN.grT = mean_grad(Ps,T_SB,pres,info.dpD);
    BIN.N2 = -9.81*mean_grad(Ps,sort_rho,pres,info.dpD)/1000;
    BIN.LT = sqrt(pres_av(Ps,displ.^2,pres,info.dpD,0));
    BIN.LTuT1 = sqrt(pres_av(Pf,displuT1.^2,pres,info.dpD,0));
    BIN.LTuT2 = sqrt(pres_av(Pf,displuT2.^2,pres,info.dpD,0));
    
    try
        Chl = (DATA.FL(iips)-info.chl_cal(1))*info.chl_cal(2);
        Turb = (DATA.NTU(iips)-info.turb_cal(1))*info.turb_cal(2);
        SLOW.Chl = Chl;
        SLOW.Turb = Turb;
        BIN.Chl = pres_av(Ps,Chl,pres,info.dp,2.7);
        BIN.Turb = pres_av(Ps,Turb,pres,info.dp,2.7);
    catch
        BIN.Chl = nan(size(BIN.T));
        BIN.Turb = nan(size(BIN.T));
    end
    %matches the FP07 to the high accuracy sensors
    T1f = calibration_FP07(timef,T1f, times, T_SB,fsf,fss);
    T2f = calibration_FP07(timef,T2f, times, T_SB, fsf,fss);

    %plots raw
    figure(1)
    clf
    subplot(4,1,1)
    plot(DATA.P_slow)
    ylabel('p (db)')
    hold on
    plot(iips,Ps)
    set(gca,'xticklabel',[])
    set(gca,'xticklabel',[])
    subplot(4,1,2)
    plot( T1f)
    hold on
    plot( T2f)
    ylabel('T fast (°C)')
    set(gca,'xticklabel',[])
    subplot(4,1,3)
    plot( sh1)
    hold on
    plot( sh2)
    ylabel('shear (1/s)')
    set(gca,'xticklabel',[])
    subplot(4,1,4)
    plot(Wf)
    ylabel('w (m/s)')
    saveas(gcf,['profile_',dataf_name,'_p','.png'])    
    

    %filters shear and highpasses microstructure
    %makes strange things at the borders (avoid?)
    [sh1, ~, ~, ~ ] =  despike(sh1, despike_sh(1), despike_sh(2), fsf, round(despike_sh(3)*fsf));
    [sh2, ~, ~, ~ ] =  despike(sh2, despike_sh(1), despike_sh(2), fsf, round(despike_sh(3)*fsf));
    mW = (max(Pf)-min(Pf))/(max(timef)-min(timef));
    if info.k_HP_cut>0
        f_HP_cut = info.k_HP_cut*mW;
        [bh,ah] = butter(1, f_HP_cut/(fsf/2), 'high');
        
        sh1_hp = filter(bh, ah, sh1);
        sh1_hp = flipud(sh1_hp);
        sh1_hp = filter(bh, ah, sh1_hp);
        sh1_hp = flipud(sh1_hp);
        
        sh2_hp = filter(bh, ah, sh2);
        sh2_hp = flipud(sh2_hp);
        sh2_hp = filter(bh, ah, sh2_hp);
        sh2_hp = flipud(sh2_hp);
        
        T1f_hp = filter(bh, ah, T1f);
        T1f_hp = flipud(T1f_hp);
        T1f_hp = filter(bh, ah, T1f_hp);
        T1f_hp = flipud(T1f_hp);

        T2f_hp = filter(bh, ah, T2f);
        T2f_hp = flipud(T2f_hp);
        T2f_hp = filter(bh, ah, T2f_hp);
        T2f_hp = flipud(T2f_hp);
    else
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
    
    %defines output variables
    
    BIN.W = nan(1,length(pres));
    
    BIN.epsSH1 = nan(1,length(pres));
    BIN.KBSH1 = nan(1,length(pres));
    BIN.MADsh1 = nan(1,length(pres));
    BIN.MADcsh1 = nan(1,length(pres));
    BIN.fit_flag_sh1 = nan(1,length(pres));
    
    BIN.epsSH2 = nan(1,length(pres));
    BIN.KBSH2 = nan(1,length(pres));
    BIN.MADsh2 = nan(1,length(pres));
    BIN.MADcsh2 = nan(1,length(pres));
    BIN.fit_flag_sh2 = nan(1,length(pres));
    
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
    BIN.MADc1 = nan(1,length(pres));
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
    BIN.MADc2 = nan(1,length(pres));
    BIN.fit_flag_T2 = nan(1,length(pres));

    
    for i = 1:length(pres)
        jp = find(Pf>=pres(i)-info.dpD/2 & Pf<=pres(i)+info.dpD/2);
        visco = viscosity(BIN.T(i));
        if isfield(info,'Nfft')
           Nfft = info.Nfft;
        else
            Nfft = round(length(jp)/2)-1;
        end
        
        if isfield(info,'overlap')
           overlap = info.overlap;
        else
           overlap = round(Nfft/2);
        end

        if length(jp)<Nfft | length(jp)<256
            continue
        end
        BIN.W(i) = mean(abs(Wf(jp)));
        
        try
            [BIN.epsSH1(i), ~, BIN.fit_flag_sh1(i), BIN.MADsh1(i), BIN.MADcsh1(i)] = dis_spec_ODAS(sh1_hp(jp), AA(jp,:), fsf, fss, Wf(jp), T_SB_fast(jp),timef(jp), Pf(jp), visco,Nfft, PLOT);
            BIN.KBSH1(i)=1/(2*pi())*(BIN.epsSH1(i)/(visco*D^2))^(1/4); 
        end
        
        
        try
            [BIN.epsSH2(i), ~, BIN.fit_flag_sh2(i), BIN.MADsh2(i), BIN.MADcsh2(i)] = dis_spec_ODAS(sh2_hp(jp), AA(jp,:), fsf, fss, Wf(jp), T_SB_fast(jp),timef(jp), Pf(jp), visco,Nfft,  PLOT);
            BIN.KBSH2(i)=1/(2*pi())*(BIN.epsSH2(i)/(visco*D^2))^(1/4); 
        end
        
        %FP07 calculations
        try
            [BIN.Xiv1(i),BIN.Xif1(i),BIN.KB1(i),BIN.fit_flag_T1(i),BIN.sXif1(i), BIN.sKB1(i), BIN.MADf1(i),BIN.MADc1(i),BIN.LKH1(i), BIN.LKHratio1(i), BIN.maxK1(i)] =Xi_spec_no_shear(Pf(jp),T1f_hp(jp),info.minKT,info.fAA,mean(abs(Wf(jp))),info.noisep_T1,  Nfft, overlap,info.Tspec,info.Tmethod,info.time_res,PLOT);
            BIN.epsT1(i) = visco*D^2*(2*pi()*BIN.KB1(i))^4;
            BIN.epsT1max(i) = visco*D^2*(2*pi()*BIN.maxK1(i)*info.kmax_factor)^4;
        end
        
        try
            [BIN.Xiv2(i),BIN.Xif2(i),BIN.KB2(i),BIN.fit_flag_T2(i), BIN.sXif2(i), BIN.sKB2(i), BIN.MADf2(i),BIN.MADc2(i),BIN.LKH1(i), BIN.LKHratio2(i), BIN.maxK2(i)] =Xi_spec_no_shear(Pf(jp),T2f_hp(jp),info.minKT,info.fAA, mean(abs(Wf(jp))),info.noisep_T2,  Nfft, overlap,info.Tspec,info.Tmethod,info.time_res, PLOT);
            BIN.epsT2(i) = visco*D^2*(2*pi()*BIN.KB2(i))^4;
            BIN.epsT2max(i) = visco*D^2*(2*pi()*BIN.maxK2(i)*info.kmax_factor)^4;
        end
        
    end


    BIN.KOsb1 = 0.2*BIN.epsSH1.*(BIN.N2).^-1;
    BIN.KOsb2 = 0.2*BIN.epsSH2.*(BIN.N2).^-1;
    BIN.KTf1 = 0.5*BIN.Xif1.*(BIN.grT).^-2;
    BIN.KTf2 = 0.5*BIN.Xif2.*(BIN.grT).^-2;
    [BIN.Reb1,~,BIN.KBB1] = Bouffard_model(BIN.epsT1,BIN.N2,BIN.T,7);
    [BIN.Reb2,~,BIN.KBB2] = Bouffard_model(BIN.epsT2,BIN.N2,BIN.T,7);
    BIN.LO1 = (BIN.epsT1./BIN.N2.^(3/2)).^(0.5);
    BIN.LO2 = (BIN.epsT2./BIN.N2.^(3/2)).^(0.5);

    pmin = info.pmin - info.dpD/2;
    pmax = pres( find(isfinite(BIN.T),1,'last')) + info.dpD/2;
    
    %plots profile
    figure(2)
    clf
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperSize', [29 20]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 29 20]);

    ax1=subplot(2,4,1);
    plot(T_SB, Ps,'linewidth',1)
    xlabel('T (°C)')
    ylabel('p (db)')
    set(gca,'YDir','reverse')
    ylim([pmin,pmax])
    ax2=axes('Position',get(ax1,'Position'));
    set(ax2,'box','off')
    plot(S_SB, Ps,'r', 'parent' , ax2,'linewidth',1)
    set(ax2,'XAxisLocation','top',...
            'YAxisLocation','right',...
            'Color','none',...
            'XColor','r','YColor','k');
    yticklabels([])
    ylim([pmin,pmax])
    xlabel('Salinity')
    set(gca,'yticklabel',[])
    set(gca,'YDir','reverse')
        
    ax3=subplot(2,4,2)
    plot(BIN.W,BIN.pres,'-','linewidth',1,'markersize',4)
    hold on
    plot(Wf, Pf, 'color','r')
    xlabel('W (db/s)')
    ylim([pmin,pmax])
    set(gca,'yticklabel',[])
    set(gca,'YDir','reverse')
    grid('on')
    
    ax4=axes('Position',get(ax3,'Position'));
    set(ax4,'box','off')
    plot(BIN.Chl, BIN.pres,'g', 'parent' , ax4,'linewidth',1)
    hold on
    plot(BIN.Turb, BIN.pres,'r', 'parent' , ax4,'linewidth',1)
    set(ax4,'XAxisLocation','top',...
            'YAxisLocation','right',...
            'Color','none',...
            'XColor','g','YColor','k');
    yticklabels([])
    ylim([pmin,pmax])
    xlabel('Chl-Turb','color','g')
    set(gca,'yticklabel',[])
    set(gca,'YDir','reverse')
    
    subplot(2,4,3)
    plot(BIN.epsT1,BIN.pres,'-','linewidth',1,'markersize',4)
    hold on
    plot(BIN.epsT2,BIN.pres,'-','linewidth',1,'markersize',4)
    plot(BIN.epsSH1,BIN.pres,'linewidth',1)
    plot(BIN.epsSH2,BIN.pres,'linewidth',1)
    plot(0.5*(BIN.epsT1max+BIN.epsT2max),BIN.pres,'-','linewidth',1,'markersize',4,'color',[0.5,0.5,0.5])
    legend('T01','T02','sh1','sh2','maxT','fontsize',7,'location','southeast')
    xlabel('\epsilon (m^2 s^{-3})')
    set(gca,'yticklabel',[])
    set(gca,'YDir','reverse')
    set(gca, 'xscale','log')
    xlim([1e-12,1e-4])
    ylim([pmin,pmax])
    set(gca,'xtick',10.^[-12:2:-4])
    grid('on')
    
    subplot(2,4,4)
    l1=plot(BIN.Xif1,BIN.pres,'-','linewidth',1,'markersize',4);
    hold on
    l2=plot(BIN.Xif2,BIN.pres,'-','linewidth',1,'markersize',4);
    plot(BIN.Xiv1,BIN.pres,'--','linewidth',1,'color',l1.Color())
    plot(BIN.Xiv2,BIN.pres,'--','linewidth',1,'color',l2.Color())
    legend('T01_f','T02_f','fontsize',7,'location','southeast')
    xlabel('\chi (K^2 s^{-1})')
    set(gca,'yticklabel',[])
    set(gca,'YDir','reverse')
    set(gca,'xscale','log')
    xlim([1e-11,1e-3])
    ylim([pmin,pmax])
    set(gca,'xtick',10.^[-11:2:-3])
    grid('on')
    
    subplot(2,4,5)
    plot(BIN.KTf1,BIN.pres,'-','linewidth',1,'markersize',4)
    hold on
    plot(BIN.KTf2,BIN.pres,'-','linewidth',1,'markersize',4)
    plot(BIN.KOsb1,BIN.pres,'-','linewidth',1,'markersize',4)
    plot(BIN.KOsb2,BIN.pres,'-','linewidth',1,'markersize',4)
    legend('T01 (O&C)','T02 (O&C)','T01 (O)', 'T02 (O)','fontsize',7,'location','southeast')
    xlabel('K (m^2 s^{-1})')
    ylabel('p (db)')
    set(gca,'YDir','reverse')
    set(gca,'xscale','log')
    xlim([1e-9,1e0])
    ylim([pmin,pmax])
     set(gca,'xtick',10.^[-7:2:0])
    grid('on')
    
    subplot(2,4,6)
%     l1=plot(BIN.KB1,BIN.pres,'-','linewidth',1,'markersize',4);
%     hold on
%     l2=plot(BIN.KB2,BIN.pres,'-','linewidth',1,'markersize',4);
%     plot(BIN.maxK1,BIN.pres,'--','linewidth',1,'color',l1.Color())
%     plot(BIN.maxK2,BIN.pres,'--','linewidth',1,'color',l2.Color())
%     xlabel('K_B (m)')
%     set(gca,'YDir','reverse')
%     %set(gca,'xscale','log')
%     grid('on')
%     ylim([pmin,pmax])
    
     plot(BIN.LT,BIN.pres,'-','linewidth',1,'markersize',4)
     hold on
     plot(BIN.LTuT1,BIN.pres,'-','linewidth',1,'markersize',4)
     plot(BIN.LTuT2,BIN.pres,'-','linewidth',1,'markersize',4)
     plot(BIN.LO1,BIN.pres,'-','linewidth',1,'markersize',4)
     plot(BIN.LO2,BIN.pres,'-','linewidth',1,'markersize',4)
     legend('L_T','L_T^{uT1}','L_T^{uT2}','L_O_1','L_O_2','fontsize',7,'location','southeast')
     xlabel('L_T, L_O (m)')
     set(gca,'YDir','reverse')
     set(gca,'xscale','log')
     grid('on')
     ylim([pmin,pmax])
    
    subplot(2,4,7)
    plot(BIN.MADf1, BIN.pres,'-','linewidth',1,'markersize',4)
    hold on
    plot(BIN.MADf2, BIN.pres,'-','linewidth',1,'markersize',4)
    plot(BIN.MADsh1, BIN.pres,'-','linewidth',1,'markersize',4)
    plot(BIN.MADsh2, BIN.pres,'-','linewidth',1,'markersize',4)
    plot(BIN.MADc1, BIN.pres,'k-')
    plot(BIN.MADc1*2, BIN.pres,'k-')
    xlabel('MAD')
    xlim([0,2])
    set(gca,'yticklabel',[])
    set(gca,'YDir','reverse')
    grid('on')
    ylim([pmin,pmax])
    
    subplot(2,4,8)
    plot(BIN.LKHratio1, BIN.pres,'-','linewidth',1,'markersize',4)
    hold on
    plot(BIN.LKHratio2, BIN.pres,'-','linewidth',1,'markersize',4)
    line([2,2],ylim, 'color','k')
    xlabel('Likelihood ratio')
    set(gca,'yticklabel',[])
    set(gca,'YDir','reverse')
    grid('on')
    ylim([pmin,pmax])
    
    saveas(gcf,[dataf_name,'_P','.png'])
    

    %hold on
    %semilogx(epsilon2,-pres)
    %semilogx(epsilonN,-pres)