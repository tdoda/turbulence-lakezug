function [BIN,SLOW] = resolve_microCTD_profile_uC(DATA, inp, info,dataf_name,PLOT)
    despike_sh  = [ 8  0.5 0.04];
    despike_A = [8 0.5000 0.0400];
    
    if ~ isfield(info,'minvel_detect')
        info.minvel_detect = 0.1;
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
    if ~ isfield(info,'HP_cut')
        info.HP_cut = 0.4;
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
    end
    if ~ isfield(info,'noisep_T2')
        info.noisep_T2 = [-10.08,-0.97,info.fAA];
    end

    %adds some more noise
    info.noisep_T1 = info.noisep_T1 + [0.25,0,0];
    info.noisep_T2 = info.noisep_T2 + [0.25,0,0];
        
    %constants
    %visco = 1e-6;
    D = 1.44e-7;
    
    %defines times
    time_fast0 = [0:1:length(DATA.P_fast)-1]/DATA.fs_fast;
    time_slow0= [0:1:length(DATA.P_slow)-1]/DATA.fs_slow;
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
    C1f = DATA.C1_fast(iipf);
    sh = DATA.sh1(iipf);
    Ax = DATA.Ax(iipf);
    Ay = DATA.Ay(iipf);
    Wf = DATA.W_fast(iipf);
    AA = [Ax,Ay];
    Chl = DATA.Chlorophyll(iipf);
    Turb = DATA.Turbidity(iipf);
  
    %pshx = Ax./Wf;
    %pshy = Ay./Wf;

    %gets slow response sensors
    times = time_slow0(iips);
    Ps = DATA.P_slow(iips);
    T_JAC = DATA.JAC_T(iips);
    C_JAC = DATA.JAC_C(iips);
    Incl_X = DATA.Incl_X(iips);
    Incl_Y = DATA.Incl_Y(iips);

    

    %sets maximum depth
    if ~ isfield(info,'pmax')
        info.pmax = round(max(Ps)-info.dpD/2)+1;
    end

    %calculates salinity and density
    [rhos,Ss,depths]=rho_salinity_Geneva(T_JAC,C_JAC,Ps);
    
    %slow data for output
    SLOW.date = date;
    SLOW.filename = dataf_name;
    SLOW.pres = Ps;
    SLOW.depth = depths;
    SLOW.T = T_JAC;
    SLOW.C = C_JAC;
    SLOW.S = Ss;
    SLOW.rho = rhos;
    SLOW.Incl_X = Incl_X;
    SLOW.Incl_Y = Incl_Y;
    
    %calculates displacements for thorpe length
    if ismember(info.prof_dir,'down')
        [sort_rho, isd] = sort(rhos);
        displ = Ps - Ps(isd);
    else
        [Ps1,is1] = sort(Ps);
        rhos1 = rhos(is1);
        [sort_rho, isd] = sort(rhos1);
        displ = Ps1 - Ps1(isd);
    end
    

    %binned temperature, salinity and density
    %defines the presure vector where to calculate
    pres = [info.pmin:info.dp:info.pmax];
    BIN.date = date;
    BIN.filename = dataf_name;
    BIN.pres = pres;
    BIN.T=pres_av(Ps,T_JAC,pres,info.dp,2.7);
    BIN.C = pres_av(Ps,C_JAC,pres,info.dp,2.7);
    BIN.S = pres_av(Ps,Ss,pres,info.dp,2.7);
    BIN.rho = pres_av(Ps,rhos,pres,info.dp,2.7);
    BIN.grT = mean_grad(Ps,T_JAC,pres,info.dpD);
    BIN.N2 = -9.81*mean_grad(Ps,rhos,pres,info.dpD)/1000;
    BIN.Chl = pres_av(Pf,Chl,pres,info.dp,2.7);
    BIN.Turb = pres_av(Pf,Turb,pres,info.dp,2.7);
    BIN.LT = sqrt(pres_av(Ps,displ.^2,pres,info.dpD,0));


    %matches the FP07 to the high accuracy sensors
    T1f = calibration_FP07(timef,T1f, times, T_JAC,fsf,fss);
    T2f = calibration_FP07(timef,T2f, times, T_JAC, fsf,fss);
 
    %plots raw
    figure(1)
    clf
    subplot(4,1,1)
    plot(DATA.P_slow)
    ylabel('p (db)')
    hold on
    plot(iips,Ps)
    set(gca,'xticklabel',[])
    subplot(4,1,2)
    plot( sh)
    ylabel('sh (1/s)')
    set(gca,'xticklabel',[])
    subplot(4,1,3)
    plot( T1f)
    hold on
    plot( T2f)
    ylabel('T fast (°C)')
    set(gca,'xticklabel',[])
    subplot(4,1,4)
    plot(Wf)
    ylabel('w (m/s)')
    
    saveas(gcf,['profile_',dataf_name,'_p',num2str(inp,'%02d'),'.png'])    
    

    %filters shear and highpasses microstructure
    [bh,ah] = butter(1, info.HP_cut/(fsf/2), 'high');
    [sh, ~, ~, ~ ] =  despike(sh, despike_sh(1), despike_sh(2), fsf, round(despike_sh(3)*fsf));
    sh_hp = filter(bh, ah, sh);
    sh_hp = flipud(sh_hp);
    sh_hp = filter(bh, ah, sh_hp);
    sh_hp = flipud(sh_hp);

    %makes strange things at the borders (avoid?)
    T1f_hp = filter(bh, ah, T1f);
    T1f_hp = flipud(T1f_hp);
    T1f_hp = filter(bh, ah, T1f_hp);
    T1f_hp = flipud(T1f_hp);

    T2f_hp = filter(bh, ah, T2f);
    T2f_hp = flipud(T2f_hp);
    T2f_hp = filter(bh, ah, T2f_hp);
    T2f_hp = flipud(T2f_hp);
    
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
    BIN.epsSH = nan(1,length(pres));
    BIN.epsSH2 = nan(1,length(pres));
    BIN.KBSH = nan(1,length(pres));
    BIN.W = nan(1,length(pres));
    
    BIN.Xic1 = nan(1,length(pres));
    BIN.Xif1 = nan(1,length(pres));
    BIN.KB1 = nan(1,length(pres));
    BIN.sXif1 = nan(1,length(pres));
    BIN.sKB1 = nan(1,length(pres));
    BIN.Xiv1 = nan(1,length(pres));
    BIN.maxK1 = nan(1,length(pres));
    BIN.epsT1 = nan(1,length(pres));
    BIN.MAD1 = nan(1,length(pres));
    BIN.MADf1 = nan(1,length(pres));
    BIN.LKH1 = nan(1,length(pres));
    BIN.LKHratio1 = nan(1,length(pres));
    BIN.MAD1 = nan(1,length(pres));

    BIN.Xic2 = nan(1,length(pres));
    BIN.Xif2 = nan(1,length(pres));
    BIN.KB2 = nan(1,length(pres));
    BIN.sXif2 = nan(1,length(pres));
    BIN.sKB2 = nan(1,length(pres));
    BIN.Xiv2 = nan(1,length(pres));
    BIN.maxK2 = nan(1,length(pres));
    BIN.epsT2 = nan(1,length(pres));
    BIN.MAD2 = nan(1,length(pres));
    BIN.MADf2 = nan(1,length(pres));
    BIN.LKH2 = nan(1,length(pres));
    BIN.LKHratio2 = nan(1,length(pres));
    BIN.MAD2 = nan(1,length(pres));

    
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
        
        %ODAS epsilon calculation
        try
            [BIN.epsSH(i), BIN.W(i)] = dis_spec_ODAS(sh_hp(jp), AA(jp,:), fsf, fss, Wf(jp), T1f_hp(jp),timef(jp), Pf(jp), visco,Nfft, overlap, PLOT);
            %BIN.W(i) = nanmean(Wf(jp));
            BIN.KBSH(i)=1/(2*pi())*(BIN.epsSH(i)/(visco*D^2))^(1/4); 
        end
        
        %with my function
        %[BIN.epsSH2(i),~]=dis_spec(Pf(jp),sh_hp(jp),Ay(jp),1,14,30,visco, Nfft, overlap, PLOT);
        
        %FP07 calculations
        try
            [BIN.Xic1(i),BIN.Xif1(i),BIN.KB1(i), BIN.sXif1(i), BIN.sKB1(i),BIN.Xiv1(i), BIN.maxK1(i), BIN.MAD1(i), BIN.MADf1(i),BIN.LKH1(i), BIN.LKHratio1(i)] =Xi_spec(Pf(jp),T1f_hp(jp),info.minKT,info.fAA,BIN.KBSH(i),mean(abs(Wf(jp))),info.noisep_T1,  Nfft, overlap,info.Tspec,info.Tmethod,PLOT);
            BIN.epsT1(i) = visco*D^2*(2*pi()*BIN.KB1(i))^4;
        end
        try
            [BIN.Xic2(i),BIN.Xif2(i),BIN.KB2(i), BIN.sXif2(i), BIN.sKB2(i),BIN.Xiv2(i),BIN.maxK2(i), BIN.MAD2(i), BIN.MADf2(i),BIN.LKH2(i), BIN.LKHratio2(i)] =Xi_spec(Pf(jp),T2f_hp(jp),info.minKT,info.fAA,BIN.KBSH(i),mean(abs(Wf(jp))),info.noisep_T2,  Nfft, overlap,info.Tspec,info.Tmethod,PLOT);
            BIN.epsT2(i) = visco*D^2*(2*pi()*BIN.KB2(i))^4;
        end
    end


    BIN.Krho = 0.2*BIN.epsSH.*(BIN.N2).^-1;
    BIN.KTf1 = 0.5*BIN.Xif1.*(BIN.grT).^-2;
    BIN.KTf2 = 0.5*BIN.Xif2.*(BIN.grT).^-2;
    BIN.KTc1 = 0.5*BIN.Xic1.*(BIN.grT).^-2;
    BIN.KTc2 = 0.5*BIN.Xic2.*(BIN.grT).^-2;
    BIN.LO = (BIN.epsSH./BIN.N2.^(3/2)).^(0.5);

    pmin = info.pmin - info.dpD/2;
    pmax = info.pmax + info.dpD/2;
    
    %plots profile
    figure(2)
    clf
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperSize', [29 20]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 29 20]);

    ax1=subplot(2,4,1);
    plot(T_JAC, Ps,'linewidth',2)
    xlabel('T (°C)')
    ylabel('p (db)')
    set(gca,'YDir','reverse')
    ylim([pmin,pmax])
    ax2=axes('Position',get(ax1,'Position'));
    set(ax2,'box','off')
    plot(rhos, Ps,'r', 'parent' , ax2,'linewidth',2)
    set(ax2,'XAxisLocation','top',...
            'YAxisLocation','right',...
            'Color','none',...
            'XColor','r','YColor','k');
    yticklabels([])
    ylim([pmin,pmax])
    xlabel('dens (kg m^{-3})')
    set(gca,'yticklabel',[])
    set(gca,'YDir','reverse')
        
    subplot(2,4,2)
    plot(BIN.W,BIN.pres,'o-','linewidth',2,'markersize',4)
    xlabel('W (db/s)')
    ylim([pmin,pmax])
    set(gca,'yticklabel',[])
    set(gca,'YDir','reverse')
    grid('on')
    
    
    subplot(2,4,3)
    plot(BIN.epsT1,BIN.pres,'o-','linewidth',2,'markersize',4)
    hold on
    plot(BIN.epsT2,BIN.pres,'o-','linewidth',2,'markersize',4)
    plot(BIN.epsSH,BIN.pres,'o-','linewidth',2,'markersize',4)
    legend('T01','T02','shear','fontsize',7,'location','southeast')
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
    l1=plot(BIN.Xif1,BIN.pres,'o-','linewidth',2,'markersize',4);
    hold on
    l2=plot(BIN.Xif2,BIN.pres,'o-','linewidth',2,'markersize',4);
    plot(BIN.Xic1,BIN.pres,'--','linewidth',1,'color',l1.Color())
    plot(BIN.Xic2,BIN.pres,'--','linewidth',1,'color',l2.Color())
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
    plot(BIN.KTf1,BIN.pres,'o-','linewidth',2,'markersize',4)
    hold on
    plot(BIN.KTc2,BIN.pres,'o-','linewidth',2,'markersize',4)
    plot(BIN.Krho,BIN.pres,'o-','linewidth',2,'markersize',4)
    legend('T01 (O&C)','T02 (O&C)','she (O)','fontsize',7,'location','southeast')
    xlabel('K (m^2 s^{-1})')
    ylabel('p (db)')
    set(gca,'YDir','reverse')
    set(gca,'xscale','log')
    xlim([1e-7,1e0])
    ylim([pmin,pmax])
     set(gca,'xtick',10.^[-7:2:0])
    grid('on')
    
    subplot(2,4,6)
    l1=plot(BIN.KB1,BIN.pres,'o-','linewidth',2,'markersize',4);
    hold on
    l2=plot(BIN.KB2,BIN.pres,'o-','linewidth',2,'markersize',4);
    plot(BIN.KBSH,BIN.pres,'o-','linewidth',2,'markersize',4)
    plot(BIN.maxK1,BIN.pres,'--','linewidth',1,'color',l1.Color())
    plot(BIN.maxK2,BIN.pres,'--','linewidth',1,'color',l2.Color())
    xlabel('K_B (m)')
    set(gca,'YDir','reverse')
    %set(gca,'xscale','log')
    grid('on')
    ylim([pmin,pmax])
    
% $$$     plot(BIN.LT,BIN.pres,'o-','linewidth',2,'markersize',4)
% $$$     hold on
% $$$     plot(BIN.LO,BIN.pres,'o-','linewidth',2,'markersize',4)
% $$$     legend('L_T','L_O','fontsize',7,'location','southeast')
% $$$     xlabel('L_T, L_O (m)')
% $$$     set(gca,'YDir','reverse')
% $$$     set(gca,'xscale','log')
% $$$     grid('on')
% $$$     ylim([pmin,pmax])
    
    subplot(2,4,7)
    plot(BIN.MADf1, BIN.pres,'o-','linewidth',2,'markersize',4)
    hold on
    plot(BIN.MADf2, BIN.pres,'o-','linewidth',2,'markersize',4)
    line([1.15,1.15],ylim, 'color','k')
    xlabel('MAD')
    xlim([0,2])
    set(gca,'yticklabel',[])
    set(gca,'YDir','reverse')
    grid('on')
    ylim([pmin,pmax])
    
    subplot(2,4,8)
    plot(BIN.LKHratio1, BIN.pres,'o-','linewidth',2,'markersize',4)
    hold on
    plot(BIN.LKHratio2, BIN.pres,'o-','linewidth',2,'markersize',4)
    line([2,2],ylim, 'color','k')
    xlabel('Likelihood ratio')
    set(gca,'yticklabel',[])
    set(gca,'YDir','reverse')
    grid('on')
    ylim([pmin,pmax])
    
    saveas(gcf,[dataf_name,'_P',num2str(inp,'%02d'),'.png'])
    

    %hold on
    %semilogx(epsilon2,-pres)
    %semilogx(epsilonN,-pres)