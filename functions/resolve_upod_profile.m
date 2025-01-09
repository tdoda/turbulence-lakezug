function BIN = resolve_upod_profile(DATA, inp, info,dataf_name,PLOT)
    DATA
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
        info.prof_dir = 'up';
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
    if ~ isfield(info,'noisep_T')
        %info.noisep_T1 = [-10.24,-0.89,info.fAA];
        info.noisep_T = [-10.5,-0.6,info.fAA];
    else
        info.noisep_T = [info.noisep_T1, info.fAA];
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

    %adds some more noise
    %info.noisep_T1 = info.noisep_T1 + [0.25,0,0];
    %info.noisep_T2 = info.noisep_T2 + [0.25,0,0];
    
    %constants
    %visco = 1e-6;
    D = 1.44e-7;
    
    %defines times
    time_fast0 = [0:1:length(DATA.P_fast)-1]/DATA.fs_fast;
        
    %gets profiles
    %gets profiles
    iPf0 = get_profile(DATA.P_fast,DATA.W_fast,0,info.minvel_detect,info.prof_dir,info.mindur_detect,DATA.fs_fast);
    NP = size(iPf0);
    NP = NP(2);
    if inp<0 > inp>NP
        fprintf('Wrong number of profiles')
        return
    end
    %gets index for the desired profile
    if NP>0
        iipf = iPf0(1,inp):iPf0(2,inp);
    else
        iipf = [1:1:length(DATA.P_fast)];
    end

    %sampling frequencies
    fsf = DATA.fs_fast;  
    
    %date
    date = datenum(DATA.Year, DATA.Month, DATA.Day, DATA.Hour, ...
                   DATA.Minute, DATA.Second );
    date = date + iPf0(1,inp)/DATA.fs_fast/60/60/24;
    datestr(date)

    %gets fast response sensors
    timef = time_fast0(iipf);
    Pf = DATA.P_fast(iipf);
    T1f = DATA.T2_fast(iipf);
    T2f = DATA.T3_fast(iipf);
    Wf = DATA.W_fast(iipf);
    

    
    %from velocity
    %acc = nan(size(Wf));
    %acc(2:end-1) = (Wf(3:end)-Wf(1:end-2))*(fsf/2);
    %acc(1) = (Wf(2)-Wf(1))*fsf;
    %acc(end) = (Wf(end)-Wf(end-1))*fsf;
    %pshear = acc./Wf;
    
    
    %pressure correction
    if isfield(info,'pres_cor')
       Pf = Pf*info.pres_cor(1) + info.pres_cor(2);
    end


    %sets maximum depth
    if ~ isfield(info,'pmax')
        info.pmax = round(max(Pf)-info.dpD/2)+1;
    end

    
    %calculates displacements for thorpe length
    if ismember(info.prof_dir,'down')
        
        [sort_uT1, ist1] = sort(T1f,'descend');
        displuT1 = Pf - Pf(ist1);
        
        [sort_uT2, ist2] = sort(T2f,'descend');
        displuT2 = Pf - Pf(ist2);
    else
        
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
    BIN.T = pres_av(Pf,T1f,pres,info.dp,2.7);
    BIN.grT = mean_grad(Pf,T1f,pres,info.dpD);
    %BIN.grT = mean_grad(Pf,sort(T1f,'descend'),pres,info.dpD);
    BIN.N2 = abs(9.81*sw_alpha(0.2*ones(size(BIN.T)),BIN.T,BIN.pres).*BIN.grT);
    BIN.LTuT1 = sqrt(pres_av(Pf,displuT1.^2,pres,info.dpD,0));
    BIN.LTuT2 = sqrt(pres_av(Pf,displuT2.^2,pres,info.dpD,0));
    
    
    %plots raw
    figure(1)
    clf
    subplot(3,1,1)
    plot(DATA.P_fast)
    ylabel('p (db)')
    hold on
    plot(iipf,Pf)
    set(gca,'xticklabel',[])
    set(gca,'xticklabel',[])
    subplot(3,1,2)
    plot( T1f)
    hold on
    plot( T2f)
    ylabel('T fast (°C)')
    set(gca,'xticklabel',[])
    subplot(3,1,3)
    plot(Wf)
    ylabel('w (m/s)')
    saveas(gcf,['profile_',dataf_name,'_p',num2str(inp,'%02d'),'.png'])    
    

    %filters shear and highpasses microstructure
    %makes strange things at the borders (avoid?)
    mW = (max(Pf)-min(Pf))/(max(timef)-min(timef));
    if info.k_HP_cut>0
        f_HP_cut = info.k_HP_cut*mW;
        [bh,ah] = butter(1, f_HP_cut/(fsf/2), 'high');
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
    
    %defines output variables
    
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
        %FP07 calculations
        try
            [BIN.Xiv1(i),BIN.Xif1(i),BIN.KB1(i),BIN.fit_flag_T1(i),BIN.sXif1(i), BIN.sKB1(i), BIN.MADf1(i),BIN.MADc1(i),BIN.LKH1(i), BIN.LKHratio1(i), BIN.maxK1(i)] =Xi_spec_no_shear(Pf(jp),T1f_hp(jp),info.minKT,info.fAA,mean(abs(Wf(jp))),info.noisep_T,  Nfft, overlap,info.Tspec,info.Tmethod,info.time_res,info.peak_rem_T1,PLOT);
            BIN.epsT1(i) = visco*D^2*(2*pi()*BIN.KB1(i))^4;
        end
        
        try
            [BIN.Xiv2(i),BIN.Xif2(i),BIN.KB2(i),BIN.fit_flag_T2(i), BIN.sXif2(i), BIN.sKB2(i), BIN.MADf2(i),BIN.MADc2(i),BIN.LKH1(i), BIN.LKHratio2(i), BIN.maxK2(i)] =Xi_spec_no_shear(Pf(jp),T2f_hp(jp),info.minKT,info.fAA, mean(abs(Wf(jp))),info.noisep_T,  Nfft, overlap,info.Tspec,info.Tmethod,info.time_res, info.peak_rem_T2,PLOT);
            BIN.epsT2(i) = visco*D^2*(2*pi()*BIN.KB2(i))^4;
        end
        
        %BIN.pseps(i) = 0.5*7.5*visco*(var(psh_x(jp)) + var(psh_y(jp)));  
        %BIN.pseps(i) = 0.5*7.5*visco*var(pshear(jp)); 
    end


    BIN.KOsb1 = 0.2*BIN.epsT1.*(BIN.N2).^-1;
    BIN.KOsb2 = 0.2*BIN.epsT2.*(BIN.N2).^-1;
    BIN.KTf1 = 0.5*BIN.Xif1.*(BIN.grT).^-2;
    BIN.KTf2 = 0.5*BIN.Xif2.*(BIN.grT).^-2;
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
    plot(T1f - nanmean(T1f), Pf,'linewidth',1)
    hold on
    plot(T2f - nanmean(T2f), Pf,'linewidth',1)
    xlabel('relative T (°C)')
    ylabel('p (db)')
    set(gca,'YDir','reverse')
    ylim([pmin,pmax])
        
    ax3=subplot(2,4,2)
    plot(BIN.W,BIN.pres,'-','linewidth',2,'markersize',4)
    hold on
    plot(abs(Wf), Pf, 'color','r')
    xlabel('W (db/s)')
    ylim([pmin,pmax])
    set(gca,'yticklabel',[])
    set(gca,'YDir','reverse')
    grid('on')
    title(datestr(date))
    
    
    subplot(2,4,3)
    plot(BIN.epsT1,BIN.pres,'-','linewidth',1,'markersize',4)
    hold on
    plot(BIN.epsT2,BIN.pres,'-','linewidth',1,'markersize',4)
    %plot(BIN.pseps,BIN.pres,'color',[0.5,0.5,0.5],'linewidth',2)
    legend('T01','T02','pseps','fontsize',7,'location','southeast')
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
    legend('T01_f','T02_f','T01_v','T02_v','fontsize',7,'location','southeast')
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
    legend('T01 (\chi)','T02 (\chi)','T01 (\epsilon)','T02 (\epsilon)','fontsize',7,'location','southeast')
    xlabel('K (m^2 s^{-1})')
    ylabel('p (db)')
    set(gca,'YDir','reverse')
    set(gca,'xscale','log')
    xlim([1e-9,1e0])
    ylim([pmin,pmax])
     set(gca,'xtick',10.^[-7:2:0])
    grid('on')
    
    subplot(2,4,6)
    
     plot(BIN.LTuT1,BIN.pres,'-','linewidth',1,'markersize',4)
     hold on
     plot(BIN.LTuT2,BIN.pres,'-','linewidth',1,'markersize',4)
     plot(BIN.LO1,BIN.pres,'-','linewidth',1,'markersize',4)
     hold on
     plot(BIN.LO2,BIN.pres,'-','linewidth',1,'markersize',4)
     legend('L_T^{1}','L_T^{2}','L_O^{1}','L_O^{2}','fontsize',7,'location','southeast')
     xlabel('L_T, L_O(m)')
     set(gca,'YDir','reverse')
     set(gca,'xscale','log')
     set(gca,'yticklabel',[])
     grid('on')
     ylim([pmin,pmax])
    
    subplot(2,4,7)
    plot(BIN.MADf1, BIN.pres,'-','linewidth',1,'markersize',4)
    hold on
    plot(BIN.MADf2, BIN.pres,'-','linewidth',1,'markersize',4)
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
    
    saveas(gcf,[dataf_name,'_P_',num2str(inp,'%02d'),'.png'])
    

    %hold on
    %semilogx(epsilon2,-pres)
    %semilogx(epsilonN,-pres)