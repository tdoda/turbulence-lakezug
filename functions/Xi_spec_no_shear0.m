function [Xiv,Xif,KBT,sXif,sKBT,maxK, MADf,MLKH, LKHratio]=Xi_spec_no_shear(pres,x0,K1,fn,W,noisep, sL, sOV,Tdis,method,plt)
    %INPUT:
    %pres: pressure vector
    %x0: temperature vector
    %K1: minimium wavenumber for  integration
    %fn: maximum frequency for calculations
    %(KB: KB wavenumber determined from shear probe)
    %W: mean velocity
    %noise p: two (four for microCTD) parameters noise function
    %sL: length of segments for fft 
    %sOV: overlap for fft
    %Tdis: type of spectrum (B: batchelor, K: Kraichnan)
    %method: (B)ieito or (O)scar
    %plt: if 1 shows the spectra
    
    %OUPUT
    %Xiv: Chi initial guess from spectral integration
    %Xif, sXif: Chi after fit and error
    %KBT, sKBT: KB from fit and error
    %maxK: maximum K used for calculation
    %MADf: MAD from fit
    %MLKH: fit (maximum) likelilood
    %LKHratio: likelihood ratio
    Xi=NaN;Xif=NaN;KBT=NaN;sXif=NaN;sKBT=NaN;MADf = NaN; MLKH=nan; LKHratio = nan;
    Pr = 7.56;
    if nargin<11
        plt=0;
    end
    if nargin<10
        %sets the method
        %B: bieito
        %O: Oscar
        method = 'B';
    end
    if nargin<9
        %sets spectrum type (B: Batchelor, K: Kraichnan)
        %only for my method
        Tdis = 'B';
    end
    if method ~= 'B' && method ~= 'O'
        method = 'B';
    end
    
    KBT=NaN;
    
    %detrend
    x = detrend(x0,'linear');
    
    I=find(isfinite(x));
    x=x(I);
    pres=pres(I);
    
    Fs=length(pres)./ (max(pres)-min(pres));
    
    if isempty(x) | sum(x==0)==length(pres)
        return;
    end
    
    dof = round( 2*((length(x)-sL)/sOV +1));
    [PSDT,K] = csd_odas(x-nanmean(x),x-nanmean(x),sL,Fs,[],sOV,'linear');
     PSD = (2*pi()*K).^2.*PSDT; %gradient spectra
    fr = W*K;
    
    %correccion varianza Goto2016/tau 10 ms
%     PSDnc = PSD;
     %tau = 0.005*W.^-0.32;
     tau = 0.010;
     H = 1./(1+(2*pi()*tau*fr).^2).^2;
     PSD = PSD./H;
     
    %correction varianza microCTD
    %F0 = 25*sqrt(W);
    %tau = (2*pi()*F0/sqrt(sqrt(2)-1))^(-1);
    %H = 1./(1 + (2*pi()*tau*fr).^2).^2;
    %PSD = PSD./H;
    
    %fits noise function
    %Sn = FP07noise0(noisep,fr);
    Sn = FP07noise(noisep,fr); %microCTD
    Sn = Sn*W;
    Sn = Sn.*(2*pi()*K).^2;
    Sn = Sn./H;
    Sn(1) = 0;

    
    mpres=mean(pres);
    Kn = fn/W;    
    iK1=find(K>K1,1,'first');
    iKn=find(K<=Kn,1,'last');


    %deletes undesired part of the spectrum
    K0 = K;
    PSD0 = PSD;
    Sn0 = Sn;
    
    K = K(1:iKn);
    maxK = K(end);
    PSD = PSD(1:iKn);
    H = H(1:iKn);
    Sn = Sn(1:iKn);
    
    %variance in the noise free part determined from the noise model
    iKnM0 = find(PSD<2*Sn | H<0.1,1,'first');
    ikcor = iK1:iKnM0;
    Xiv =6*1.4e-7*sum( (PSD(ikcor(2:end))+PSD(ikcor(1:end-1))).* (K(ikcor(2:end))-K(ikcor(1:end-1))) )/2;
    Xiv = Xiv - 6*1.4e-7*sum( (Sn(ikcor(2:end))+Sn(ikcor(1:end-1))).* (K(ikcor(2:end))-K(ikcor(1:end-1))) )/2;
    %continues only if Xiv is detectable
    if Xiv>0

         %fits parameters in the noise free region
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
       if method == 'B'  
           KF = K(iKn);
           Ktest = linspace(max([5*K1,KF/20]),KF*5,40);
           dK = Ktest(2)-Ktest(1);
           for i = 1:length(Ktest)
                ks = 0.04*Pr^(-0.5)*Ktest(i);
                Kf1 = max([ks,K1]);
                iKf1 = find(K>=Kf1,1,'first');
                ikfit = iKf1:iKn;
                [cost(i), Xif0(i)] =  cost_T_fit(K, PSD, Ktest(i), ikfit, Tdis, dof, Sn);
           end
           LKHtest = - cost;
           ML = max(LKHtest);
           iML = find(LKHtest == ML);
           KB00 = Ktest(iML);
           MLp = -cost_T_fit(K, PSD, KB00+dK, ikfit, Tdis, dof, Sn);
           MLm = -cost_T_fit(K, PSD, KB00-dK, ikfit, Tdis, dof, Sn);
           deltak0=abs((2*dK)/sqrt(2*ML-MLm-MLp));
           deltak = max([deltak0,dK]);

           Ktest = linspace(KB00-deltak,KB00+deltak,40);
           clear cost Xif0
           for i = 1:length(Ktest)
               ks = 0.04*Pr^(-0.5)*Ktest(i);
               Kf1 = max([ks,K1]);
               iKf1 = find(K>=Kf1,1,'first');
               ikfit = iKf1:iKn;
               [cost(i), Xif0(i)] =  cost_T_fit(K, PSD, Ktest(i), ikfit, Tdis, dof, Sn);
           end
           LKHtest = - cost;
           MLKH = max(LKHtest);
           iML = find(LKHtest == MLKH);
           Xif = Xif0(iML);
           KBT = Ktest(iML);        
        
           BATf = Tspec(Tdis,Xif,KBT,K);
           ks = 0.04*Pr^(-0.5)*KBT;
           Kf1 = max([ks,K1]);
           iKf1 = find(K>=Kf1,1,'first');
           ikfit = iKf1:iKn;
           MADf = meanabsdev( PSD(ikfit), BATf(ikfit), Sn(ikfit) );
           
            %calculates uncertainties in the fitting parameters
           [MLp,Xip] = cost_T_fit(K, PSD, KBT+deltak, ikfit, Tdis, dof, Sn);
           [MLm,Xim] = cost_T_fit(K, PSD, KBT-deltak, ikfit, Tdis, dof, Sn);
           MLp = - MLp;
           MLm = - MLm;
           sKBT=abs((2*deltak)/sqrt(2*MLKH-MLm-MLp));
           sXif = abs((Xip-Xim)/sqrt(2*MLKH-MLm-MLp));
        
           %fits to polynom (avoiding the noisy part)
           ikfitA = ikfit(1):min([ikfit(end),find(BATf<2*Sn,1,'first')]);
           logK = log(K(ikfitA));
           logS = log(PSD(ikfitA));
           pp=polyfit(logK, logS,1);
           Sm = exp(polyval(pp, log(K)));
           LKHpol = -cost_MLE(PSD(ikfit), Sm(ikfit), dof, Sn(ikfit));
        
           LKHratio = MLKH - LKHpol;
       else
            %EPFL method          
           [LKHratio,MLKH,KBT,~,QUAL]=Fit_kB_OSS_Carpenter(K,PSD,Sn,1);
           Xif = QUAL.X_T1;
           BATf = Tspec(Tdis,Xif,KBT,K);
           MADf = QUAL.MAD1;
           ikfit = iK1:iKn;
       end
           
        
        %plots
        iKBT=find(K0>=KBT,1,'first');if isempty(iKBT), iKBT=length(K); end
        
        if plt~=0
            figure(1)
            clf
            set(gcf, 'PaperUnits', 'centimeters');
            set(gcf, 'PaperSize', [18 8]); 
            set(gcf, 'PaperPositionMode', 'manual') ;
            set(gcf, 'PaperPosition', [0 0 17 8]);
            axes('position',[0.1,0.12,0.2,0.8])
            plot(x0-mean(x0),pres,'linestyle','--','color',[0.5,0.5,0.5])
            hold on
            plot(x,pres,'k')
            grid('on')
            xlabel('T" (^oC)')
            ylabel('p (db)')
 
            axes('position',[0.4,0.12,0.5,0.8])
            loglog(K0,PSD0, 'color' ,[0.5,0.5,0.5])
            hold on
            loglog(K,PSD,'k')
            loglog(K(ikfit),PSD(ikfit),'ko', 'markersize',4)     
            ylim([10^-9,1e3])
            loglog(K,Sn+BATf, 'color',[0.5,0.5,0.5], 'linewidth',2);
            loglog(K0,Tspec(Tdis, Xif, KBT,K0), 'color',[0.5,0.5,0.5], 'linewidth',1,'linestyle','--');
            loglog(K0,Sn0, 'color','r', 'linewidth',1,'linestyle','-');
            line([K0(iKn) K0(iKn)], ylim,'color','k','linestyle','--')
            line([K0(iK1) K0(iK1)], ylim,'color','k','linestyle','--')
            line([KBT KBT], ylim,'color','k')
            text(1.1*KBT,0.1,'K_B')
            if Tdis=='B'
                text(1.2,10^2,'Batchelor Spectrum','horizontalalignment','left')
            elseif Tdis=='K'
                 text(1.2,10^2,'Kraichnan Spectrum','horizontalalignment','left')
            end           

            text(1.2,10^1,['\chi_{var.} = ', num2str(Xiv,'%1.3e'),' K^2/s'],'horizontalalignment','left')
            text(1.2,10^0.25,['\chi_{fit} = ', num2str(Xif,'%1.3e'),'\pm',num2str(sXif,'%1.3e'),' K^2/s'],'horizontalalignment','left')
            text(1.2,10^-2,['K_{B}= ', num2str(KBT,'%1.0f'),'\pm',num2str(sKBT,'%1.0f'),' cpm'],'horizontalalignment','left')
            text(120, 1e-7, ['MADf = ', num2str(MADf, '%1.2f')],'horizontalalignment','left')
            text(120, 1e-6, ['LKHratio = ', num2str(LKHratio, '%1.1f')],'horizontalalignment','left')
            title([num2str(mpres), 'm'],'Fontsize',12)
            xlabel('K (cpm)')
            ylabel('PSD (K^2 s^{-2} cpm^{-1})')
            xlim([1,1000])

            %text(100,10^-5,['K_B^{Batchelor} = ', num2str(KBT,'%1.0f'),' cpm'],'horizontalalignment','left')
            %saveas(gcf,'Tspec.png')
            pause()
            

        end

end
