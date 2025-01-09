function [Xiv,Xif,KBT,fit_flag, sXif,sKBT,MADf,MADc,MLKH,LKHratio,maxK]=Xi_spec_no_shear(pres,x0,K1,fn,W,noisep, sL, sOV,Tdis,method,time_res,peak_rem,plt)
    %INPUT:
    %pres: pressure vector
    %x0: temperature vector
    %K1: minimium wavenumber for  integration
    %fn: maximum frequency for calculations
    %(KB: KB wavenumber determined from shear probe, if 0 does not calculate)
    %W: mean velocity
    %noise p: two (four for microCTD) parameters noise function
    %sL: length of segments for fft 
    %sOV: overlap for fft
    %Tdis: type of spectrum (B: batchelor, K: Kraichnan)
    %method: (B)ieito or (O)scar
    %time_res: time response in ms
    %plt: if !=0 shows the spectra
    
    %OUPUT
    %Xiv: Xi from spectral integration in the well resolved part  
    %XiC: Chi corrected with epsilon from shear
    %Xif: Chi after fit
    %KBT: KB from fit
    %fit_flag: flag for Batchelor fit, 1 if good
    %sXif, sKBT: fit uncertainties
    %MAD: MAD of the theoretical spectrum using KB from shear    
    %MADf: MAD from fit
    %MLKH: fit (maximum) likelilood
    %LKHratio: likelihood ratio
    %MADc: theoretical MAD    
    %maxK: maximum K considered good
    Xiv = NaN;XiC=NaN;Xif=NaN;KBT=NaN;sXif = nan; sKBT = nan; MAD = NaN; MADf = NaN;MLKH=nan; LKHratio = nan; maxK = nan; fit_flag = 0;
    Pr = 7.56;
    if nargin<12
        plt=0;
    end
    if nargin<11
        time_res = nan;
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
    
    %detrend
    x = detrend(x0,'linear');
    
    I=find(isfinite(x));
    x=x(I);
    pres=pres(I);
    
    Fs=length(pres)./ (max(pres)-min(pres));
    
    if isempty(x) | sum(x==0)==length(pres)
        return;
    end
    
    %degrees of freedom acording to me
    %Nseg = floor((length(x)-sOV)/(sL-sOV));
    %dof = 2*Nseg;
    %more complex (methods in oceanography book)
    %NsegM = floor(length(x)/(sL/2));
    %dof = 0.92*8/3*NsegM;
    %according to ODAS
    NsegM = floor(length(x)/(sL/2));
    %dof = 2*2*(9/11)*NsegM; %according to ODAS
    dof = 1.9*NsegM; %ODAS 4.4
    MADc = sqrt(2/dof);
    [PSDT,K] = csd_odas(x-nanmean(x),x-nanmean(x),sL,Fs,[],sOV,'linear');
    fr = W*K;
    ii = find(K>peak_rem(1) & K<peak_rem(2));
    fr(ii)=[];
    K(ii)=[];
    PSDT(ii) = [];
    PSD = (2*pi()*K).^2.*PSDT; %gradient spectra
    

    %correccion varianza Goto2016/ tau = 10ms
     %tau = 0.005*W.^-0.32;
%     tau = 0.010;
%     H = 1./(1+(2*pi()*tau*fr).^2).^2;
%     PSD = PSD./H;
%     
    %correction varianza microCTD
    F0 = 25*sqrt(W);
    tau = (2*pi()*F0/sqrt(sqrt(2)-1))^(-1);
    if ~isnan(time_res)
        tau = time_res;
    end
    H = 1./(1 + (2*pi()*tau*fr).^2).^2;
    PSD = PSD./H;
    
    %noise function
    Sn = FP07noise(noisep,fr); %for microCTD with up and down!
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
    iKnM0 = find(PSD<2*Sn | H<0.1 | K==max(K),1,'first');
    ikcor = iK1:iKnM0;
    Xiv = 6*1.4e-7*sum( (PSD(ikcor(2:end))+PSD(ikcor(1:end-1))).* (K(ikcor(2:end))-K(ikcor(1:end-1))) )/2;
    Xin = 6*1.4e-7*sum( (Sn(ikcor(2:end))+Sn(ikcor(1:end-1))).* (K(ikcor(2:end))-K(ikcor(1:end-1))) )/2;
    cont = false;
    if Xiv>1.3*Xin
        cont = true;
    end
    Xiv = Xiv -Xin;

    %continues only if Xiv is detectable
    if cont

         %fits parameters in the noise free region
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
       if method == 'B'  
           KF = K(iKnM0); %iKn
           Ktest = linspace(max([5*K1,KF/5]),KF*5,40); %5*KF
           cost = nan(size(Ktest));
           Xif0 = nan(size(Ktest));
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
           ks = 0.04*Pr^(-0.5)*KB00;
           Kf1 = max([ks,K1]);
           iKf1 = find(K>=Kf1,1,'first');
           ikfit = iKf1:iKn;
           MLp = -cost_T_fit(K, PSD, KB00+dK, ikfit, Tdis, dof, Sn);
           MLm = -cost_T_fit(K, PSD, KB00-dK, ikfit, Tdis, dof, Sn);
           deltak0=abs((2*dK)/sqrt(2*ML-MLm-MLp));
           deltak = max([deltak0,dK]);
           
           kmin2 = max([KB00-deltak,K(iK1+1)]);
           kmax2 = KB00+deltak;%min([KB00+deltak, (K(iKn)-3*dK)/0.04*Pr^0.5 ]);
           Ktest = linspace(kmin2,kmax2,40);
           clear cost Xif0
           cost = nan(size(Ktest));
           Xif0 = nan(size(Ktest));
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
           LKHratio = log10(exp(1))*LKHratio;
         %likelyhood ratio is log10(Pteo/Pexp) and C = log(Pteo), to converto
           
           if LKHratio>2 && MADf<MADc*2 && KBT<3*maxK && abs(sKBT)<0.5*abs(KBT)
             %quality flag for the fit
             fit_flag =1;
           end
       else
            %EPFL method          
           [LKHratio,MLKH,KBT,~,QUAL]=Fit_kB_OSS_Carpenter(K,PSD,Sn,1);
           Xif = QUAL.X_T1;
           BATf = Tspec(Tdis,Xif,KBT,K);
           MADf = QUAL.MAD1;
           ikfit = iK1:iKn;
           Sm = nan(size(K));
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
            %plot(x0-mean(x0),pres,'linestyle','--','color',[0.5,0.5,0.5])
            %hold on
            %plot(x,pres,'k')
            plot(x0-mean(x0),pres,'k')
            axis ij
            grid('on')
            xlabel('T" (^oC)')
            ylabel('p (db)')
 
            axes('position',[0.4,0.12,0.5,0.8])
            loglog(K0,PSD0, 'color' ,[0.5,0.5,0.5])
            hold on
            loglog(K,PSD,'k')
            %loglog(K2,PSD2,'-b')
            if cont
                loglog(K(ikfit),PSD(ikfit),'ko', 'markersize',4)   
                loglog(K(ikfit),Sm(ikfit)+Sn(ikfit),'b')   
                loglog(K,Sn+BATf, 'color',[0.5,0.5,0.5], 'linewidth',2);
                loglog(K0,Tspec(Tdis, Xif, KBT,K0), 'color',[0.5,0.5,0.5], 'linewidth',1,'linestyle','--');
            end
            loglog(K0,Sn0, 'color','r', 'linewidth',1,'linestyle','-');
            ylim([10^-9,1e3])
            line([K0(iKn) K0(iKn)], ylim,'color','k','linestyle','--')
            line([K0(iK1) K0(iK1)], ylim,'color','k','linestyle','--')
            line([KBT KBT], ylim,'color','k')
            text(1.1*KBT,0.1,'K_B')
            if Tdis=='B'
                text(1.2,10^2.5,'Batchelor Spectrum','horizontalalignment','left')
            elseif Tdis=='K'
                 text(1.2,10^2.5,'Kraichnan Spectrum','horizontalalignment','left')
            end           

            text(1.2,10^0.75,['\chi_{fit} = ', num2str(Xif,'%1.3e'),'\pm',num2str(sXif,'%1.3e'),' K^2/s'],'horizontalalignment','left')
            text(1.2,10^0,['\chi_{var} = ', num2str(Xiv,'%1.3e'),' K^2/s'],'horizontalalignment','left')

            text(1.2,10^-1.5,['K_{B}= ', num2str(KBT,'%1.0f'),'\pm',num2str(sKBT,'%1.0f'),' cpm'],'horizontalalignment','left')
            text(100, 1e-7, ['MADf = ', num2str(MADf, '%1.2f'),'(', num2str(2*MADc, '%1.2f'),')'],'horizontalalignment','left')
            text(100, 1e-6, ['LKHratio = ', num2str(LKHratio, '%1.1f')],'horizontalalignment','left')
            title([num2str(mpres), 'm'],'Fontsize',12)
            xlabel('K (cpm)')
            ylabel('PSD (K^2 m^{-2} cpm^{-1})')
            xlim([1,3000])

            fit_flag
            %saveas(gcf,'Tspec.png')
            pause()
            

        end

end
