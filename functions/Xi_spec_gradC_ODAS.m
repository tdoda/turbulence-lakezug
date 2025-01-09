function [Xiv,XiC,Xif,KBT,fit_flag, sXif,sKBT,MAD, MADf,MADc,MLKH,LKHratio,maxK,rangeK]=Xi_spec_gradC(pres,x0,K1,fn,KB,W,noisep, sL, sOV,Tdis,q,method,time_res,time_corr,npoles,peak_rem,plt,time,D,visco,T_dT,T_string,setupstr,int_range,folder_out)
%INPUT:
%pres: pressure vector
%x0: grad temperature vector
%K1: minimium wavenumber for  integration
%fn: maximum frequency for calculations
%(KB: KB wavenumber determined from shear probe, if 0 does not calculate)
%W: mean velocity
%noise p: two (four for microCTD) parameters noise function
%sL: length of segments for fft
%sOV: overlap for fft
%Tdis: type of spectrum (B: batchelor, K: Kraichnan)
%q: turbulent parameter
%method: (B)ieito or (O)scar
%time_res: response time in miliseconds
%time_corr: time correction approach: RSI or as in Kocsis et al., 1999 (See also Goto et al., 2016)
%peak_rem: [fmin,fmax]: to frequencies to remove a noise peak (for T1221,T1228)
%plt: if !=0 shows the spectra
%time: time vector
%D: thermal diffusivity
%visco: viscosity
%T_dT: raw pre-emphasized signal
%T_string: name of the channel 'T1_dT1' or 'T2_dT2'
%setupstr: configuration file
%int_range: Integration range according to Steinbuck 2009 ('S') or Luketina and Imberger 2001 ('L')
%folder_out: output folder for saving figures

%OUTPUT
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
Pr = visco/D;
if nargin<14
    plt=0;
end
if nargin<13
    time_res = nan;
end
if nargin<12
    peak_rem = [0,0];
end
if nargin<11
    %sets the method
    %B: bieito
    %O: Oscar
    method = 'B';
end
if nargin<10
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
dof = 1.9*NsegM; %according to ODAS 4.4
MADc = sqrt(2/dof);

scalar_info.fft_length      = sL;
scalar_info.spec_length     = length(x);
scalar_info.overlap         = length(x)/2;
scalar_info.fs              = 512;
scalar_info.gradient_method = 'high_pass';
scalar_info.f_AA            = 98;
%     sp = get_scalar_spectra_odas(x, [], pres, time', W, scalar_info);
%     fr = sp.F; K=sp.K; PSDT=sp.scalar_spec;

[PSDT,fr] = csd_odas(x-nanmean(x),x-nanmean(x),sL,scalar_info.fs,[],sOV,'linear');
% [PSDT,fr] = pwelch(x-nanmean(x),hann(sL),sOV,sL,scalar_info.fs,'onesided');

PSDT=PSDT*W;
K = fr/W;

ii = find(fr>peak_rem(1) & fr<peak_rem(2));
fr(ii)=[];
K(ii)=[];
PSDT(ii) = [];
PSD=PSDT;

% %correction variance Goto2016/ tau = 10ms
% if strcmp(time_corr,'KOC')
%     tau = time_res*W.^-0.5;
% elseif strcmp(time_corr,'RSI')
%     F0 = 25*sqrt(W);
%     tau = (2*pi()*F0/sqrt(sqrt(2)-1))^(-1);
% elseif strcmp(time_corr,'NAS')
%     tau = time_res*W.^-0.12;
% else
%     error('Error: time_corr can be RSI, KOC or NAS')
% end
% %tau = 0.010;
% %tau=4.1E-3*W^(-0.5);
% 
% % Single or double correction
% if strcmpi(npoles,'single')
%     H = 1./(1 + (2*pi()*tau*fr).^2);
%     H_lim=0.1;
% elseif strcmpi(npoles,'double')
%     H = 1./(1 + (2*pi()*tau*fr).^2).^2;
%     H_lim=0.1;
% end
% 
% PSD = PSD./H;
% %iKmax = find(H>0.1,1,'last');
% 
% %noise function
% noise_info=gradT_noise_odas;
% noise_info.gamma_RSI = 1;
% noise_info.E_n = 3e-9;
% Sn_ODAS = gradT_noise_odas(T_dT, T_string, W, fr, setupstr,noise_info);
% Sn_ODAS = Sn_ODAS*W;
% Sn_ODAS = Sn_ODAS./H;
% Sn_ODAS(1) = 0;
% 
% Snold = FP07noise(noisep,fr); %for microCTD with up and down!
% Snold = Snold*W;
% Snold = Snold.*(2*pi()*K).^2;
% Snold = Snold./H;
% Snold(1) = 0;
% 
% Sn = Sn_ODAS;
% 
mpres=mean(pres);
% Kn = fn/W;
% iK1=find(K>K1,1,'first');
% iKn = find(K<=Kn,1,'last');
% %iKn=min([find(K<=Kn,1,'last'),iKmax]); %if I want to remove
% %strong correction
% iKB=find(K>=KB,1,'first'); if isempty(iKB), iKB=length(K); end
% 
% %deletes undesired part of the spectrum
% K0 = K;
% PSD0 = PSD;
% Sn0 = Sn;
% H0 = H;
% 
% K = K(1:iKn);
% maxK = K(end);
% PSD = PSD(1:iKn);
% H = H(1:iKn);
% Sn = Sn(1:iKn);
% 
% %variance in the noise free part determined from the noise
% %model
% ksfact=0.04; % (0.04 in Steinbuck 2009: https://journals.ametsoc.org/doi/pdf/10.1175/2009JTECHO611.1)
% Snfact=1.55; %Snfact=1.55;   (1.55 as in Goto 2016 DOI: 10.1175/JTECH-D-15-0220.1)
% iKnM0 = find(PSD<Snfact*Sn | H<H_lim | K==max(K),1,'first');
% % iKnM0 = find(PSD<Snfact*Sn,1,'first');
% 
% Ks = ksfact*Pr^(-0.5)*KB;
% K11 = max([Ks,K1]);
% iK11 = find(K>=K11,1,'first'); %if I do like that the fit fails for high KB
% ikcor = iK11:iKnM0;
% Xiv = 6*D*sum( (PSD(ikcor(2:end))+PSD(ikcor(1:end-1))).* (K(ikcor(2:end))-K(ikcor(1:end-1))) )/2;
% Xin = 6*D*sum( (Sn(ikcor(2:end))+Sn(ikcor(1:end-1))).* (K(ikcor(2:end))-K(ikcor(1:end-1))) )/2;
% cont = false;
% if Xiv>1.3*Xin
%     cont = true;
%     Xiv = Xiv -Xin;
% else
%     Xiv = 0;
%     Xif = 0;
% end
% 
% 
% %continues only if Xiv is detectable
% if cont
%     
%     %calculates Xi using theoretical spectrum to correct for unresolved
%     if KB>0
%         BAT=Tspec(Tdis,Xiv,KB,K,D,q);
%         XiT=6*D*sum( (BAT(ikcor(2:end))+BAT(ikcor(1:end-1))).*(K(ikcor(2:end))-K(ikcor(1:end-1))) )/2;
%         XiC=Xiv.*Xiv./XiT;
%         
%         BAT=Tspec(Tdis,XiC,KB,K,D,q);
%         MAD = meanabsdev( PSD(ikcor), BAT(ikcor), Sn(ikcor) );
%     end
%     
%     %fits parameters in the noise free region
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     if method == 'B'
%         KF = K(iKnM0); %iKn
%         Ktest = linspace(max([5*K1,KF/5]),KF*5,40); %5*KF
%         cost = nan(size(Ktest));
%         Xif0 = nan(size(Ktest));
%         dK = Ktest(2)-Ktest(1);
%         for i = 1:length(Ktest)
%             Ks = ksfact*Pr^(-0.5)*Ktest(i);
%             Kf1 = max([Ks,K1]);
%             iKf1 = find(K>=Kf1,1,'first');
%             ikfit = iKf1:iKnM0;
%             [cost(i), Xif0(i),Steo] =  cost_T_fit(K, PSD, Ktest(i), ikfit, Tdis, q, dof, Sn,D,Snfact);
%         end
%         
%         LKHtest = - cost;
%         ML = max(LKHtest);
%         iML = find(LKHtest == ML);
%         KB00 = Ktest(iML);
%         BATf0 = Tspec(Tdis,Xif0(iML),KB00,K,D,q);
%         Kp=KB00/sqrt(6*q); % K peak
%         
%         %looks for Luketina low-wavenumber cut-off
%         if strcmp(int_range,'S')
%             KL=K1; islope=2;
%         elseif strcmp(int_range,'L')
%             % Find the intersection between observed and theoretical spectrum
%             diffS = (PSD-BATf0)./abs(PSD-BATf0);
%             iKL = find(diffS<0,1,'first');
%             if ~isempty(iKL) && abs(PSD(iKL)-BATf0(iKL))>abs(PSD(iKL-1)-BATf0(iKL-1))
%                 iKL=iKL-1;
%             end
%             % This was not mentioned in Luketina, it is a double check on
%             % the existence of an actual finestructure pattern.
%             % Calculate the slope of the first part of the spectrum using an increasing window
%             % Note: K(1) and PSD(1) = 0 --> the window starts from the 2nd point.
%             % The iteration starts from the 4th point to have at least 3 points to calculate the slope
%             iKp = find(K<=Kp,1,'last');
%             slope=NaN(1,iKnM0);
%             for i=4:iKnM0
%                 X=log10(K(2:i)); Y=log10(PSD(2:i));
%                 X=[ones(length(X),1) (X)];
%                 b=X\Y;     % \ computes the least square fit
%                 slope(i)=b(2);
%             end
%             if ~isempty(iKL)                
%                 [vslope,islope]=nanmin(slope(1:min(iKp,iKL))); % find the minimum of the slope
%                 if vslope>0 || nanmean(slope(1:islope))>0 % if the minimum slope is > 0, then
%                     islope=2;
%                 end
%             else
%                 islope=2;
%             end
%             KL = min([K(islope),K(iKL),Kp]);
%         end
%         
%         if  plt(1)~=0 % Before running the 2nd round, plot the results of the first
%             Sbest=Tspec(Tdis,Xif0(iML),KB00,K0,D,q); % just for the plot
%             Ks = ksfact*Pr^(-0.5)*KB00;
%             Kf1 = max([Ks,K1]);
%             iKf1 = find(K>=Kf1,1,'first');
%             ikfit = iKf1:iKnM0;
%             
%             fh=figure(10);
%             set(fh,'Units', 'Inches', 'Position', [0,0,3,6],...
%             'PaperUnits', 'Inches', 'PaperSize', [3,6]);
%             subplot(2,1,1)
%             loglog(K0,PSD0,'-k'); hold on
% %             loglog(K,Sn+BAT,'color','k', 'linewidth',2);
%             loglog(K0,Sn_ODAS+Sbest,'-b', 'linewidth',2)
%             loglog(K0,Sn_ODAS.*H0+PSD0.*H0,'-','color',[0.5 0.5 0.5])
%             loglog(K(ikfit),PSD(ikfit),'ob','markersize',4);
%             loglog(K(islope),PSD(islope),'sr')
%             grid on
%             hold off
%             xlabel('K (cpm)')
%             ylabel('PSD (K^2 m^{-2} cpm^{-1})')
%             title(['z= ', num2str(mpres), ' m, sensor= ', T_string])
%             xlim([1,1000]); ylim([1e-9,1e3]);
%         end
%         
%         %second rough search
%         cost = nan(size(Ktest));
%         Xif0 = nan(size(Ktest));
%         dK = Ktest(2)-Ktest(1);
%         for i = 1:length(Ktest)
%             Ks = ksfact*Pr^(-0.5)*Ktest(i);
%             Kf1 = max([Ks,KL]);    % Note, if int_range=='L' here we consider KL (Luketina & Imberger, 2001)
%             iKf1 = find(K>=Kf1,1,'first');
%             ikfit = iKf1:iKnM0;
%             [cost(i), Xif0(i)] =  cost_T_fit(K, PSD, Ktest(i), ikfit, Tdis, q, dof, Sn,D,Snfact);
%         end
%         LKHtest = - cost;
%         ML = max(LKHtest);
%         iML = find(LKHtest == ML);
%         KB00 = Ktest(iML);
%         BATf0 = Tspec(Tdis,Xif0(iML),KB00,K,D,q);
%         Kp=KB00/sqrt(6*q); % K peak
%         
%         %looks for Luketina low-wavenumber cut-off
%         if strcmp(int_range,'S')
%             KL=K1;
%         elseif strcmp(int_range,'L')
%            
%             diffS = (PSD-BATf0)./abs(PSD-BATf0);
%             iKL = find(diffS<0,1,'first');
%             iKp = find(K<=Kp,1,'last');
%             if ~isempty(iKL) && abs(PSD(iKL)-BATf0(iKL))>abs(PSD(iKL-1)-BATf0(iKL-1))
%                 iKL=iKL-1;
%             end
%             if ~isempty(iKL)
%                 [vslope,islope]=nanmin(slope(1:min(iKp,iKL))); % find the minimum of the slope
%                 if vslope>0 || nanmean(slope(1:islope))>0 % if the minimum slope is > 0, then
%                     islope=2;
%                 end
%             else
%                 islope=2;
%             end
%             KL = min([K(islope),K(iKL),Kp]);
%         end
%         iKL = find(K>=KL,1,'first');
% 
%         %refines the search with the second derivative
%         MLp = -cost_T_fit(K, PSD, KB00+dK, ikfit, Tdis, q, dof, Sn,D,Snfact);
%         MLm = -cost_T_fit(K, PSD, KB00-dK, ikfit, Tdis, q, dof, Sn,D,Snfact);
%         deltak0=abs((2*dK)/sqrt(2*ML-MLm-MLp));
%         deltak = max([deltak0,dK]);
%         
%         kmin2 = max([KB00-deltak,K(iKL+1)]);
%         kmax2 = KB00+deltak;%min([KB00+deltak, (K(iKn)-3*dK)/0.04*Pr^0.5 ]);
%         Ktest = linspace(kmin2,kmax2,40);
%         clear cost Xif0
%         cost = nan(size(Ktest));
%         Xif0 = nan(size(Ktest));
%         %         Kf1 = max([KL,K1]);
%         %         iKf1 = find(K>=Kf1,1,'first');
%         %         ikfit = iKf1:iKnM0;
%         for i = 1:length(Ktest)
%             Ks = ksfact*Pr^(-0.5)*Ktest(i);
%             Kf1 = max([Ks,KL]);    % Note, if int_range=='L' here we consider KL (Luketina & Imberger, 2001)
%             iKf1 = find(K>=Kf1,1,'first');
%             ikfit = iKf1:iKnM0;
%             [cost(i), Xif0(i)] =  cost_T_fit(K, PSD, Ktest(i), ikfit, Tdis, q, dof, Sn,D,Snfact);
%         end
%         
%         LKHtest = - cost;
%         MLKH = max(LKHtest);
%         iML = find(LKHtest == MLKH);
%         Xif = Xif0(iML);
%         KBT = Ktest(iML);
%         Kp=KBT/sqrt(6*q); % K peak
%         
%         BATf = Tspec(Tdis,Xif,KBT,K,D,q);
%         MADf = meanabsdev( PSD(ikfit), BATf(ikfit), Sn(ikfit) );
%         % re-compute the integration range (for figure and fit_flag)
%         Ks = ksfact*Pr^(-0.5)*KBT;
%         Kf1 = max([Ks,KL]);
%         iKf1 = find(K>=Kf1,1,'first');
%         ikfit = iKf1:iKnM0;
%                 
%         %calculates uncertainties in the fitting parameters
%         [MLp,Xip,~] = cost_T_fit(K, PSD, KBT+deltak, ikfit, Tdis, q, dof, Sn,D,Snfact);
%         [MLm,Xim,~] = cost_T_fit(K, PSD, KBT-deltak, ikfit, Tdis, q, dof, Sn,D,Snfact);
%         MLp = - MLp;
%         MLm = - MLm;
%         sKBT=abs((2*deltak)/sqrt(2*MLKH-MLm-MLp));
%         sXif = abs((Xip-Xim)/sqrt(2*MLKH-MLm-MLp));
%         
%         %fits to polynom (avoiding the noisy part)
%         ikfitA = ikfit(1):min([ikfit(end),find(BATf<Snfact*Sn,1,'first')]);
%         logK = log(K(ikfitA));
%         logS = log(PSD(ikfitA));
%         pp=polyfit(logK, logS,1);
%         Sm = exp(polyval(pp, log(K)));
%         LKHpol = -cost_MLE(PSD(ikfit), Sm(ikfit), dof, Sn(ikfit));
%         
%         LKHratio = MLKH - LKHpol;
%         LKHratio = log10(exp(1))*LKHratio;
%         %likelyhood ratio is log10(Pteo/Pexp) and C = log(Pteo), to converto
%         maxK=K(iKnM0);
%         rangeK = log10(K(iKnM0)) - log10(K(iKf1));
%         if LKHratio>2 && MADf<MADc*2 ... % && KBT<maxK ...
%                 && abs(sKBT)<0.5*abs(KBT)  && rangeK> 0.8 ...
%                 && K(iKnM0)>2*Kp && Kf1<Kp
%             %quality flag for the fit
%             fit_flag =1;
%         end
        
%         if  plt(1)~=0
%             figure(10)
%             subplot(2,1,2)
%             Sbest=Tspec(Tdis,Xif,KBT,K0,D,q);
%             
%             loglog(K0,PSD0,'-k'); hold on
% %             loglog(K,Sn+BAT,'color','k', 'linewidth',2);
%             loglog(K0,Sn_ODAS+Sbest,'-b', 'linewidth',2)
%             loglog(K0,Sn_ODAS.*H0+PSD0.*H0,'-','color',[0.5 0.5 0.5])
%             loglog(K(ikfit),PSD(ikfit),'ob','markersize',4);
%             loglog(K(islope),PSD(islope),'sr')
%             grid on
%             xlabel('K (cpm)')
%             ylabel('PSD (K^2 m^{-2} cpm^{-1})')
%             title(['z= ', num2str(mpres), ' m, sensor= ', T_string])
%             xlim([1,1000]); ylim([1e-9,1e3]);
%             text(1.2,10^1.5,['\epsilon_{T} = ', num2str(visco*D^2*(2*pi()*KBT)^4,'%1.3e'),' m^2/s^{-3}'],'horizontalalignment','left','fontsize',7)
%             text(1.2,10^0.5,['\epsilon_{sh} = ', num2str(visco*D^2*(2*pi()*KB)^4,'%1.3e'),' m^2/s^{-3}'],'horizontalalignment','left','fontsize',7)
%             text(1.2,10^-5,['Quality flag = ', num2str(fit_flag)],'horizontalalignment','left','fontsize',7)
%             text(1.2,10^-6, ['MADf = ', num2str(MADf, '%1.2f'),'(<', num2str(2*MADc, '%1.2f'),')'],'horizontalalignment','left','fontsize',7)
%             text(1.2,10^-7, ['KBT = ', num2str(KBT, '%1.0f'),'(<',num2str(maxK, '%1.0f'),')'],'horizontalalignment','left','fontsize',7)
%             text(1.2,10^-8, ['LKHratio = ', num2str(LKHratio, '%1.1f'),'(>2)'],'horizontalalignment','left','fontsize',7)
%             saveas(gcf,[folder_out '/Spec_z= ', num2str(mpres), ' m, sensor= ', T_string, '.png'])
%             hold off
%         end
%         
%     else
%         %EPFL method
%         [LKHratio,MLKH,KBT,~,QUAL]=Fit_kB_OSS_Carpenter(K,PSD,Sn,1);
%         Xif = QUAL.X_T1;
%         BATf = Tspec(Tdis,Xif,KBT,K,D,q);
%         MADf = QUAL.MAD1;
%         ikfit = iK1:iKn;
%     end
%     
%     
%     %plots
%     iKBT=find(K0>=KBT,1,'first');if isempty(iKBT), iKBT=length(K); end
    
    if plt(2)~=0
        figure(11)
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
        xlabel('gradC (mS cm^{-1} m^{-1})')
        ylabel('p (db)')
        
        axes('position',[0.4,0.12,0.5,0.8])
%         loglog(K0,PSD0, 'color' ,[0.5,0.5,0.5])
%         hold on
        loglog(K,PSD,'k')
%         %loglog(K2,PSD2,'-b')
%         if cont
%             loglog(K(ikfit),PSD(ikfit),'ko', 'markersize',4)
%             loglog(K(ikfit),Sm(ikfit)+Sn(ikfit),'b')
%             
%             loglog(K,Sn+BAT,'color','k', 'linewidth',2);
%             loglog(K,Sn+BATf, 'color',[0.5,0.5,0.5], 'linewidth',2);
%             loglog(K0,Tspec(Tdis, Xif, KBT,K0,D,q), 'color',[0.5,0.5,0.5], 'linewidth',1,'linestyle','--');
%             loglog(K0,Tspec(Tdis, XiC, KB, K0,D,q),'color','k','linewidth',1,'linestyle','--');
%         end
%         loglog(K0,Sn0, 'color','r', 'linewidth',1,'linestyle','-');
%         ylim([10^-9,1e3])
%         line([K0(iKn) K0(iKn)], ylim,'color','k','linestyle','--')
%         line([K0(iK1) K0(iK1)], ylim,'color','k','linestyle','--')
%         line([KBT KBT], ylim,'color','k')
%         text(1.1*KBT,0.1,'K_B')
%         if Tdis=='B'
%             text(1.2,10^2.5,'Batchelor Spectrum','horizontalalignment','left')
%         elseif Tdis=='K'
%             text(1.2,10^2.5,'Kraichnan Spectrum','horizontalalignment','left')
%         end
%         
%         text(1.2,10^1.5,['\chi_{var. cor.} = ', num2str(XiC,'%1.3e'),' K^2/s'],'horizontalalignment','left')
%         text(1.2,10^0.75,['\chi_{fit} = ', num2str(Xif,'%1.3e'),'\pm',num2str(sXif,'%1.3e'),' K^2/s'],'horizontalalignment','left')
%         text(1.2,10^0,['\chi_{var} = ', num2str(Xiv,'%1.3e'),' K^2/s'],'horizontalalignment','left')
%         text(1.2,10^-0.75,['K_B^{sh}= ', num2str(KB,'%1.0f'),' cpm'],'horizontalalignment','left')
%         text(1.2,10^-1.5,['K_{B}= ', num2str(KBT,'%1.0f'),'\pm',num2str(sKBT,'%1.0f'),' cpm'],'horizontalalignment','left')
%         text(100, 1e-7, ['MADf = ', num2str(MADf, '%1.2f'),'(', num2str(2*MADc, '%1.2f'),')'],'horizontalalignment','left')
%         text(100, 1e-8, ['MAD = ', num2str(MAD, '%1.2f')],'horizontalalignment','left')
%         text(100, 1e-6, ['LKHratio = ', num2str(LKHratio, '%1.1f')],'horizontalalignment','left')
        title([num2str(mpres), 'm'],'Fontsize',12)
        xlabel('K (cpm)')
        ylabel('PSD (mS^2 cm^{-2} m^{-2} cpm^{-1})')
        xlim([1,1000])
        saveas(gcf,[folder_out '/Spec_z= ', num2str(mpres), ' m, sensor=C1.png'])
        fit_flag
%         loglog(K0,Snold,'-g')
        %saveas(gcf,'Tspec.png')
%         pause()
        
        
    end
    Xiv=NaN;
    XiC=NaN;
    Xif=NaN;
    KBT=NaN;
    fit_flag=NaN;
    sXif=NaN;
    sKBT=NaN;
    MAD=NaN;
    MADf=NaN;
    MADc=NaN;
    MLKH=NaN;
    LKHratio=NaN;
    maxK=NaN;
    rangeK=NaN;
    
end
