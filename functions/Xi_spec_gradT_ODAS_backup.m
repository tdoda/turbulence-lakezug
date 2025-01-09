function [Xiv,XiC,Xif,KBT,fit_flag, sXif,sKBT,MAD, MADf,MADc,MLKH,LKHratio,maxK]=Xi_spec_gradT(pres,x0,K1,fn,KB,W,noisep, sL, sOV,Tdis,q,method,time_res,peak_rem,plt,time,D,visco,T_dT,T_string,setupstr,int_range)
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
%peak_rem: [fmin,fmax]: to frequencies to remove a noise peak (for T1221,T1228)
%plt: if !=0 shows the spectra
%time: time vector
%D: thermal diffusivity
%visco: viscosity
%T_dT: raw pre-emphasized signal
%T_string: name of the channel 'T1_dT1' or 'T2_dT2'
%setupstr: configuration file
%int_range: Integration range according to Steinbuck 2009 ('S') or Luketina and Imberger 2001 ('L')

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

scalar_info.fft_length      = sL;
scalar_info.spec_length     = length(x);
scalar_info.overlap         = length(x)/2;
scalar_info.fs              = 512;
scalar_info.gradient_method = 'high_pass';
scalar_info.f_AA            = 98;

MADc = sqrt(2/dof);
%     sp = get_scalar_spectra_odas(x, [], pres, time', W, scalar_info);
%     fr = sp.F; K=sp.K; PSDT=sp.scalar_spec;

[PSDT,fr] = csd_odas(x-nanmean(x),x-nanmean(x),sL,scalar_info.fs,[],sOV,'linear');
PSDT=PSDT*W;
K = fr/W;

ii = find(fr>peak_rem(1) & fr<peak_rem(2));
fr(ii)=[];
K(ii)=[];
PSDT(ii) = [];
PSD=PSDT;

%correction variance Goto2016/ tau = 10ms
tau = 0.002*W.^-0.32;
%tau = 0.010;
%tau=4.1E-3*W^(-0.5);
%correction varianza microCTD

    F0 = 25*sqrt(W);
    tau = (2*pi()*F0/sqrt(sqrt(2)-1))^(-1);
%     if ~isnan(time_res)
%         tau = time_res;
%     end
H = 1./(1 + (2*pi()*tau*fr).^2).^2;
PSD = PSD./H;
%iKmax = find(H>0.1,1,'last');

%noise function
noise_info=gradT_noise_odas;
noise_info.gamma_RSI = 1;
noise_info.E_n = 3e-9;
Sn_ODAS = gradT_noise_odas(T_dT, T_string, W, fr, setupstr,noise_info);
Sn_ODAS = Sn_ODAS*W;
Sn_ODAS = Sn_ODAS./H;
Sn_ODAS(1) = 0;

Snold = FP07noise(noisep,fr); %for microCTD with up and down!
Snold = Snold*W;
Snold = Snold.*(2*pi()*K).^2;
Snold = Snold./H;
Snold(1) = 0;
%     figure(10)
%     loglog(K,PSD); hold on
%     loglog(K,Sn_ODAS)
%     loglog(K,Sn)

Sn = Sn_ODAS;

mpres=mean(pres);
Kn = fn/W;
iK1=find(K>K1,1,'first');
iKn = find(K<=Kn,1,'last');
%iKn=min([find(K<=Kn,1,'last'),iKmax]); %if I want to remove
%strong correction
iKB=find(K>=KB,1,'first'); if isempty(iKB), iKB=length(K); end

%deletes undesired part of the spectrum
K0 = K;
PSD0 = PSD;
Sn0 = Sn;
H0 = H;

K = K(1:iKn);
maxK = K(end);
PSD = PSD(1:iKn);
H = H(1:iKn);
Sn = Sn(1:iKn);

%variance in the noise free part determined from the noise
%model
ksfact=0.04; % (0.04 in Steinbuck 2009: https://journals.ametsoc.org/doi/pdf/10.1175/2009JTECHO611.1)
Snfact=2; %Snfact=1.55;   (1.55 as in Goto 2016 DOI: 10.1175/JTECH-D-15-0220.1)
% iKnM0 = find(PSD<Snfact*Sn | H<0.01 | K==max(K),1,'first');
iKnM0 = find(PSD<Snfact*Sn,1,'first');

% Calculate the slope of the first part of the spectrum using an increasing window
% Note: K(1) and PSD(1) = 0 --> the window starts from the 2nd point. 
% The iteration starts from the 4th point to have at least 3 points to calculate the slope
for i=4:iKnM0  
    X=log10(K(2:i)); Y=log10(PSD(2:i));
    X=[ones(length(X),1) (X)];
    b=X\Y;     % \ computes the least square fit
    slope(i)=b(2);
end
[vslope,islope]=min(slope); % find the minimum of the slope
if vslope>0  % if the minimum slope is > 0, then 
    islope=1;
end

ks = ksfact*Pr^(-0.5)*KB;
K11 = max([ks,K1]);
iK11 = find(K>=K11,1,'first'); %if I do like that the fit fails for high KB
ikcor = iK11:iKnM0;
Xiv = 6*D*sum( (PSD(ikcor(2:end))+PSD(ikcor(1:end-1))).* (K(ikcor(2:end))-K(ikcor(1:end-1))) )/2;
Xin = 6*D*sum( (Sn(ikcor(2:end))+Sn(ikcor(1:end-1))).* (K(ikcor(2:end))-K(ikcor(1:end-1))) )/2;
cont = false;
if Xiv>1.3*Xin
    cont = true;
    Xiv = Xiv -Xin;
else
    Xiv = 0;
    Xif = 0;
end


%continues only if Xiv is detectable
if cont
    
    %calculates Xi using theoretical spectrum to correct for unresolved
    if KB>0
        BAT=Tspec(Tdis,Xiv,KB,K,D,q);
        XiT=6*D*sum( (BAT(ikcor(2:end))+BAT(ikcor(1:end-1))).*(K(ikcor(2:end))-K(ikcor(1:end-1))) )/2;
        XiC=Xiv.*Xiv./XiT;
            
        BAT=Tspec(Tdis,XiC,KB,K,D,q);
        MAD = meanabsdev( PSD(ikcor), BAT(ikcor), Sn(ikcor) );
    end
    
    %fits parameters in the noise free region
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if method == 'B'
        KF = K(iKnM0); %iKn
        Ktest = linspace(max([5*K1,KF/5]),KF*5,40); %5*KF
        cost = nan(size(Ktest));
        Xif0 = nan(size(Ktest));
        dK = Ktest(2)-Ktest(1);
        Ktest=Ktest;
        best=1e10; %large number, initialization
        for i = 1:length(Ktest)
            ks = ksfact*Pr^(-0.5)*Ktest(i);
            if strcmp(int_range,'S') % Steinerbuck
                Kf1 = max([ks,K1]);
            else % Luketina and Imberger 2001
                %Prior to the first iteration, ?L is set to the lowest wavenumber of the
                %measured spectra D(?). For each subsequent set of trials, ?L is the lowest
                %wavenumber at which the measured spectra D(?) and the best-fit theoretical
                %spectra S(?) (from the set of fits done within the preceding iteration) intersect
                if i==1
                    Kf1 = max([ks,K1]);
                else
                    if (cost(i-1)<best)
                        best=cost(i-1);
                        idxs=find(PSD<Steo,1,'first');
                        if ~isempty(idxs) && abs(PSD(idxs)-Steo(idxs))>abs(PSD(idxs-1)-Steo(idxs-1))
                            idxs=idxs-1;
                        end
                    end
                    
                    
                    
                    
% Calculate the slope of the first part of the spectrum using an increasing window
% Note: K(1) and PSD(1) = 0 --> the window starts from the 2nd point
for i=3:iKnM0  
    X=log10(K(2:i)); Y=log10(PSD(2:i));
    X=[ones(length(X),1) (X)];
    b=X\Y;     % \ computes the least square fit
    slope(i)=b(2);
end
[vslope,islope]=min(slope); % find the minimum of the slope
if vslope>0  % if the minimum slope is > 0, then 
    islope=1;
end



                    if ~isempty(idxs)
                        kp = Ktest(i)/sqrt(6*q);          % Kpeak (See e.g., Goto et al., 2016)
                        Kf1=min( max(K(idxs),ks) ,kp);
                    else
                        Kf1=max([ks,K1]);
                    end
                end
            end
            iKf1 = find(K>=Kf1,1,'first');
            ikfit = iKf1:iKnM0;
            [cost(i), Xif0(i),Steo] =  cost_T_fit(K, PSD, Ktest(i), ikfit, Tdis, q, dof, Sn,D,Snfact);
            
            if plt(1)~=0
                figure(10)
                title(['z= ', num2str(mpres), ' m, sensor= ', T_string])
                subplot(2,1,1)
%                 loglog(K,Sn+Steo,'-','color',[0.5 0.5 0.5], 'linewidth',2); hold on
%                 loglog(K,Sn+BAT,'color','k', 'linewidth',2);
%                 loglog(K0,PSD0,'-k');
%                 loglog(K(ikfit),PSD(ikfit),'ok');
%                 loglog(K0,Sn_ODAS,'-r')
%                 loglog(K0,Snold,'-g')
%                 hold off
            end
        end
       
        LKHtest = - cost;
        ML = max(LKHtest);
        iML = find(LKHtest == ML);
        KB00 = Ktest(iML);
        
        ks = ksfact*Pr^(-0.5)*KB00;        
        if strcmp(int_range,'S') % Steinerbuck
            Kf1 = max([ks,K1]);
        elseif strcmp(int_range,'L') % Luketina and Imberger 2001
            if ~isempty(idxs)
                kp = KB00/sqrt(6*q);          % Kpeak (See e.g., Goto et al., 2016)
                Kf1=min( max(K(idxs),ks) ,kp);
            else
                Kf1=max([ks,K1]);
            end
        end
        iKf1 = find(K>=Kf1,1,'first');
        ikfit = iKf1:iKnM0;

        if plt(1)~=0
            Sbest=Tspec(Tdis,Xif0(iML),KB00,K0,D,q);
            loglog(K0,PSD0,'-k'); hold on
            loglog(K,Sn+BAT,'color','k', 'linewidth',2);
            loglog(K0,Sn_ODAS+Sbest,'-b', 'linewidth',2)
            loglog(K0,Sn_ODAS.*H0+PSD0.*H0,':k')
            loglog(K(ikfit),PSD(ikfit),'ob','markersize',4);
            hold off
        end
        
        MLp = -cost_T_fit(K, PSD, KB00+dK, ikfit, Tdis, q, dof, Sn,D,Snfact);
        MLm = -cost_T_fit(K, PSD, KB00-dK, ikfit, Tdis, q, dof, Sn,D,Snfact);
        deltak0=abs((2*dK)/sqrt(2*ML-MLm-MLp));
        deltak = max([deltak0,dK]);
        
        kmin2 = max([KB00-deltak,K(iK1+1)]);
        kmax2 = KB00+deltak;%min([KB00+deltak, (K(iKn)-3*dK)/0.04*Pr^0.5 ]);
        Ktest = linspace(kmin2,kmax2,40);
        clear cost Xif0
        cost = nan(size(Ktest));
        Xif0 = nan(size(Ktest));
        best=1e10;
        for i = 1:length(Ktest)
            ks = ksfact*Pr^(-0.5)*Ktest(i);
            if strcmp(int_range,'S') % Steinerbuck
                Kf1 = max([ks,K1]);
            else % Luketina and Imberger 2001
                %Prior to the first iteration, ?L is set to the lowest wavenumber of the
                %measured spectra D(?). For each subsequent set of trials, ?L is the lowest
                %wavenumber at which the measured spectra D(?) and the best-fit theoretical
                %spectra S(?) (from the set of fits done within the preceding iteration) intersect
                if i==1
                    Kf1 = max([ks,K1]);
                else
                    if (cost(i-1)<best)
                        best=cost(i-1);
                        idxs=find(PSD<Steo,1,'first');
                        if ~isempty(idxs) && abs(PSD(idxs)-Steo(idxs))>abs(PSD(idxs-1)-Steo(idxs-1))
                            idxs=idxs-1;
                        end
                    end
                    if ~isempty(idxs)
                        kp = Ktest(i)/sqrt(6*q);          % Kpeak (See e.g., Goto et al., 2016)
                        Kf1=min( max(K(idxs),ks) ,kp);
                    else
                        Kf1=max([ks,K1]);
                    end
                end
            end
            iKf1 = find(K>=Kf1,1,'first');
            ikfit = iKf1:iKnM0;
            [cost(i), Xif0(i),Steo] =  cost_T_fit(K, PSD, Ktest(i), ikfit, Tdis, q, dof, Sn,D,Snfact);
            
            if plt(1)~=0
                figure(10)
                subplot(2,1,2)
%                 loglog(K,Sn+Steo,'-','color',[0.5 0.5 0.5], 'linewidth',2); hold on
%                 loglog(K,Sn+BAT,'color','k', 'linewidth',2);
%                 loglog(K0,PSD0,'-k');
%                 loglog(K(ikfit),PSD(ikfit),'ok');
%                 loglog(K0,Sn_ODAS,'-r')
%                 loglog(K0,Snold,'-g')
%                 hold off
            end
        end
        
        LKHtest = - cost;
        MLKH = max(LKHtest);
        iML = find(LKHtest == MLKH);
        Xif = Xif0(iML);
        KBT = Ktest(iML);
        
        BATf = Tspec(Tdis,Xif,KBT,K,D,q);
        ks = ksfact*Pr^(-0.5)*KB00;        
        if strcmp(int_range,'S') % Steinerbuck
            Kf1 = max([ks,K1]);
        elseif strcmp(int_range,'L') % Luketina and Imberger 2001
            if ~isempty(idxs)
                kp = KB00/sqrt(6*q);          % Kpeak (See e.g., Goto et al., 2016)
                Kf1=min( max(K(idxs),ks) ,kp);
            else
                Kf1=max([ks,K1]);
            end
        end
        iKf1 = find(K>=Kf1,1,'first');
        ikfit = iKf1:iKnM0;

        if plt(1)~=0
            Sbest=Tspec(Tdis,Xif,KBT,K0,D,q);
            loglog(K0,PSD0,'-k'); hold on
            loglog(K,Sn+BAT,'color','k', 'linewidth',2);
            loglog(K0,Sn_ODAS+Sbest,'-b', 'linewidth',2)
            loglog(K0,Sn_ODAS.*H0+PSD0.*H0,':k')
            loglog(K(ikfit),PSD(ikfit),'ob','markersize',4);
            hold off
            saveas(gcf,['Spec_z= ', num2str(mpres), ' m, sensor= ', T_string, '.png'])
        end
                    
        iMADf = meanabsdev( PSD(ikfit), BATf(ikfit), Sn(ikfit) );
        
        %calculates uncertainties in the fitting parameters
        [MLp,Xip,~] = cost_T_fit(K, PSD, KBT+deltak, ikfit, Tdis, q, dof, Sn,D,Snfact);
        [MLm,Xim,~] = cost_T_fit(K, PSD, KBT-deltak, ikfit, Tdis, q, dof, Sn,D,Snfact);
        MLp = - MLp;
        MLm = - MLm;
        sKBT=abs((2*deltak)/sqrt(2*MLKH-MLm-MLp));
        sXif = abs((Xip-Xim)/sqrt(2*MLKH-MLm-MLp));
        
        %fits to polynom (avoiding the noisy part)
        ikfitA = ikfit(1):min([ikfit(end),find(BATf<Snfact*Sn,1,'first')]);
        logK = log(K(ikfitA));
        logS = log(PSD(ikfitA));
        pp=polyfit(logK, logS,1);
        Sm = exp(polyval(pp, log(K)));
        LKHpol = -cost_MLE(PSD(ikfit), Sm(ikfit), dof, Sn(ikfit));
        
        LKHratio = MLKH - LKHpol;
        LKHratio = log10(exp(1))*LKHratio;
        %likelyhood ratio is log10(Pteo/Pexp) and C = log(Pteo), to converto
        
        if LKHratio>2 && MADf<MADc*2 && KBT<maxK/1.66 && abs(sKBT)<0.5*abs(KBT)
            %quality flag for the fit
            fit_flag =1;
        end
    else
        %EPFL method
        [LKHratio,MLKH,KBT,~,QUAL]=Fit_kB_OSS_Carpenter(K,PSD,Sn,1);
        Xif = QUAL.X_T1;
        BATf = Tspec(Tdis,Xif,KBT,K,D,q);
        MADf = QUAL.MAD1;
        ikfit = iK1:iKn;
    end
    
    
    %plots
    iKBT=find(K0>=KBT,1,'first');if isempty(iKBT), iKBT=length(K); end
    
    if plt(2)~=0
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
            
            loglog(K,Sn+BAT,'color','k', 'linewidth',2);
            loglog(K,Sn+BATf, 'color',[0.5,0.5,0.5], 'linewidth',2);
            loglog(K0,Tspec(Tdis, Xif, KBT,K0,D,q), 'color',[0.5,0.5,0.5], 'linewidth',1,'linestyle','--');
            loglog(K0,Tspec(Tdis, XiC, KB, K0,D,q),'color','k','linewidth',1,'linestyle','--');
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
        
        text(1.2,10^1.5,['\chi_{var. cor.} = ', num2str(XiC,'%1.3e'),' K^2/s'],'horizontalalignment','left')
        text(1.2,10^0.75,['\chi_{fit} = ', num2str(Xif,'%1.3e'),'\pm',num2str(sXif,'%1.3e'),' K^2/s'],'horizontalalignment','left')
        text(1.2,10^0,['\chi_{var} = ', num2str(Xiv,'%1.3e'),' K^2/s'],'horizontalalignment','left')
        text(1.2,10^-0.75,['K_B^{sh}= ', num2str(KB,'%1.0f'),' cpm'],'horizontalalignment','left')
        text(1.2,10^-1.5,['K_{B}= ', num2str(KBT,'%1.0f'),'\pm',num2str(sKBT,'%1.0f'),' cpm'],'horizontalalignment','left')
        text(100, 1e-7, ['MADf = ', num2str(MADf, '%1.2f'),'(', num2str(2*MADc, '%1.2f'),')'],'horizontalalignment','left')
        text(100, 1e-8, ['MAD = ', num2str(MAD, '%1.2f')],'horizontalalignment','left')
        text(100, 1e-6, ['LKHratio = ', num2str(LKHratio, '%1.1f')],'horizontalalignment','left')
        title([num2str(mpres), 'm'],'Fontsize',12)
        xlabel('K (cpm)')
        ylabel('PSD (K^2 m^{-2} cpm^{-1})')
        xlim([1,1000])
        
        fit_flag
        loglog(K0,Snold,'-g')
        %saveas(gcf,'Tspec.png')
        pause()
        
        
    end
    
end
