function [Xiv,Xi_ST,Xi_T,kB,eps_T,MAD_ST,MAD_T,MADc,LR,kL,kU,krange,kP,fit_flag_T,SPECTRUM]=gradT_dis_spec(pres,x0,k1,fn,kB_S,W,sL,sOV,Tdis,q,tau_0,time_corr,npoles,int_range,D,visco,T_dT,T_string,setupstr,plt,presplot,Tplot,ksfact, Snfact)

% GOAL
% Estimate of TKE and temperature variance dissipation rates eps and Xi by 
% integration of the measured temperature gradient wavenumber spectrum

% INPUT:
% pres: pressure vector (db)
% x0: grad temperature vector (°C/m)
% k1: minimium wavenumber for integration (small number, e.g., 0.1 cpm)
% fn: maximum frequency for calculations (90% of anti-aliasing filter f_AA. E.g., if f_AA=98 Hz, fn=0.9*98=88.2 Hz)
% kB_S: Batchelor wavenumber determined from shear probe used to calculate Xi_ST (if 0 does not calculate. When two sh are available and accepted, the average of kB_S is used.)
% W: mean profiling speed (m/s)
% sL: length of segments for fft (scans)
% sOV: overlap for fft (scans)
% Tdis: type of theoretical spectrum (B: Batchelor, K: Kraichnan)
% q: turbulent parameter
% tau_0: nominal response time (s)
% time_corr: time correction approach: KOC, RSI, NAS, SOM (see below for details)
% npoles: transfer function for time response correction 'single' or 'double' pole pole 
% plt: if !=0 shows the spectra
% int_range: lower limit of the integration range according to Steinbuck et al 2009 ('S') or Luketina and Imberger 2001 ('L')
% D: thermal diffusivity (m^2/s)
% visco: kinematic viscosity of water (m^2/s)
% T_dT (*): raw pre-emphasized temperature signal
% T_string (*): name of the FP07 channel 'T1_dT1' or 'T2_dT2'
% setupstr (*): setupfilestr configuration file provided by ODAS libraries
% plt: flag for making the figure (0 no figure, ~=0 make figure)
% presplot: pressure vector for making the plot
% Tplot: temperature vector for making the plot
% ksfact (optional, defaut: 0.04): upper bound of the inertial-convective
% subrange, as in Steinbuck et al. (2009) [added by T. Doda]
% Snfact (optional, default: 1.55): maximum acceptable signal to noise
% ratio 1.55, as in Goto et al. (2016) [added by T. Doda]

% (*) these entries are required for calculating the noise spectrum according 
% to the function profided by RSI in the ODAS libraries.
% They are obtained reading the .p files with the ODAS libraries. If the
% ODAS libraries are not accessible by the user, an alternative for
% caluclating the noise spectrum is suggested below.

% OUTPUT
% Xiv: Xi from spectral integration in the well resolved part (°C^2/s)
% Xi_ST: Xi corrected with kB_S obtained from shear probes (°C^2/s)
% Xi_T: Xi after MLE spectral fitting  (°C^2/s)
% kB: Batchelor wavenumber after MLE spectral fitting (cpm) 
% eps_T: TKE dissipation rate after MLE spectral fitting (m^2/s^3)
% MAD_ST: Mean Absolute Deviation between observed and empirical spectra using kB_S from shear
% MAD_T: Mean Absolute Deviation between observed and empirical spectra using MLE fitting
% MADc: Threshold for the Mean Absolute Deviation between observed and empirical spectra
% LR: likelihood ratio
% kL: Lower integration wavenumber (cpm)
% kU: Upper integration wavenumber (cpm)
% krange: wavenumber range used for spectral integration 
% kP: wavenumber corresponding to the peak of the fitted theoretical spectrum (cpm)
% fit_flag_ST: acceptance flag of Xi_ST according to quality metrics: 0 rejected, 1 accepted
% fit_flag_T: acceptance flag of eps_T and Xi_T according to quality metrics: 0 rejected, 1 accepted
% SPECTRUM: structure array containing all information to plot the spectra
% [added by T. Doda, 06.02.25]

%% 
if nargin<23
    ksfact=0.04;
end

if nargin<24
    Snfact=1.55;
end
%% Initialization
Xiv = NaN; Xi_ST=NaN; Xi_T=NaN; kB=NaN; sXif = nan; sKBT = nan; 
MAD_ST = NaN; MAD_T = NaN;MLKH=nan; LR = nan; kU = nan; 
fit_flag_ST = NaN; fit_flag_T = NaN;
Pr = visco/D; % Prandtl number

%% Detrend
x = detrend(x0,'linear');
I=find(isfinite(x));
x=x(I);
pres=pres(I);

%% Spatial resolution of the time series (cpm)
Fs=length(pres)./ (max(pres)-min(pres));
if isempty(x) | sum(x==0)==length(pres)
    return;
end

%% Degrees of freedom and MAD threshold
NsegM = 2 * floor(length(x(:,1))/sL)-1;
dof = 1.9*NsegM; %according to ODAS 4.4
MADc = sqrt(2/dof); 

%% Compute raw spectrum PSDT
scalar_info.fft_length      = sL;
scalar_info.spec_length     = length(x);
scalar_info.overlap         = length(x)/2;
scalar_info.fs              = 512;
scalar_info.gradient_method = 'high_pass';
scalar_info.f_AA            = 98;
% This function is provided by RSI in the ODAS libraries.
% Alternatively, pwelch can be used obtaining comparable results.
[PSDT,fr] = csd_odas(x-nanmean(x),x-nanmean(x),sL,scalar_info.fs,[],sOV,'linear');
% [PSDT,fr] = pwelch(x-nanmean(x),hann(sL),sOV,sL,scalar_info.fs,'onesided');

%% convert the frequency spectrum into the corresponding wavenumber spectrum (Taylor's frozen turbulence hypothesis)
PSDT=PSDT*W;
k = fr/W;
PSD=PSDT;

%% Time-response correction
if strcmp(time_corr,'KOC') % as in Kocsis et al. (1999)
    tau = tau_0*W.^-0.5;
elseif strcmp(time_corr,'RSI') % as in Vachon and Lueck (1984)
    F0 = 25*sqrt(W);
    tau = (2*pi()*F0/sqrt(sqrt(2)-1))^(-1);
elseif strcmp(time_corr,'NAS') % as in Nash et al. (1999)
    tau = tau_0*W.^-0.12;
elseif strcmp(time_corr,'SOM') % as in Sommer et al. (2013)
    tau = 0.010;
else
    error('Error: time_corr can be RSI, KOC or NAS')
end

%% Single or double correction
if strcmpi(npoles,'single')
    H = 1./(1 + (2*pi()*tau*fr).^2);
    H_lim=0.1;
elseif strcmpi(npoles,'double')
    H = 1./(1 + (2*pi()*tau*fr).^2).^2;
    H_lim=0.1;
end
PSD_raw=PSD;  % Raw spectrum
PSD = PSD./H; % Corrected spectrum

%% Noise function (according to the function provided by RSI in the ODAS libraries)
noise_info=gradT_noise_odas;
noise_info.gamma_RSI = 1;
noise_info.E_n = 3e-9;
Sn = gradT_noise_odas(T_dT, T_string, W, fr, setupstr,noise_info); % Noise spectrum [(K/m)^2/Hz].
Sn = Sn*W; % [K^2/m]
Sn = Sn./H; % Correction for time response
Sn(1) = 0;

% Alternatively, the noise can be computed as follows
% noisep = sensor dependent parameters. E.g. [-10.24,-0.89,fAA]
% where fAA = frequency of the anti-aliasing filter
% Sn = FP07noise(noisep,fr);
% Sn = Sn*W;
% Sn = Sn.*(2*pi()*k).^2;
% Sn = Sn./H;
% Sn(1) = 0;

mpres=mean(pres);
kn = fn/W;
ik1=find(k>k1,1,'first');
ikn = find(k<=kn,1,'last');
ikB=find(k>=kB_S,1,'first'); if isempty(ikB), ikB=length(k); end

%% Delete undesired parts of the spectrum
k0 = k;
PSD0 = PSD;
Sn0 = Sn;
H0 = H;
k = k(1:ikn);
kU = k(end);
PSD = PSD(1:ikn);
H = H(1:ikn);
Sn = Sn(1:ikn);

%% Variance in the noise free part determined from the noise model
% ksfact=0.04; % Upper bound of the inertial-convective subrange, as in Steinbuck et al. (2009)
% Snfact=1.55; % Maximum acceptable signal to noise ratio 1.55, as in Goto et al. (2016)


iknM0 = find(PSD<Snfact*Sn | H<H_lim | k==max(k),1,'first'); % Upper wavenumber limit of the fitting range

kS = ksfact*Pr^(-0.5)*kB_S; % Lower wavenumber limit according to Steinbuck et al. (2009)
k11 = max([kS,k1]);
ik11 = find(k>=k11,1,'first');
ikcor = ik11:iknM0;
Xiv = 6*D*sum( (PSD(ikcor(2:end))+PSD(ikcor(1:end-1))).*(k(ikcor(2:end))-k(ikcor(1:end-1))) )/2;
Xin = 6*D*sum( (Sn(ikcor(2:end))+Sn(ikcor(1:end-1))).*(k(ikcor(2:end))-k(ikcor(1:end-1))) )/2;
cont = false;

%% A signal-to-noise ratio larger than 1.3 was required in the integration range
if Xiv>1.3*Xin
    cont = true;
    Xiv = Xiv -Xin;
else
    Xiv = 0;
    Xi_T = 0;
end
if cont
    
    %% Calculates Xi using kB_S from shear probes and using the theoretical spectrum to correct for unresolved variance
    % See eqs. (6) and (7) in the main text
    if kB_S>0
        BAT=Tspec(Tdis,Xiv,kB_S,k,D,q);
        Xi_ref=6*D*sum( (BAT(ikcor(2:end))+BAT(ikcor(1:end-1))).*(k(ikcor(2:end))-k(ikcor(1:end-1))) )/2;
        Xi_ST=Xiv.*Xiv./Xi_ref;  % See eq. (7) in the main text.
        BAT=Tspec(Tdis,Xi_ST,kB_S,k,D,q);
        %% Evaluation of MAD
        MAD_ST = meanabsdev( PSD(ikcor), BAT(ikcor), Sn(ikcor) );
        fit_flag_ST=0;
        if MAD_ST<2*MADc
            fit_flag_ST=1;
        end
    end
    
    %% Fit the parameters in the noise free region
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    kF = k(iknM0);
    ktest = linspace(max([5*k1,kF/5]),kF*5,40);
    cost = nan(size(ktest));
    Xif0 = nan(size(ktest));
    dk = ktest(2)-ktest(1);
    
    %% First search
    for i = 1:length(ktest)
        kS = k1;   % The Steinbuck et al. (2009) lower limit is not considered during the first round
        kL = max([kS,k1]);
        ikL = find(k>=kL,1,'first');
        ikfit = ikL:iknM0;
        [cost(i), Xif0(i),Steo] =  cost_T_fit(k, PSD, ktest(i), ikfit, Tdis, q, dof, Sn,D,Snfact);
    end
    
    LKHtest = - cost;
    ML = max(LKHtest);
    iML = find(LKHtest == ML);
    kB00 = ktest(iML);
    BATf0 = Tspec(Tdis,Xif0(iML),kB00,k,D,q);
    kP=kB00/sqrt(6*q); % k at the peak
    
    %% Look for the lower integration limit kL
    if strcmp(int_range,'S')     % If only Steinbuck et al. (2009) condition
        kLuk=k1; islope=2;
    elseif strcmp(int_range,'L') % If only Steinbuck et al. (2009) + Luketina and Imberger (2001) conditions
        % Find the intersection between observed and theoretical spectra
        diffS = (PSD-BATf0)./abs(PSD-BATf0);
        ikL = find(diffS<0,1,'first');
        if ~isempty(ikL) && abs(PSD(ikL)-BATf0(ikL))>abs(PSD(ikL-1)-BATf0(ikL-1))
            ikL=ikL-1;
        end
        % This was not mentioned in Luketina, it is a double check on
        % the existence of an actual finestructure pattern.
        % Calculate the slope of the first part of the spectrum using an increasing window
        % Note: k(1) and PSD(1) = 0 --> the window starts from the 2nd point.
        % The iteration starts from the 4th point to have at least 3 points to calculate the slope
        ikP = find(k<=kP,1,'last');
        slope=NaN(1,iknM0);
        for i=4:iknM0
            X=log10(k(2:i)); Y=log10(PSD(2:i));
            X=[ones(length(X),1) (X)];
            b=X\Y;     % \ computes the least square fit
            slope(i)=b(2);
        end
        if ~isempty(ikL)
            [vslope,islope]=nanmin(slope(1:min(ikP,ikL))); % find the minimum of the slope
            if vslope>0 || nanmean(slope(1:islope))>0 % if the minimum slope is > 0, then
                islope=2;
            end
        else
            islope=2;
        end
        kLuk = min([k(islope),k(ikL),kP]);
    end
        
    %% Second search
    cost = nan(size(ktest));
    Xif0 = nan(size(ktest));
    dk = ktest(2)-ktest(1);
    for i = 1:length(ktest)
        kS = ksfact*Pr^(-0.5)*kB00;    % Note: kS is based on the previous best
        kL = max([kS,kLuk]);           % Note, if int_range=='L' here we consider the Luketina and Imberger (2001) condition
        ikL = find(k>=kL,1,'first');
        ikfit = ikL:iknM0;
        [cost(i), Xif0(i)] =  cost_T_fit(k, PSD, ktest(i), ikfit, Tdis, q, dof, Sn,D,Snfact);
    end
    LKHtest = - cost;
    ML = max(LKHtest);
    iML = find(LKHtest == ML);
    kB00 = ktest(iML);
    BATf0 = Tspec(Tdis,Xif0(iML),kB00,k,D,q);
    kP=kB00/sqrt(6*q); % k at the peak
    
    %% Update the lower integration limit
    if strcmp(int_range,'S')     % If only Steinbuck et al. (2009) condition
        kLuk=k1;
    elseif strcmp(int_range,'L') % If only Steinbuck et al. (2009) + Luketina and Imberger (2001) conditions
        diffS = (PSD-BATf0)./abs(PSD-BATf0);
        ikL = find(diffS<0,1,'first');
        ikP = find(k<=kP,1,'last');
        if ~isempty(ikL) && abs(PSD(ikL)-BATf0(ikL))>abs(PSD(ikL-1)-BATf0(ikL-1))
            ikL=ikL-1;
        end
        if ~isempty(ikL)
            [vslope,islope]=nanmin(slope(1:min(ikP,ikL))); % find the minimum of the slope
            if vslope>0 || nanmean(slope(1:islope))>0 % if the minimum slope is > 0, then
                islope=2;
            end
        else
            islope=2;
        end
        kLuk = min([k(islope),k(ikL),kP]);
    end
    ikL = find(k>=kLuk,1,'first');
    
    %% Third, refined search
    MLp = -cost_T_fit(k, PSD, kB00+dk, ikfit, Tdis, q, dof, Sn,D,Snfact);
    MLm = -cost_T_fit(k, PSD, kB00-dk, ikfit, Tdis, q, dof, Sn,D,Snfact);
    deltak0=abs((2*dk)/sqrt(2*ML-MLm-MLp));
    deltak = max([deltak0,dk]);
    
    kmin2 = max([kB00-deltak,k(ikL+1)]);
    kmax2 = kB00+deltak;
    ktest = linspace(kmin2,kmax2,40);
    clear cost Xif0
    cost = nan(size(ktest));
    Xif0 = nan(size(ktest));
    for i = 1:length(ktest)
        kS = ksfact*Pr^(-0.5)*kB00;    % Note: kS is based on the previous best
        kL = max([kS,kLuk]);            % Note, if int_range=='L' here we consider the Luketina and Imberger (2001) condition
        ikL = find(k>=kL,1,'first');
        ikfit = ikL:iknM0;
        [cost(i), Xif0(i)] =  cost_T_fit(k, PSD, ktest(i), ikfit, Tdis, q, dof, Sn,D,Snfact);
    end
    
    LKHtest = - cost;
    MLKH = max(LKHtest);
    iML = find(LKHtest == MLKH);
    Xi_T = Xif0(iML);
    kB = ktest(iML);
    kP=kB/sqrt(6*q); % k peak
    BATf = Tspec(Tdis,Xi_T,kB,k,D,q);
    
    %% Compute quality metrics
    MAD_T = meanabsdev( PSD(ikfit), BATf(ikfit), Sn(ikfit) );
    % Re-compute the integration range (for figure and fit_flag_T)
    kS = ksfact*Pr^(-0.5)*kB;
    kL = max([kS,kLuk]);
    ikL = find(k>=kL,1,'first');
    kL = k(ikL);
    ikfit = ikL:iknM0;    
    % Calculates uncertainties in the fitting parameters
    [MLp,Xip,~] = cost_T_fit(k, PSD, kB+deltak, ikfit, Tdis, q, dof, Sn,D,Snfact);
    [MLm,Xim,~] = cost_T_fit(k, PSD, kB-deltak, ikfit, Tdis, q, dof, Sn,D,Snfact);
    MLp = - MLp;
    MLm = - MLm;
    sKBT=abs((2*deltak)/sqrt(2*MLKH-MLm-MLp));
    sXif = abs((Xip-Xim)/sqrt(2*MLKH-MLm-MLp));    
    % Fits to polynom (avoiding the noisy part)
    ikfitA = ikfit(1):min([ikfit(end),find(BATf<Snfact*Sn,1,'first')]);
    logK = log(k(ikfitA));
    logS = log(PSD(ikfitA));
    pp=polyfit(logK, logS,1);
    Sm = exp(polyval(pp, log(k)));
    LKHpol = -cost_MLE(PSD(ikfit), Sm(ikfit), dof, Sn(ikfit));    
    LR = MLKH - LKHpol;
    LR = log10(exp(1))*LR;   % Likelihood ratio, see eq. (16)
    kU=k(iknM0);
    krange = log10(k(iknM0)) - log10(k(ikL));
    fit_flag_T =0;
    if LR>2 && MAD_T<MADc*2 ...
            && krange> 0.8 ...
            && kU>2*kP && kL<kP
        fit_flag_T =1;
    end
    eps_T=visco*D^2*(2*pi()*kB)^4; % TKE dissipation rate

    %% Export spectra information (T. Doda, 05.02.2025)
    SPECTRUM.k0=k0; % All wavenumbers
    SPECTRUM.k=k; % Wavenumbers up to the maximum wavenumber specified from the anti-aliasing filter f_AA (fn)
    SPECTRUM.PSD_raw=PSD_raw; % Initial spectrum (all wavenumbers)
    SPECTRUM.PSD0=PSD0; % Spectrum corrected for time response (all wavenumbers)
    SPECTRUM.PSD=PSD; % Spectrum corrected for time response over k
    SPECTRUM.ind_fit=ikfit; % Indices of wavenumbers used for fitting, from max([kS, kuK]) to upper wavenumber defined from signal-to-noise ratio
    SPECTRUM.PSD_theo=Tspec(Tdis, Xi_T, kB,k0,D,q); % Theoretica spectrum (all wavenumbers)
    SPECTRUM.Sm=Sm; % Exponential fitted PSD over k
    SPECTRUM.Sn0=Sn0; % Noise spectrum corrected for time response (all wavenumbers)
    SPECTRUM.Sn=Sn; % Noise spectrum corrected for time response over k
    SPECTRUM.H=H; % Time response correction parameters over k
    SPECTRUM.H_lim=H_lim; % Maximal value of H

    %% Plots
    if plt~=0         
        fh=figure;
        clf
        set(fh,'color','white','Units', 'Centimeters', 'Position', [0,0,16,8],...
            'PaperUnits', 'Centimeters', 'PaperSize', [16,8]);
        ax1=axes('position',[0.08,0.175,0.2,0.72]);
        plot(x0-mean(x0),pres,'k')
        axis ij
        grid('on')
        xlabel('\partial T^{\prime}/\partial z [^\circ C m^{-1}]','fontsize',10)
        ylabel('Pressure [dbar]','fontsize',10)
        ax2=axes('Position',ax1.Position,'XAxisLocation','top',...
            'YAxisLocation','right','color','none',...
            'xColor',[0,114,178]/255,'yColor','k');
        set(ax1,'box','off');  set(ax2,'YtickLabels',[]);
        line(movavg(Tplot,'linear',100),presplot,'color',[0,114,178]/255); set(gca,'ydir','reverse');
        line(Tplot,presplot,'color',[0,114,178]/255); set(gca,'ydir','reverse');
        xlabel('Temperature [°C]','fontsize',10);
        
        ax=axes('position',[0.40,0.175,0.55,0.72]);
        s_raw=loglog(k0,PSD_raw, '-k','linewidth',0.5); hold on
        s_corr=loglog(k0,PSD0,'-','color',[75 97 209]/255,'linewidth',0.75);
        if cont
            loglog(k(ikfit),PSD(ikfit),'ok','markerfacecolor',[75 97 209]/255, 'markersize',5)
            s_fit=loglog(k0,Tspec(Tdis, Xi_T, kB,k0,D,q), 'color',[75 97 209]/255, 'linewidth',2);
            powerfit=loglog(k(ikfit),Sm(ikfit)+Sn(ikfit),'-','color',[174  38  43]/255,'linewidth',2);
        end
        s_noise=loglog(k0,Sn0, ':k', 'linewidth',0.7);
        ylim([10^-9,1e3])

        text(1.2,1e-4,['\chi_{T} = ', num2str(Xi_T,'%1.3e'), '^{\circ}C^2 s^{-1}'],'horizontalalignment','left','fontsize',8)
        text(1.2,10^-4.9,['k_B = ', num2str(kB,'%3.0f'),' cpm'],'horizontalalignment','left','fontsize',8)
        text(1.2,10^-5.7,['\epsilon_T = ', num2str(eps_T,'%1.2e'),' m^2 s^{-3}'],'horizontalalignment','left','fontsize',8)
        text(1.2,10^-6.6,['LR = ', num2str(LR,'%1.1f'),'(2)'],'horizontalalignment','left','fontsize',8)
        text(1.2,10^-7.45,['MAD = ', num2str(MAD_T,'%1.2f'),' (', num2str(2*MADc, '%1.2f'),')'],'horizontalalignment','left','fontsize',8)
        
        plot(kB,1.4E-9,'^','markeredgecolor','k','markerfacecolor','k','markersize',5)
        plot(fn/W,1.4E-9,'^','markeredgecolor',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5],'markersize',5)
        itmp=find(PSD0<Snfact*Sn0 ,1,'first');
        if isempty(itmp)
            itmp=length(k0);
        end
        plot(k0(itmp),1.4E-9,'^','markerfacecolor',[0.75 0.75 0.75],'markeredgecolor',[0.75 0.75 0.75],'markersize',5)
        plot(k0(find(H<H_lim ,1,'first')),1.4E-9,'^','markeredgecolor','k','markerfacecolor',[1 1 1],'markersize',5)
        plot(kP,1.4E-9,'^','markeredgecolor',[75 97 209]/255,'markerfacecolor',[75 97 209]/255,'markersize',5)
        xlabel('k [cpm]','fontsize',10)
        ylabel('\Psi_T [(^{\circ}C m^{-1})^{-2} cpm^{-1}]','fontsize',10)
        xlim([1,1000])
        set(gca,'ytick',[1e-9 1e-7 1e-5 1e-3 1e-1 1e1 1e3]);
        legend([s_raw,s_corr,s_fit,powerfit,s_noise],{'Raw spec.','Corrected spec.','Kraichnan spec. (fit)','Power law (fit)','Noise spec.'},...
            'location','North','orientation','horizontal','NumColumns',2)
        grid on; ax.GridAlpha=0.1;
        set(gca, 'xminorgrid', 'on')
        set(gca,'MinorGridLineStyle','-','MinorGridAlpha',0.1)
        set(gca,'GridLineStyle','-')
%         saveas(gcf,[folder_out '/Spec_z= ', num2str(mpres), ' m, sensor= ', T_string, '.pdf'])
%         saveas(gcf,[folder_out '/Spec_z= ', num2str(mpres), ' m, sensor= ', T_string, '.png'])  
    end
end
