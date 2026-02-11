function [eps_S, MAD_S, MADc, fit_flag_S, kL, kU, SPECTRUM]=TKE_dis_spec(pres,x0,px0,k1,k2,fn,visco,W,sL,sOV,noise_corr,sh_number,plt,presplt,Tplt)

% GOAL
% Estimate of TKE dissipation rate eps_S by integration of the measured
% velocity shear wavenumber spectrum

% INPUT:
% pres: pressure vector (db)
% x0: shear matrix [S1 S2] (1/s) (both are required for Goodman correction)
% px0: vibrations matrix [Ax Ay] (counts) (both are required for Goodman correction) 
% k1: minimium wavenumber for integration (small number, e.g., 0.1 cpm)
% k2: smaller wavenumber for which the observed spectrum should be defined to proceed with the eps_S estimate (e.g., 14 cpm) 
% fn: maximum frequency for calculations (90% of anti-aliasing filter f_AA. E.g., if f_AA=98 Hz, fn=0.9*98=88.2 Hz)
% visco: kinematic viscosity of water (m^2/s)
% W: mean profiling speed (m/s)
% sL: length of segments for fft (scans)
% sOV: overlap for fft (scans)
% noise_corr: type of denoising correction 'none' or 'Goodman'
% sh_number: shear probe for which eps_S is calculated 'sh1' or 'sh2'
% plt: flag for making the figure (0 no figure, ~=0 make figure)
% presplt: pressure vector for making the plot
% Tplt: temperature vector for making the plot

% OUTPUT
% eps_S: TKE dissipation rate (m^2/s^3)
% MAD_S: Mean Absolute Deviation between observed and empirical spectra
% MADc: Threshold for the Mean Absolute Deviation between observed and empirical spectra
% fit_flag_S: acceptance flag according to quality metrics: 0 rejected, 1 accepted
% kL: Lower integration wavenumber (cpm)
% kU: Upper integration wavenumber (cpm)
% SPECTRUM: structure array containing all information to plot the spectra
% [added by T. Doda, 10.02.2026]


% Last version: 10.02.2026

%% Initialization
eps_S=NaN;kK=NaN;epsN=NaN;kKN=NaN;

[~,ncsh]=size(x0);
[~,ncA]=size(px0);
if (ncA~=ncsh)
    disp('Error: number of columns of AA and sh is not the same')
    pause
else
    nc=ncsh; % Number of shear sensors
end

%% Delete NaN and detrend
I = isfinite(x0); I=find(I(:,1)==1 & I(:,nc)==1); % generalized for the case with 2 sh
x0 = x0(I,:);
px0 = px0(I,:);
pres = pres(I);
x = detrend(x0,'linear');
px = detrend(px0,'linear');

I = isfinite(x); I=find(I(:,1)==1 & I(:,nc)==1);
x = x(I,:);
px = px(I,:);
pres = pres(I);

if isempty(x) |  sum(x==0)==length(pres)
    return;
end

%% Spatial resolution of the time series (cpm), FFT length and overlapping
Fs=length(pres)./ (max(pres)-min(pres)); % Mean sampling rate in pressure space [number of samples/dbar] 
% FFT segment length (number of samples) if not specified (half the bin by
% default):
if sL == 0
    sL = round(length(x(:,1))/2)-1; 
end
% Number of samples shared between overlapping segments if not specified
% (50 % overlap by default)
if sOV == 0
    sOV = round(sL/2); 
end

%% Degrees of freedom and MAD threshold
NsegM = 2 * floor(length(x(:,1))/sL)-1;
dof = 1.9*NsegM; %according to ODAS 4.4
MADc = sqrt(2/dof);

%% Compute raw spectrum PSD0
if nc==1 % Only one shear sensor
    % This function is provided by RSI in the ODAS libraries.
    % Alternatively, pwelch can be used obtaining comparable results.
    [PSD0,k] = csd_odas(x,x,sL,Fs,[],sOV,'linear');
    %     [PSD0,k] = pwelch(x-nanmean(x),hann(sL),sOV,sL,Fs,'onesided');
elseif nc==2 % Two shear sensors
    if strcmp(sh_number,'sh1')
        column_index=1;
    elseif strcmp(sh_number,'sh2')
        column_index=2;
    end
    % This function is provided by RSI in the ODAS libraries.
    % Alternatively, pwelch can be used obtaining comparable results.
    [PSD0,k,~,~] = csd_matrix_odas(x(:,column_index),x(:,column_index),sL,Fs,[],sOV,'linear');
    %     [PSD0, k] = cpsd(x(:,column_index)-nanmean(x(:,column_index)),x(:,column_index)-nanmean(x(:,column_index)),hann(sL),sOV,sL,Fs,'onesided');
end

%% Estimate epsilon only when there are points for k<k2
% (may not be the case if fft_length is short, e.g., at the beginning
% or end of the profile.
if k(2)<k2
    %% Compute denoised spectrum PSD
    if strcmpi(noise_corr,'Goodman')
        if nc==1
            if sum(abs(px))>0
                % Goodman denoising if only 1 sh is available.
                % These functions are provided by RSI in the ODAS libraries.
                % Alternatively, pwelch can be used obtaining comparable results.

                % Auto-spectrum of vibration signal:
                [PSDps, ~] = csd_odas(px,px,sL,Fs,[],sOV,'linear');

                % Cross-spectrum of vibration and shear;
                [CPS, ~, ~, ~] = csd_matrix_odas(x,px,sL,Fs,[],sOV,'linear');
                %                 [PSDps, ~] = pwelch(px-nanmean(px),hann(sL),sOV,sL,Fs,'onesided');
                %                 [CPS, ~] = cpsd(x-nanmean(x),px-nanmean(px),hann(sL),sOV,sL,Fs,'onesided');
                H = CPS./PSDps; % Spectral density ratio
                PSDcont = abs(H).^2.*PSDps;
                PSD = PSD0 - PSDcont;
            else
                PSD = PSD0;
            end
        elseif nc==2
            % Goodman denoising if 2 sh are available.
            % This function is provided by RSI in the ODAS libraries.
            [PSD, ~,~,~,~] = clean_shear_spec(px, x, sL, Fs);
            PSD = permute(PSD, [3 1 2]);
            PSD = squeeze(PSD(:, column_index, column_index));
        end
    elseif strcmpi(noise_corr,'none')
        PSD = PSD0;
    end
    
    if nc==2
        x=x(:,column_index);
    end
    
    %% Apply spatial response correction
    Hsr = 1./(1+(k./48).^2);
    PSD = PSD./Hsr;
    PSD0 = PSD0./Hsr;
    
    dk = k(2) - k(1);
    mpres=mean(pres);
    
    ik1=find(k>=k1,1,'first');
    ik2=find(k>=k2,1,'first');
    if ik2 == ik1
        ik2 = ik2+2;
    end
    
    %% Adjust kn = upper integration wavenumber
    kn=fn/W; kn=min(150,kn);       % it is not wise to extend the correction beyond 150 cpm because, at that wavenumber the spectrum is boosted by a factor of >10.
    ikn=find(k>=kn,1,'first');
    
    %% EPSILON ITERATIVE CALCULATION
    ik3=ik2;
    k3=k(ik3);
    kK0=k2;
    flag=0;
    while flag==0
        eps_S=7.5*visco*sum( (PSD(ik1+1:ik3)+PSD(ik1:ik3-1)).* (k(ik1+1:ik3)-k(ik1:ik3-1)) )/2;
        %%%%
        kK=0.121*(eps_S./visco^3).^(1/4);  % 0.121=0.76/(2*pi) (0.76 of kK). kK=Kolmogorov wavenumber
        %%%%
        if abs(kK-kK0)<=2*dk
            flag = 1;
            %             ik3 = find(k>=kK,1,'first');
            ik3 = find(k<=kK,1,'last');  % More cautelative to avoid the noisy part
        elseif kK<=kn
            inc = (kK-k3)/abs(kK-k3);
            ik3= ik3+inc; %find(k >= kK,1,'first');
            k3=k(ik3);
            kK0 = kK;
        else
            ik3=ikn;
            flag = 2;
        end
        
    end
    kU=k(ik3);   % Upper integration limit
    kL=k(ik1);   % Lower integration limit
    if ik3<ik1+1
        disp('Shear epsilon (%s) cannot be calculated because of too low Kolmogorov scale',sh_number)
        % Leave it returning an error in that case
        % eps_S=NaN;
        % MAD_S=NaN;
        % MADc=NaN;
        % fit_flag_S=NaN;
        return
    end
    eps_S=7.5*visco*sum( (PSD(ik1+1:ik3)+PSD(ik1:ik3-1)).* (k(ik1+1:ik3)-k(ik1:ik3-1)) )/2;
    kK=1/(2*pi())*(eps_S./visco^3).^(1/4); % Kolmogorov wavenumber
    
    %% Correct for lost variance, see eqs. (6) and (7) in the main text
    % This function is provided by RSI. It contains the definition of the
    % Nasmyth empirical spectrum. See eqs.(11) and (12) in the main text.
    NAS=nasmyth(eps_S,visco,k);
    varianceN=sum( (NAS(ik1+1:ik3)+NAS(ik1:ik3-1)).* (k(ik1+1:ik3)-k(ik1:ik3-1)) )/2;
    eps_ref=7.5*visco*varianceN;
    eps_S=eps_S*eps_S/eps_ref; % See eq. (7) in the main text.
    kK=1/(2*pi())*(eps_S./visco^3).^(1/4);
    
    %% EPSILON CALCULATION BY FITTING TO NASMYTH 
    % (not available in the database)
    epsN = nan;
    kKN = nan;
    try
        epsN=nlinfit(k(ik1:ik3),log(PSD(ik1:ik3)),@(e,kk)log(nasmyth(e,visco,kk)),10^-9);
        kKN=1/(2*pi())*(eps_S./visco^3).^(1/4); % Kolmogorov wavenumber
    end
    
    NAS2=nasmyth(epsN,visco,k);
    NAS=nasmyth(eps_S,visco,k);
    
    %% Evaluation of MAD
    MAD_S = mean( abs( (PSD(ik1:ik3)./NAS(ik1:ik3)) - mean(PSD(ik1:ik3)./NAS(ik1:ik3)) )  );
    MADf = mean( abs( (PSD(ik1:ikn)./NAS(ik1:ikn)) - mean(PSD(ik1:ikn)./NAS(ik1:ikn)) )  );
    fit_flag_S = 0;
    if MAD_S<2*MADc
        fit_flag_S = 1;
    end

    %% Export spectra information (T. Doda, 10.02.2026)
    SPECTRUM.k      = k;        % Wavenumbers used for the shear spectrum
    SPECTRUM.PSD0   = PSD0;     % Raw / corrected-for-response shear spectrum
    SPECTRUM.PSD    = PSD;      % Denoised shear spectrum
    SPECTRUM.ind_fit=1:ik3; % Indices of wavenumbers used for fitting
    SPECTRUM.PSD_theo    = NAS;      % Nasmyth spectrum (iteration result)
    SPECTRUM.kK     = kK;       % Kolmogorov wavenumber
    SPECTRUM.kn   = fn/W;     % 90% of the anti-aliasing cutoff wavenumber
    SPECTRUM.k_lim  = 150;      % Hard upper k-limit used in plot
    
    %% Plots
    if plt~=0        
        fh=figure;
        clf
        set(fh,'color','white','Units', 'Centimeters', 'Position', [0,0,16,8],...
            'PaperUnits', 'Centimeters', 'PaperSize', [16,8]);
        ax1=axes('position',[0.08,0.175,0.2,0.72]);
        plot(x,pres,'k')
        hold on
        axis ij
        grid('on')
        xlabel('Shear [s^{-1}]','fontsize',10)
        ylabel('Pressure [dbar]','fontsize',10)
        ax2=axes('Position',ax1.Position,'XAxisLocation','top',...
            'YAxisLocation','right','color','none',...
            'xColor',[0,114,178]/255,'yColor','k');
        set(ax1,'box','off'); set(ax2,'YtickLabels',[]);
        line(movavg(Tplt,'linear',100),presplt,'color',[0,114,178]/255); set(gca,'ydir','reverse');
        line(Tplt,presplt,'color',[0,114,178]/255); set(gca,'ydir','reverse');
        xlabel('Temperature [°C]','fontsize',10);
                
        ax=axes('position',[0.40,0.175,0.55,0.72]);
        loglog(logspace(0,3,100),nasmyth(logspace(-11,-3,9),visco,logspace(0,3,100)),'color',[0.60 0.60 0.60])
        hold on
        text(1.1,0.5,'\epsilon_{th} [m^2s^{-3}]','color',[0.60 0.60 0.60],'fontsize',8,'horizontalalignment','left')
        for i=1:9
            text(1.1,min(nasmyth(10^(-12+i),visco,[1.5 2])),['10^{',num2str(-12+i),'}'],'color',[0.60 0.60 0.60],'fontsize',7)
        end
        s_raw=loglog(k,PSD0,'-k','linewidth',0.5);
        s_corr= loglog(k,PSD,'-','color',[75 97 209]/255,'linewidth',0.75);
        loglog(k(1:ik3),PSD(1:ik3),'ok','markerfacecolor',[75 97 209]/255, 'markersize',5)
        s_fit = loglog(k,NAS,'-','color',[75 97 209]/255,'linewidth',2);
        ylim([10^-9,2])
        xlim([1 10^3])
        plot(kK,1.27E-9,'^','markeredgecolor','k','markerfacecolor','k','markersize',5)
        plot(fn/W,1.27E-9,'^','markeredgecolor',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5],'markersize',5)
        plot(150,1.27E-9,'^','markeredgecolor','k','markerfacecolor',[1 1 1],'markersize',5)
 
        text(1.2,10^-7.2,['\epsilon_{S} = ', num2str(eps_S,'%1.2e'),' m^2 s^{-3}'],'horizontalalignment','left','fontsize',8)
        text(1.2, 10^-7.8, ['MAD_S = ', num2str(MAD_S, '%1.2f'),' (', num2str(2*MADc, '%1.2f'),')'],'horizontalalignment','left','fontsize',8)
        xlabel('k [cpm]','fontsize',10)
        yl=ylabel('\Psi_S [s^{-2} cpm^{-1}]','fontsize',10); yl.Position(1) = 0.55;
        set(gca,'ytick',[1e-10 1e-8 1e-6 1e-4 1e-2 1e0]);
        legend([s_raw,s_corr,s_fit],{'Corrected spec.','Denoised spec.','Nasmyth spec. (iteration)'})
        grid on; ax.GridAlpha=0.1;
        set(gca, 'xminorgrid', 'on')
        set(gca,'MinorGridLineStyle','-','MinorGridAlpha',0.1)
        set(gca,'GridLineStyle','-')
        
%         saveas(gcf,['/Spec_z= ', num2str(mpres), ' db, sensor = ' sh_number '.pdf'])
%         saveas(gcf,['/Spec_z= ', num2str(mpres), ' db, sensor = ' sh_number '.png'])
    end
else
    eps_S=NaN; W=NaN; fit_flag_S=NaN; MAD_S=NaN; MADc=NaN;
end
