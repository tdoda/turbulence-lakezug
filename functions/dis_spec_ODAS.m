function [epsilon, W, fit_flag, MAD, MADc] = dis_spec_ODAS(sh, A, fsf, fss, W, Tf,timef, Pf, visco,Nfft,  plt,folder_out, noise_corr)
%info for epsilon ODAS calculation
        % - set up info struct
       
        infoD.fft_length    = Nfft;
        infoD.diss_length   = length(sh);
        infoD.overlap       = length(sh);
        infoD.fs_fast       = fsf;
        infoD.fs_slow       = fss;
        infoD.speed         = abs(W);
        infoD.T             = Tf;
        infoD.t             = timef;
        infoD.P             = Pf;
        infoD.viscosity     = visco;            % Added for being used in get_diss_odas_EPFL
        infoD.noise_corr    = noise_corr;       % Added to choose the noise algorithm

        diss = get_diss_odas_EPFL(sh, A, infoD);
        MAD = diss.mad(1);
        dof = diss.dof_spec;
        MADc = sqrt(2/dof);
        fit_flag = 0;
        if MAD<2*MADc
            fit_flag = 1;
        end
        
        epsilon = diss.e;
        mpres = diss.P;
        W = diss.speed;

        Kn = diss.K_max;
        K = diss.K;        
        iKn = find(K<=Kn,1,'last');
        Kc=1/(2*pi())*(epsilon./visco^3).^(1/4);
        iKc = find(K<=Kc,1,'last');
        PSD = diss.sh_clean;
        PSDnc = diss.sh;
        NAS = diss.Nasmyth_spec;
        
        if plt(1) ~=0
            
            figure(11)
            clf
            set(gcf, 'PaperUnits', 'centimeters');
            set(gcf, 'PaperSize', [18 8]); 
            set(gcf, 'PaperPositionMode', 'manual') ;
            set(gcf, 'PaperPosition', [0 0 18 8]);
            axes('position',[0.1,0.12,0.2,0.8])
            plot(sh,Pf,'k')
            axis ij
            grid('on')
            xlabel('shear (s^{-1})')
            ylabel('p (db)')
            
            axes('position',[0.4,0.12,0.5,0.8])
            loglog(logspace(0,3,100),nasmyth(logspace(-11,-3,9),visco,logspace(0,3,100)),'color',[0.60 0.60 0.60])
            hold on
            for i=1:9
               text(1.5,min(nasmyth(10^(-12+i),visco,[1.5 2])),['10^{',num2str(-12+i),'} W/kg'],'color',[0.60 0.60 0.60])
            end
            loglog(K,PSDnc,'color',[0.5,0.5,0.5])
            loglog(K,PSD,'k')    
            loglog(K(1:iKn),PSD(1:iKn),'ko', 'markersize',4)
            loglog(K,NAS,'k','linewidth',2)
            ylim([10^-9,2])
            xlim([1 10^3])
            line([Kn,Kn],ylim,'color','k','linestyle','--')
            line([Kc,Kc],ylim,'color','k')
            text(1.1*Kc,0.1,'K_c')
            text(1.2,10^-0.5,['\epsilon = ', num2str(epsilon,'%1.3e'),' W/kg'],'horizontalalignment','left')

            text(2, 10^-8, ['MAD = ', num2str(MAD, '%1.2f'),'(', num2str(MADc*2, '%1.2f'),')'],'horizontalalignment','left')
            title(['ODAS epsilon calculation at ',num2str(mean(Pf)), 'db'])
            xlabel('k (cpm)')
            ylabel('S_{sh} (s^{-2} cpm^{-1})')
            saveas(gcf,[folder_out '/Spec_z= ', num2str(mpres), ' m, SH sensor.png'])

%             pause()
        
        end