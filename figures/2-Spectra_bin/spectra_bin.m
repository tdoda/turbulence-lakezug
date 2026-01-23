% Check turbulence spectra in a given bin

close all
clear
clc

%% Parameters

savefig=true;
campaign_name='20250211';
datafolder=['..\..\data\VMP\',campaign_name,'\Level2\'];
% folder_param_name='down2.0_K5.3_L_NAS0.0058_singlepole_nfft3_EPFL_Goodman';
% folder_param_name='down2.0_K5.3_L_NAS0.0058_singlepole_nfft5_EPFL_Goodman';
folder_param_name='down1.0_K5.3_L_NAS0.0058_singlepole_nfft3_EPFL_Goodman';
% folder_param_name='down0.5_K5.3_L_NAS0.0058_singlepole_nfft5_EPFL_Goodman';
profname='VMP003';
depth_bin=[119]; % [m]

addpath("..\..\functions\figures\") % Add functions

%% Load the data
disp('Loading the data...')

DATA=load([datafolder,profname,'_',folder_param_name,'\results_',profname,'_down.mat']);
DATA.filename=profname;

%% Figure spectrum temperature

for kbin=1:length(depth_bin)
    % Get the longest profile 
    depth_range=NaN(1,length(DATA.BINNED)); % [m]
    for kprof=1:length(DATA.BINNED)
        depth_range(kprof)=max(DATA.BINNED{kprof}.depth)-min(DATA.BINNED{kprof}.depth);
    end
    [~,indprof]=max(depth_range);
    data_prof_bin=DATA.BINNED{indprof};

    [~,ind_bin]=min(abs(data_prof_bin.depth-depth_bin(kbin)));
   
    fig1=plot_spectrum_FP07(data_prof_bin,DATA.BINNED{indprof}.SPECTRUMT1,DATA.param,ind_bin,'T1');
    fig2=plot_spectrum_FP07(data_prof_bin,DATA.BINNED{indprof}.SPECTRUMT2,DATA.param,ind_bin,'T2');

    if savefig
        saveas(fig1,['spectrum_',campaign_name,'_',folder_param_name,'_',profname,'_',num2str(ind_bin),'_T1.fig'])
        exportgraphics(fig1,['spectrum_',campaign_name,'_',folder_param_name,'_',profname,'_',num2str(ind_bin),'_T1.png'],'Resolution',400)

        saveas(fig2,['spectrum_',campaign_name,'_',folder_param_name,'_',profname,'_',num2str(ind_bin),'_T2.fig'])
        exportgraphics(fig2,['spectrum_',campaign_name,'_',folder_param_name,'_',profname,'_',num2str(ind_bin),'_T2.png'],'Resolution',400)
        disp('Figures saved!')
    end
end



