% Check turbulence spectra in a given bin

close all
clear
clc

%% Parameters

savefig=false;
campaign_name='20241205';
datafolder=['..\..\data\VMP\',campaign_name,'\Level2\'];
folder_param_name='down1.0_K5.3_L_NAS0.0058_singlepole_nfft3_EPFL_Goodman';
profname='VMP001';
depth_bin=[50]; % [m]

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
   
    fig1=plot_spectrum_FP07(data_prof_bin,DATA.FAST{indprof},ind_bin);
end

% if savefig
%     saveas(fig1,['turbulence_profiles',filenames{:},'.fig'])
%     exportgraphics(fig1,['turbulence_profiles',filenames{:},'.png'],'Resolution',400)
% end

