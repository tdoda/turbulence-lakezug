% Compare turbulence  profiles

close all
clear
clc

%% Parameters

savefig=false;
campaign_name='20241204';
datafolder=['..\..\data\',campaign_name,'\Level2\'];
folder_param_name='down1.0_K5.3_L_NAS0.0058_singlepole_nfft3_EPFL_Goodman';
filenames={'VMP004','VMP006','VMP003'};

addpath("..\..\functions\figures\") % Add functions

%% Load the data
disp('Loading the data...')


for kf=1:length(filenames)
    datavar=load([datafolder,filenames{kf},'_',folder_param_name,'\results_',filenames{kf},'_down.mat']);
    datavar.filename=filenames{kf};
    DATA(kf)=datavar;
end


%% Figures turbulence
col_day1=winter(3);
col_day1(3,:)=[0,0,0];
[fig1,fig2]=plot_comparison_turbulence(DATA,filenames,'all',col_day1,{'PC1','PC2','PSB'});

if savefig
    saveas(fig1,['turbulence_profiles',filenames{:},'.fig'])
    exportgraphics(fig1,['turbulence_profiles',filenames{:},'.png'],'Resolution',400)
end

