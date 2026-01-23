% Compare turbulence  profiles

close all
clear
clc

%% Parameters

savefig=false;
% campaign_name='20211110';
campaign_name='20250211';
datafolder=['..\..\data\VMP\',campaign_name,'\Level2\'];
folder_param_name='down1.0_K5.3_L_NAS0.0058_singlepole_nfft3_EPFL_Goodman';
filenames={'VMP004','VMP007','VMP010','VMP014'};

addpath("..\..\functions\figures\") % Add functions

%% Load the data
disp('Loading the data...')


for kf=1:length(filenames)
    datavar=load([datafolder,filenames{kf},'_',folder_param_name,'\results_',filenames{kf},'_down.mat']);
    datavar.filename=filenames{kf};
    DATA(kf)=datavar;

end


%% Figures turbulence
col_day1=winter(length(filenames));
%col_day1(3,:)=[0,0,0];
[fig1,fig2,fig3,fig4]=plot_comparison_turbulence(DATA,filenames,'longest',col_day1);

if savefig
    saveas(fig1,['CTD_profiles_',campaign_name,'_',filenames{:},'.fig'])
    exportgraphics(fig1,['CTD_profiles_',campaign_name,'_',filenames{:},'.png'],'Resolution',400)

    saveas(fig2,['fast_profiles_',campaign_name,'_',filenames{:},'.fig'])
    exportgraphics(fig2,['fast_profiles_',campaign_name,'_',filenames{:},'.png'],'Resolution',400)

    saveas(fig3,['turbulence_profiles_',campaign_name,'_',filenames{:},'.fig'])
    exportgraphics(fig3,['turbulence_profiles_',campaign_name,'_',filenames{:},'.png'],'Resolution',400)
end

