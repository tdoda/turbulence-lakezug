% Perform bench test

close all
clear 
clc

%% Parameters
instrument='microCTD'; % Options: 'VMP' or 'microCTD'

if strcmp(instrument,'VMP')
    %profiler_folder='..\Profiler VMP';
    SNname='028';
    odas_folder='..\Profiler VMP\2018_05_15_USB_Stick\odas_v4.3.03\odas';
elseif strcmp(instrument,'microCTD')
    SNname='310';
    %profiler_folder='..\Profiler MicroCTD';
    odas_folder='..\Profiler MicroCTD\Scripts\Sebastiano_private_scripts\odas_v4.4\odas';
else
    error('Wrong instrument name')
end

addpath(odas_folder)

%foldername=[profiler_folder,'\BenchTest\',date_bench_test];

%% Open file selection dialog
[filename, foldername] = uigetfile('*.P', 'Select a bench test data file','..\');
                              
% Check if user selected a file (not canceled)
if filename == 0
    error('No file was selected.');
end

%% Bench test
scriptpath=pwd;
cd(foldername)
quick_bench(filename, SNname,false,true)
data=odas_p2mat(filename);% Default data extraction
cd(scriptpath)


