%% odas_version_info 
% Return the version of the current ODAS library.
%%
% <latex>\index{Functions!odas\_version\_info}</latex>
%
%%% Syntax
%   version = odas_version_info( )
%
% * []
% * [version] Numeric value of the version.
%
%%% Description
% Function called by other functions to determin the current version of
% the ODAS Libarary.  Used when checking for outdated data files.

% Version History
%
% * 2015-06-03 (WID) initial version 4.0
% * 2016-08-02 (WID) updated to version 4.1.  Added support for some JAC
%                    sensors.  Minor bug fixes.  Notable inclusion of
%                    show_ch fuction.
% * 2017-01-27 (WID) updated to version 4.2.  
%                    addition of AEM1-G analog and digital types
% * 2018-01-05 (JMM) updated to version 4.3.
%                    Bug fixes, including despiking of piezo-accelerometers
%                    This version will be the last major update to the ODAS
%                    library. Minor improvements and bug fixes will be
%                    made, but major changes will be released as part of
%                    the Zissou software package.
% * 2018-04-14 (JMM) updated to version 4.3.02
%                    Fixed some small plotting bugs and annoyances (they
%                    arose during training at Lake Tahoe).
% * 2018-04-25 (JMM) updated to version 4.3.03
%                    Fixed some small plotting bugs and annoyances (they
%                    arose during training at UNAL). Most significant
%                    change is that a plot of pressure will be generated
%                    even if no profiles are detected. 
% * 2018-05-18 (JMM) updated to version 4.3.04
%                    Fixed call to salinity_JAC to include sampling rate
%                    and mean speed. Fixed comments so manual could be
%                    updated. 
% * 2018-08-28 (JMM) updated to version 4.3.05
%                    Added a warning message in read_odas if there is no
%                    channel section in the setup file for an id number 
%                    that is in the address matrix.
% * 2018-09-27 (JMM) updated to version 4.3.06 
%                    Small fix for accelerometer plot for seagliders
% * 2018-09-27 (JMM) updated to version 4.3.07 
%                    Removed two lines that overwrote default speed and 
%                    cutoff frequency in salinity_JAC.m
% * 2019-04-19 (JMM) updated to version 4.3.08 
%                    Added support for type 'raw'. Fixed conflicting
%                    variable name in convert_odas.
% * 2019-06-03 (JMM) Updated to version 4.4
%                    - Added capability of handling horizontal profilers.
%                    - DOF equation for spectra updated to Nuttall formulation
%                    - Figure of Merit estimates now outputted from
%                       get_diss_odas.
%                    - Significant updates to cal_FP07_in_situ (better
%                    control of plots, warning message for small
%                    temperature ranges, consistency of inputs with
%                    quick_look)
%                    - Added functions to estimate the electronic noise for
%                      thermistor and shear channels
%                      from shear and thermistor channels.
%                    - show_ch: Warning messages removed and "convert" option
%                           modified to also prevent a deconvolve from being 
%                           performed.
%                    - convert_odas now has support for type "raw"
%                    - Other minor annotation corrections and documentation
%                    updates. 
%                    

function version = odas_version_info( )

version = 4.4;

end

