# Convert turbulence profiles from a campaign to netCDF files

# -*- coding: utf-8 -*-
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), r'..\..\functions\figures'))
from functions_figures import array_bands, get_cmap_discrete, compute_logmedian
sys.path.append(os.path.join(os.path.dirname(__file__), r'..\..\functions\ctd'))
from functions_ctd import thorpe_scale
sys.path.append(os.path.join(os.path.dirname(__file__), r'..\..\functions\general_functions'))
from general_functions import matstruct2dict, export_netCDF
import netCDF4
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 12
plt.rcParams['font.family'] = 'Arial'  # or any installed font
from scipy.io import loadmat


#%% Specify field campaign here:
campaign_name='20260113'

# Relative path (doesn't work):
#data_folder=r'../../../data/microCTD/'+campaign_name+'/Level2/'

# Other option without "..":
script_path=os.path.dirname(os.path.abspath(__file__))
data_folder=os.path.join(os.path.dirname(os.path.dirname(script_path)),'data','microCTD',campaign_name,'Level2')

folder_param_name='down1.0_K5.3_L_NAS0.0058_singlepole_nfft3_EPFL_Goodman'
filenames=['DAT_053','DAT_055','DAT_057','DAT_059'] # all files
profnames=['VMPS','VMPC','VMPW','VMPC']

# Thorpe scale: remove the bottom and top to avoid getting large overturns if temperature peaks:
dz_botrem=1 # m
dz_toprem=10 # m


#%% Load the data from the matfile and convert it to dictionary 

DATA=[None]*len(filenames)
for kf,file in enumerate(filenames):
    datavar=loadmat(os.path.join(data_folder,file+'_'+folder_param_name,'L2_'+file+'_down.mat'),squeeze_me=True, struct_as_record=False)

    if type(datavar["BINNED"])==np.ndarray: # Several profiles
        nprof=len(datavar["BINNED"])        
    else:
        nprof=1
        datavar["BINNED"]=np.array([datavar["BINNED"]])
        datavar["SLOW"]=np.array([datavar["SLOW"]])
        datavar["FAST"]=np.array([datavar["FAST"]])
        datavar["DISS_QL"]=np.array([datavar["DISS_QL"]])
        
    for kprof in range(nprof):
        dict_slow=matstruct2dict(datavar["SLOW"][kprof],'SLOW_')
        dict_fast=matstruct2dict(datavar["FAST"][kprof],'FAST_')
        if len(datavar["DISS_QL"])>0:
            dict_ql=matstruct2dict(datavar["DISS_QL"][kprof],'DISS_QL_')
        else:
            dict_ql=dict()
        dict_binned=matstruct2dict(datavar["BINNED"][kprof],'BINNED_')
        dict_param=matstruct2dict(datavar["param"],'PARAM_') # Could become arguments of netCDF
        
  
        # Thorpe scale with same method as for CTD (based on mixsea)
        depthval=dict_fast["FAST_depth"].copy()
        tempval=dict_fast["FAST_fast_T1"].copy()
        tempval[(depthval<np.min(depthval)+dz_toprem)|(depthval>np.max(depthval)-dz_botrem)]=np.nan
        ind_nan=np.isnan(tempval)
        
        Lt_nan=np.full(depthval.shape,np.nan)
        thorpe_disp_nan=np.full(depthval.shape,np.nan)
        temp_sorted_nan=np.full(depthval.shape,np.nan)
          
        Lt, thorpe_disp, temp_sorted, _, idx_overturns, idx_sorted=thorpe_scale(depthval[~ind_nan],tempval[~ind_nan],"decreasing",res=0.005)
        
        Lt_nan[~ind_nan]=Lt
        thorpe_disp_nan[~ind_nan]=thorpe_disp
        temp_sorted_nan[~ind_nan]=temp_sorted
        idx_overturns=np.where(~ind_nan)[0][idx_overturns]
        idx_sorted=np.where(~ind_nan)[0][idx_sorted]
        
        overturn_bool=np.zeros((len(depthval),))
        for ko in range(len(idx_overturns)):
            overturn_bool[idx_overturns[ko][0]:idx_overturns[ko][1]+1]=1 
        
    
        dict_thorpe={"FAST_temp_sorted":temp_sorted_nan,"FAST_idx_overturns":idx_overturns,
                            "FAST_overturn_bool":overturn_bool,"FAST_thorpe_disp":thorpe_disp_nan,"FAST_thorpe_disp_abs":np.abs(thorpe_disp_nan),"FAST_Lt":Lt_nan}
        
        # Combine dictionaries
        data_prof=dict_slow|dict_fast|dict_thorpe|dict_ql|dict_binned 
        
        breakpoint()
        # Remove quantities that do not have consistent dimensions: 
        # could be done by specifying in a class/function which variables are exported (implying that export netCDF function should also be in that class)

