# Plot VMP profiles from 11 February 2025

# -*- coding: utf-8 -*-
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), r'..\..\functions\figures'))
from functions_figures import array_bands, get_cmap_discrete, compute_logmedian
sys.path.append(os.path.join(os.path.dirname(__file__), r'..\..\functions\ctd'))
from functions_ctd import thorpe_scale
sys.path.append(os.path.join(os.path.dirname(__file__), r'..\..\functions\general_functions'))
from general_functions import matstruct2dict
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
plt.close('all')
campaign_name='20250211'
data_folder='../../data/VMP/'+campaign_name+'/Level2/'
folder_param_name='down1.0_K5.3_L_NAS0.0058_singlepole_nfft3_EPFL_Goodman'
filenames=['VMP003','VMP004','VMP007','VMP008','VMP009','VMP010','VMP012','VMP013','VMP014'] # all files
profnames=['PSB01','PC01','PC02','PC03','PC04','PC05','PC06','PC07','PSB02']

savefig=False

#%% Load the data from the matfile and computes Thorpe scale with method based on mixsea

DATA=[None]*len(filenames)
for kf,file in enumerate(filenames):
    datavar=loadmat(data_folder+file+'_'+folder_param_name+'/results_'+file+'_down.mat',squeeze_me=True, struct_as_record=False)
    DATA[kf]=dict()
    DATA[kf]["filename"]=file
    if type(datavar["BINNED"])==np.ndarray: # Several profiles
        DATA[kf]["BINNED"]=[None]*len(datavar["BINNED"])
        for kprof in range(len(datavar["BINNED"])):
            DATA[kf]["BINNED"][kprof]=matstruct2dict(datavar["BINNED"][kprof])  
    else:
            DATA[kf]["BINNED"]=matstruct2dict(datavar["BINNED"])
            
    if type(datavar["DISS_QL"])==np.ndarray: # Several profiles
        DATA[kf]["DISS_QL"]=[None]*len(datavar["DISS_QL"])
        for kprof in range(len(datavar["DISS_QL"])):
            DATA[kf]["DISS_QL"][kprof]=matstruct2dict(datavar["DISS_QL"][kprof])  
    else:
            DATA[kf]["DISS_QL"]=matstruct2dict(datavar["DISS_QL"])
            
    if type(datavar["FAST"])==np.ndarray: # Several profiles
        DATA[kf]["FAST"]=[None]*len(datavar["FAST"])
        for kprof in range(len(datavar["FAST"])):
            DATA[kf]["FAST"][kprof]=matstruct2dict(datavar["FAST"][kprof])  
    else:
            DATA[kf]["FAST"]=matstruct2dict(datavar["FAST"])
            
    if type(datavar["SLOW"])==np.ndarray: # Several profiles
        DATA[kf]["SLOW"]=[None]*len(datavar["SLOW"])
        for kprof in range(len(datavar["SLOW"])):
            DATA[kf]["SLOW"][kprof]=matstruct2dict(datavar["SLOW"][kprof])  
    else:
            DATA[kf]["SLOW"]=matstruct2dict(datavar["SLOW"])
            
    DATA[kf]["param"]=matstruct2dict(datavar["param"])
    
    # Remove the bottom and top to avoid getting large overturns if temperature peaks:
    dz_botrem=1
    dz_toprem=10
    depthval=DATA[kf]["FAST"]["depth"]
    tempval=DATA[kf]["FAST"]["fast_T1"]
    tempval[(depthval<np.min(depthval)+dz_toprem)|(depthval>np.max(depthval)-dz_botrem)]=np.nan
    ind_nan=np.isnan(tempval)
    
    Lt_nan=np.full(depthval.shape,np.nan)
    thorpe_disp_nan=np.full(depthval.shape,np.nan)
    temp_sorted_nan=np.full(depthval.shape,np.nan)
    
    
    Lt, thorpe_disp, temp_sorted, _, idx_overturns, idx_sorted=thorpe_scale(depthval[~ind_nan],tempval[~ind_nan],"decreasing",res=0.005)
    # Lt, thorpe_disp, temp_sorted, _, idx_overturns, idx_sorted=thorpe_scale(DATA[kf]["SLOW"]["depth"],DATA[kf]["SLOW"]["temperature"],"decreasing")
    
    Lt_nan[~ind_nan]=Lt
    thorpe_disp_nan[~ind_nan]=thorpe_disp
    temp_sorted_nan[~ind_nan]=temp_sorted
    idx_overturns=np.where(~ind_nan)[0][idx_overturns]
    idx_sorted=np.where(~ind_nan)[0][idx_sorted]
    
    overturn_bool=np.zeros((len(depthval),))
    for ko in range(len(idx_overturns)):
        overturn_bool[idx_overturns[ko][0]:idx_overturns[ko][1]+1]=1 
    

    DATA[kf]["THORPE"]={"depth":depthval,"temp_sorted":temp_sorted_nan,"idx_overturns":idx_overturns,
                        "overturn_bool":overturn_bool,"thorpe_disp":thorpe_disp_nan,"thorpe_disp_abs":np.abs(thorpe_disp_nan),"Lt":Lt_nan}

#%% --> COULD BE CONVERTED TO NETCDF HERE (AND MOVING THE SCRIPT TO THE ANALYSIS FOLDER)
# THE FOLLOWING PART IS RELATED TO FIGURES THAT COULD BE DONE AFTER IMPORTING THE NETCDF FILE

#%% Figure
colC=[1, 0.4, 0.15]
colC_light=[1, 0.75, 0.6]
colS=[50/255,100/255,214/255]
colS_light=[140/255,170/255,230/255]
# Create a figure with specified size in centimeters
fig,ax = plt.subplots(1,4,figsize=(18 / 2.54, 9 / 2.54),sharey=True)  # Convert from cm to inches

# 1st subplot: Eps-T1
depthC_all, epsC_all, depthC_cells, logmed_epsC = compute_logmedian(DATA, 5, 'eps_T1', np.array([1, 2, 5]))
depthS_all, epsS_all, depthS_cells, logmed_epsS = compute_logmedian(DATA, 5, 'eps_T1', np.array([8]))

# Plot the data
ax[0].plot(epsC_all, depthC_all, 'o', markersize=1, markeredgecolor=colC_light, markerfacecolor=colC_light,rasterized=True)
ax[0].plot(10**logmed_epsC, depthC_cells, '-', color=colC)
ax[0].plot(epsS_all, depthS_all, 'o', markersize=1, markeredgecolor=colS_light, markerfacecolor=colS_light,rasterized=True)
ax[0].plot(10**logmed_epsS, depthS_cells, '-',color=colS)
ax[0].set_ylabel('Depth (m)')
ax[0].set_xlabel(r'$\epsilon$ (m$^2$ s$^{-3}$)')
ax[0].set_xscale('log')
xlimval=ax[0].get_xlim()
ax[0].set_xticks(10.0**np.arange(-14, 1, 2))
ax[1].set_xlim(xlimval)
ax[0].tick_params(axis='x', labelrotation=45)
ax[0].set_xlim([10**(-11), 10**(-4)])

# 2nd subplot: X-T1
depthC_all, XiC_all, depthC_cells, logmed_XiC = compute_logmedian(DATA, 5, 'Xi_T1', np.array([1, 2, 5]))
depthS_all, XiS_all, depthS_cells, logmed_XiS = compute_logmedian(DATA, 5, 'Xi_T1', np.array([8]))

ax[1].plot(XiC_all, depthC_all, 'o', markersize=1, markeredgecolor=colC_light, markerfacecolor=colC_light,rasterized=True)
ax[1].plot(10**logmed_XiC, depthC_cells, '-', color=colC)
ax[1].plot(XiS_all, depthS_all, 'o', markersize=1, markeredgecolor=colS_light, markerfacecolor=colS_light,rasterized=True)
ax[1].plot(10**logmed_XiS, depthS_cells, '-',color=colS)
ax[1].set_xlabel(r'$\chi$ (K$^2$ s$^{-1}$)')
ax[1].set_xscale('log')
xlimval=ax[1].get_xlim()
ax[1].set_xticks(10.0**np.arange(-14, 1, 2))
ax[1].set_xlim(xlimval)
ax[1].tick_params(axis='x', labelrotation=45)
ax[1].set_xlim([10**(-13), 10**(-6)])

# 3rd subplot: Kz-T1
depthC_all, KzC_all, depthC_cells, logmed_KzC = compute_logmedian(DATA, 5, 'KOsbornCox_T1', np.array([1, 2, 5]))
depthS_all, KzS_all, depthS_cells, logmed_KzS = compute_logmedian(DATA, 5, 'KOsbornCox_T1', np.array([8]))

ax[2].plot(KzC_all, depthC_all, 'o', markersize=1, markeredgecolor=colC_light, markerfacecolor=colC_light,rasterized=True)
ax[2].plot(10**logmed_KzC, depthC_cells, '-', color=colC)
ax[2].plot(KzS_all, depthS_all, 'o', markersize=1, markeredgecolor=colS_light, markerfacecolor=colS_light,rasterized=True)
ax[2].plot(10**logmed_KzS, depthS_cells, '-',color=colS)
ax[2].set_xlabel(r'$K_{z}$ (m$^2$ s$^{-1}$)')
ax[2].set_xscale('log')
xlimval=ax[2].get_xlim()
ax[2].set_xticks(10.0**np.arange(-14, 1, 2))
ax[2].set_xlim(xlimval)
ax[2].tick_params(axis='x', labelrotation=45)

# 4th subplot: Thorpe-T1
# depthC_all, LtC_all, depthC_cells, logmed_LtC = compute_logmedian(DATA, 5, 'Lt', np.array([1, 2, 5]),fieldname='THORPE')
# depthS_all, LtS_all, depthS_cells, logmed_LtS = compute_logmedian(DATA, 5, 'Lt', np.array([8]),fieldname='THORPE')
# depthC_all, LtC_all, depthC_cells, logmed_LtC = compute_logmedian(DATA, 5, 'LTuT1', np.array([1, 2, 5]))
# depthS_all, LtS_all, depthS_cells, logmed_LtS = compute_logmedian(DATA, 5, 'LTuT1', np.array([8]))
depthC_all, LtC_all, depthC_cells, logmed_LtC = compute_logmedian(DATA, 5, 'thorpe_disp_abs', np.array([1, 2, 5]),fieldname='THORPE')
depthS_all, LtS_all, depthS_cells, logmed_LtS = compute_logmedian(DATA, 5, 'thorpe_disp_abs', np.array([8]),fieldname='THORPE')

ax[3].plot(LtC_all, depthC_all, 'o', markersize=1, markeredgecolor=colC_light, markerfacecolor=colC_light,rasterized=True)
ax[3].plot(10**logmed_LtC, depthC_cells, '-', color=colC)
ax[3].plot(LtS_all, depthS_all, 'o', markersize=1, markeredgecolor=colS_light, markerfacecolor=colS_light,rasterized=True)
ax[3].plot(10**logmed_LtS, depthS_cells, '-',color=colS)
ax[3].set_xlabel(r'$d_{T}$ (m)')
ax[3].set_xscale('log')
xlimval=ax[3].get_xlim()
ax[3].set_xticks(10.0**np.arange(-14, 1, 2))
ax[3].set_xlim(xlimval)
ax[3].tick_params(axis='x', labelrotation=45)

ax[0].set_ylim((0,175))
ax[0].invert_yaxis()  # Reverse the y-axis (like MATLAB's `set(ax1,'ydir','reverse')`)

fig.set_tight_layout(True)

#%% Save figure

if savefig:
    fig.savefig('comparison_turbulence_Python.png',dpi=400)
    fig.savefig('comparison_turbulence_Python.svg',dpi=1000)
