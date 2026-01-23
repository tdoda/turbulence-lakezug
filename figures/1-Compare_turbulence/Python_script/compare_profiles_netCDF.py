# Plot VMP or microCTD profiles from a specific campaign

# -*- coding: utf-8 -*-
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), r'..\..\..\functions\figures'))
from functions_figures import array_bands, get_cmap_discrete, compute_logmedian
sys.path.append(os.path.join(os.path.dirname(__file__), r'..\..\..\functions\ctd'))
from functions_ctd import thorpe_scale
sys.path.append(os.path.join(os.path.dirname(__file__), r'..\..\..\functions\general_functions'))
from general_functions import matstruct2dict, read_netCDF_xr, read_netCDF
import netCDF4
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 12
plt.rcParams['font.family'] = 'Arial'  # or any installed font
from scipy.io import loadmat
plt.close('all')

#%% Specify field campaign here:
campaign_name='20260113'
datalevel='Level2'

# Relative path:
data_folder_rel=os.path.join('..','..','..','data','microCTD',campaign_name,datalevel)

#folder_param_name='down1.0_K5.3_L_NAS0.0058_singlepole_nfft3_EPFL_Goodman'
fileprof=['DAT_059_down_prof1','DAT_059_down_prof2'] # all files
profnames=['VMPC1','VMPC2']

savefig=True


#%% Load the data from netCDF files

dataprof=[]
for kf in range(len(fileprof)):
    if datalevel=='Level1':
        filename_full='L1_'+fileprof[kf]+'.nc'
    elif datalevel=='Level2':
        filename_full='L2_'+fileprof[kf]+'.nc'
    else:
        raise Exception('Wrong level')
        
    dataprof.append(read_netCDF_xr(os.path.join(data_folder_rel,filename_full)))




#%% Check Thorpse length calculation in compare_profiles_mat

#%% Figure profiles
# Create a figure with specified size in centimeters
fig,ax = plt.subplots(1,4,figsize=(18 / 2.54, 9 / 2.54),sharey=True)  # Convert from cm to inches


for kf in range(len(dataprof)):
    #ax[0].plot(dataprof["temperature"],dataprof_slow["depth"],linewidth=1)
    
    ax[1].plot(dataprof[kf].FAST_fast_T1.data,dataprof[kf].FAST_depth,linewidth=1)
    
    ax[2].plot(dataprof[kf].FAST_grad_T1.data+2*kf,dataprof[kf].FAST_depth,linewidth=1)
    
    ax[3].plot(dataprof[kf].FAST_fast_S1.data+kf,dataprof[kf].FAST_depth,linewidth=1)


# Labels and limits:
ax[0].set_ylabel('Depth [m]')
ax[0].set_xlabel('CTD Temp [°C]')
ax[1].set_xlabel('FP07 T1 [°C]')
ax[2].set_xlabel('dT1/dz [°C/m]')
ax[3].set_xlabel('Shear S1 [s$^{-1}$]')

ax[2].set_xlim(-1,2*len(dataprof))
ax[3].set_xlim(-1,len(dataprof))
ax[0].invert_yaxis()

fig.set_tight_layout(True)

if savefig:
    fig.savefig('comparison_prof_'+campaign_name+'.png',dpi=400)
    fig.savefig('comparison_prof_'+campaign_name+'.svg',dpi=400)


#%% Figure turbulence
# Create a figure with specified size in centimeters
fig,ax = plt.subplots(1,4,figsize=(18 / 2.54, 9 / 2.54),sharey=True)  # Convert from cm to inches


for kf in range(len(dataprof)):
    indkeep_T1=np.where(dataprof[kf].BINNED_flag_T1.data==0)[0]
    indkeep_S1=np.where(dataprof[kf].BINNED_flag_S1.data==0)[0]
    
    ax[0].plot(dataprof[kf].BINNED_LTuT1.data,dataprof[kf].BINNED_depth.data,'-',linewidth=1)
    
    ax[1].plot(dataprof[kf].BINNED_Xi_T1.data,dataprof[kf].BINNED_depth.data,'k.',linewidth=1)
    ax[1].plot(dataprof[kf].BINNED_Xi_T1.data[indkeep_T1],dataprof[kf].BINNED_depth.data[indkeep_T1],'-',linewidth=1)
    
    ax[2].plot(dataprof[kf].BINNED_eps_T1.data[indkeep_T1],dataprof[kf].BINNED_depth.data[indkeep_T1],'-',linewidth=1)
    
    ax[3].plot(dataprof[kf].BINNED_eps_S1.data[indkeep_S1],dataprof[kf].BINNED_depth.data[indkeep_S1],'-',linewidth=1)
    
    
# Labels and limits:
ax[0].set_ylabel('Depth [m]')
ax[0].set_xlabel('$L_T$ [m]')
ax[1].set_xlabel('$\\chi$ [°C$^2$ s$^{-1}$]')
ax[2].set_xlabel('$\\varepsilon_{T1}$ [W/kg]')
ax[3].set_xlabel('$\\varepsilon_{S1}$ [W/kg]')

xlimval=ax[0].get_xlim()
ax[0].set_xticks(np.arange(0, 100, 10))
ax[0].set_xlim(xlimval)

ax[1].set_xscale('log')
xlimval=ax[1].get_xlim()
ax[1].set_xticks(10.0**np.arange(-14, -5, 1))
ax[1].set_xlim(xlimval)
ax[1].tick_params(axis='x', labelrotation=45)

ax[2].set_xscale('log')
xlimval=ax[2].get_xlim()
ax[2].set_xticks(10.0**np.arange(-14, -5, 1))
ax[2].set_xlim(xlimval)
ax[2].tick_params(axis='x', labelrotation=45)

ax[3].set_xscale('log')
xlimval=ax[3].get_xlim()
ax[3].set_xticks(10.0**np.arange(-14, -5, 1))
ax[3].set_xlim(xlimval)
ax[3].tick_params(axis='x', labelrotation=45)

ax[0].invert_yaxis()

fig.set_tight_layout(True)

if savefig:
    fig.savefig('comparison_turbulence_'+campaign_name+'.png',dpi=400)
    fig.savefig('comparison_turbulence_'+campaign_name+'.svg',dpi=400)

#%% Figure averaged turbulence
# colC=[1, 0.4, 0.15]
# colC_light=[1, 0.75, 0.6]
# colS=[50/255,100/255,214/255]
# colS_light=[140/255,170/255,230/255]
# # Create a figure with specified size in centimeters
# fig,ax = plt.subplots(1,4,figsize=(18 / 2.54, 9 / 2.54),sharey=True)  # Convert from cm to inches

# # 1st subplot: Eps-T1
# depthC_all, epsC_all, depthC_cells, logmed_epsC = compute_logmedian(DATA, 5, 'eps_T1', np.array([1, 2, 5]))
# depthS_all, epsS_all, depthS_cells, logmed_epsS = compute_logmedian(DATA, 5, 'eps_T1', np.array([8]))

# # Plot the data
# ax[0].plot(epsC_all, depthC_all, 'o', markersize=1, markeredgecolor=colC_light, markerfacecolor=colC_light,rasterized=True)
# ax[0].plot(10**logmed_epsC, depthC_cells, '-', color=colC)
# ax[0].plot(epsS_all, depthS_all, 'o', markersize=1, markeredgecolor=colS_light, markerfacecolor=colS_light,rasterized=True)
# ax[0].plot(10**logmed_epsS, depthS_cells, '-',color=colS)
# ax[0].set_ylabel('Depth (m)')
# ax[0].set_xlabel(r'$\epsilon$ (m$^2$ s$^{-3}$)')
# ax[0].set_xscale('log')
# xlimval=ax[0].get_xlim()
# ax[0].set_xticks(10.0**np.arange(-14, 1, 2))
# ax[1].set_xlim(xlimval)
# ax[0].tick_params(axis='x', labelrotation=45)
# ax[0].set_xlim([10**(-11), 10**(-4)])

# # 2nd subplot: X-T1
# depthC_all, XiC_all, depthC_cells, logmed_XiC = compute_logmedian(DATA, 5, 'Xi_T1', np.array([1, 2, 5]))
# depthS_all, XiS_all, depthS_cells, logmed_XiS = compute_logmedian(DATA, 5, 'Xi_T1', np.array([8]))

# ax[1].plot(XiC_all, depthC_all, 'o', markersize=1, markeredgecolor=colC_light, markerfacecolor=colC_light,rasterized=True)
# ax[1].plot(10**logmed_XiC, depthC_cells, '-', color=colC)
# ax[1].plot(XiS_all, depthS_all, 'o', markersize=1, markeredgecolor=colS_light, markerfacecolor=colS_light,rasterized=True)
# ax[1].plot(10**logmed_XiS, depthS_cells, '-',color=colS)
# ax[1].set_xlabel(r'$\chi$ (K$^2$ s$^{-1}$)')
# ax[1].set_xscale('log')
# xlimval=ax[1].get_xlim()
# ax[1].set_xticks(10.0**np.arange(-14, 1, 2))
# ax[1].set_xlim(xlimval)
# ax[1].tick_params(axis='x', labelrotation=45)
# ax[1].set_xlim([10**(-13), 10**(-6)])

# # 3rd subplot: Kz-T1
# depthC_all, KzC_all, depthC_cells, logmed_KzC = compute_logmedian(DATA, 5, 'KOsbornCox_T1', np.array([1, 2, 5]))
# depthS_all, KzS_all, depthS_cells, logmed_KzS = compute_logmedian(DATA, 5, 'KOsbornCox_T1', np.array([8]))

# ax[2].plot(KzC_all, depthC_all, 'o', markersize=1, markeredgecolor=colC_light, markerfacecolor=colC_light,rasterized=True)
# ax[2].plot(10**logmed_KzC, depthC_cells, '-', color=colC)
# ax[2].plot(KzS_all, depthS_all, 'o', markersize=1, markeredgecolor=colS_light, markerfacecolor=colS_light,rasterized=True)
# ax[2].plot(10**logmed_KzS, depthS_cells, '-',color=colS)
# ax[2].set_xlabel(r'$K_{z}$ (m$^2$ s$^{-1}$)')
# ax[2].set_xscale('log')
# xlimval=ax[2].get_xlim()
# ax[2].set_xticks(10.0**np.arange(-14, 1, 2))
# ax[2].set_xlim(xlimval)
# ax[2].tick_params(axis='x', labelrotation=45)

# # 4th subplot: Thorpe-T1
# # depthC_all, LtC_all, depthC_cells, logmed_LtC = compute_logmedian(DATA, 5, 'Lt', np.array([1, 2, 5]),fieldname='THORPE')
# # depthS_all, LtS_all, depthS_cells, logmed_LtS = compute_logmedian(DATA, 5, 'Lt', np.array([8]),fieldname='THORPE')
# # depthC_all, LtC_all, depthC_cells, logmed_LtC = compute_logmedian(DATA, 5, 'LTuT1', np.array([1, 2, 5]))
# # depthS_all, LtS_all, depthS_cells, logmed_LtS = compute_logmedian(DATA, 5, 'LTuT1', np.array([8]))
# depthC_all, LtC_all, depthC_cells, logmed_LtC = compute_logmedian(DATA, 5, 'thorpe_disp_abs', np.array([1, 2, 5]),fieldname='THORPE')
# depthS_all, LtS_all, depthS_cells, logmed_LtS = compute_logmedian(DATA, 5, 'thorpe_disp_abs', np.array([8]),fieldname='THORPE')

# ax[3].plot(LtC_all, depthC_all, 'o', markersize=1, markeredgecolor=colC_light, markerfacecolor=colC_light,rasterized=True)
# ax[3].plot(10**logmed_LtC, depthC_cells, '-', color=colC)
# ax[3].plot(LtS_all, depthS_all, 'o', markersize=1, markeredgecolor=colS_light, markerfacecolor=colS_light,rasterized=True)
# ax[3].plot(10**logmed_LtS, depthS_cells, '-',color=colS)
# ax[3].set_xlabel(r'$d_{T}$ (m)')
# ax[3].set_xscale('log')
# xlimval=ax[3].get_xlim()
# ax[3].set_xticks(10.0**np.arange(-14, 1, 2))
# ax[3].set_xlim(xlimval)
# ax[3].tick_params(axis='x', labelrotation=45)

# ax[0].set_ylim((0,175))
# ax[0].invert_yaxis()  # Reverse the y-axis (like MATLAB's `set(ax1,'ydir','reverse')`)

# fig.set_tight_layout(True)

# # Save figure

# if savefig:
#     fig.savefig('comparison_turbulence_Python.png',dpi=400)
#     fig.savefig('comparison_turbulence_Python.svg',dpi=1000)
