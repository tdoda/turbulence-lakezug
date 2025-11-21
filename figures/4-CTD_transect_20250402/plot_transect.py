# Plot CTD quantities along transect from 2nd April 2025

# -*- coding: utf-8 -*-
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), r'..\..\functions\figures'))
from functions_figures import array_bands, get_cmap_discrete
sys.path.append(os.path.join(os.path.dirname(__file__), r'..\..\functions\general_functions'))
from general_functions import netCDF2dict, bathy_transect
import netCDF4
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import cmocean
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 12
plt.rcParams['font.family'] = 'Arial'  # or any installed font

#%% Specify field campaign here:
plt.close('all')
date_campaign='20250402'
data_folder='../../analysis/3-Process_transect'
filename='CTD_transect_'+date_campaign+'.nc'

savefig=False
#%% Import the CTD data
nc=netCDF4.Dataset(os.path.join(data_folder,filename))
CTD_data,CTD_genatt,CTD_var=netCDF2dict(nc)
nc.close()

#%% Plot all profiles
fig,ax=plt.subplots(1,2,sharey=True)

for k in range(len(CTD_data["time"])):
    ax[0].plot(CTD_data["temp"][:,k],CTD_data["depth"])
    
    ax[1].plot(CTD_data["rho"][:,k],CTD_data["depth"])
    
ax[0].invert_yaxis()
ax[0].set_xlabel("Temp [°C]")
ax[0].set_ylabel("Depth [m]")
ax[1].set_xlabel(r"$\rho$ [kg/m$^3$]")

#%% Plot Thorpe for a given profile
k_time=11

fig,ax=plt.subplots(1,2,sharey=True)
hp1,=ax[0].plot(CTD_data["temp"][:,k_time],CTD_data["depth"])
hp2,=ax[0].plot(CTD_data["temp_sorted"][:,k_time],CTD_data["depth"],'-k')
show_overturn=np.copy(CTD_data["temp"][:,k_time])
show_overturn[CTD_data["overturn_bool"][:,k_time]==0]=np.nan
hp3,=ax[0].plot(show_overturn,CTD_data["depth"],'-r')
ax[0].set_xlabel("Temp [°C]")
ax[0].set_ylabel("Depth [m]")
ax_top = ax[0].twiny()
hp4,=ax_top.plot(CTD_data["DO"][:,k_time],CTD_data["depth"])
ax_top.set_xlabel("DO [mg/l]")
ax[0].legend([hp2,hp3,hp4],("Sorted Temp","Overturns","DO"))

ax[1].plot(np.abs(CTD_data["thorpe_disp"][:,k_time]),CTD_data["depth"],'-k')
ax[1].plot(CTD_data["Lt"][:,k_time],CTD_data["depth"],'-r')
#ax[1].set_xscale('log')
ax[1].set_xlabel("Displacements [m]")
ax[1].legend(("Thorpe displ","Thorpe scale"))

fig.suptitle("Prof {} (dist {} m)".format(CTD_data["profnames"][k_time],CTD_data["disttrans"][k_time]))
    
ax[0].invert_yaxis()

#%% Load bathymetry data
nc=netCDF4.Dataset('../../data/Bathymetry_DHM25/bathy_Zugersee_5m.nc')
bathy_data,_,_=netCDF2dict(nc)
nc.close()
bathy_data["GRID5m_depth"]=bathy_data["GRID5m_depth"].transpose()
xbathy,ybathy,dist_bathy,depth_bathy=bathy_transect(bathy_data["GRID5m_x"],bathy_data["GRID5m_y"],bathy_data["GRID5m_depth"],xtrans=np.nan,ytrans=np.nan,monotonic_dir='y',dx=1)

#%% Figure map
#indprof=np.arange(12)# transect 1
# indprof=np.concatenate((np.arange(12,15),np.arange(-4,0))) # transect 2
#indprof=np.arange(15,21) # perpendicular transect
# indprof=np.array([-4,7,-3])
indprof=np.concatenate((np.arange(4),np.arange(6,15),np.arange(-4,0))) # transect 1+2 (without 4)


fig_map,ax=plt.subplots(subplot_kw={'projection': '3d'})
X,Y=np.meshgrid(bathy_data["GRID5m_x"],bathy_data["GRID5m_y"])
#surf = ax.plot_surface(X,Y,-bathy_data["GRID5m_depth"], cmap='viridis', alpha=0.6)
ax.plot(bathy_data["contour"][0,:],bathy_data["contour"][1,:],np.zeros(bathy_data["contour"][1,:].shape),'-k',linewidth=0.5)
ax.view_init(elev=20, azim=-90-45)

# Set same scale:
xlim = ax.get_xlim()
ylim = ax.get_ylim()
# Find the maximum range to use for all axes
max_range = max(np.ptp(xlim), np.ptp(ylim))
# Set the new limits
ax.set_xlim(np.mean(xlim) - max_range / 2, np.mean(xlim) + max_range / 2)
ax.set_ylim(np.mean(ylim) - max_range / 2, np.mean(ylim) + max_range / 2)
for k in range(len(indprof)):
    ax.plot(CTD_data["xcoord"][indprof[k]],CTD_data["ycoord"][indprof[k]],0,'o',markersize=3,markeredgecolor='none',markerfacecolor='k')



# 2D plot:
# plt.plot(bathy_data["contour"][0,:],bathy_data["contour"][1,:],'-k')
# plt.scatter(x=CTD_data["xcoord"][indprof],y=CTD_data["ycoord"][indprof],c=CTD_data["time"][indprof])
# plt.axis('equal')

# Export the figures
if savefig:  
    fig_map.savefig('map_3d.svg',dpi=400)

#%% Find oxycline
DO_limit=2 # [mg/l]
depth_oxycline=np.full(CTD_data["time"].shape,np.nan)
for kp in range(len(CTD_data["time"])):
    ind_oxycline=np.where(CTD_data["DO"][:,kp]<=DO_limit)[0]
    if len(ind_oxycline)>0:
        depth_oxycline[kp]=CTD_data["depth"][ind_oxycline[0]]
#%% Plot grid as a function of distance
dx=400 # Band width [m]
# indprof=np.arange(12)# transect 1
# indprof=np.concatenate((np.arange(12,15),np.arange(-4,0))) # transect 2
indprof=np.concatenate((np.arange(4),np.arange(6,15),np.arange(-4,0))) # transect 1+2 (without 4)
# indprof=np.arange(15,21) # perpendicular transect

ind_sort_dist=np.argsort(CTD_data["disttrans"][indprof])
distsort=CTD_data["disttrans"][indprof][ind_sort_dist]
tempmat=CTD_data["temp"][:,indprof][:,ind_sort_dist]
Ltmat=CTD_data["Lt"][:,indprof][:,ind_sort_dist]
DOmat=CTD_data["DO"][:,indprof][:,ind_sort_dist]
distmat,tempmat=array_bands(distsort,tempmat,dx)
_,Ltmat=array_bands(distsort,Ltmat,dx)
_,DOmat=array_bands(distsort,DOmat,dx)

fig_grid,ax=plt.subplots(2,1,sharey=True,sharex=True,figsize=(18/2.54,9/2.54))


#ax[0].pcolormesh(CTD_data["disttrans"][indprof][ind_sort_dist],CTD_data["depth"],CTD_data["temp"][:,indprof][:,ind_sort_dist])
# hp1=ax[0].pcolormesh(distmat/1000,CTD_data["depth"],tempmat)
# ax[0].plot(dist_bathy/1000,depth_bathy,'k-')
# fig.colorbar(hp1)

hp2=ax[0].pcolormesh(distmat/1000,CTD_data["depth"],np.log10(Ltmat),rasterized=True)
ax[0].plot(dist_bathy/1000,depth_bathy,'k-')
ax[0].plot(distsort[np.array([0,-1])]/1000,[110,110],'k-') # Depth limit for profile plot
ax[0].set_ylabel('Depth (m)')
cb1=fig_grid.colorbar(hp2)
cb1.set_label("$\log_{10}(L_T)$ (-)")

DOmol=DOmat/32*1000
# hp3=ax[1].pcolormesh(distmat/1000,CTD_data["depth"],DOmat,cmap=cmocean.cm.oxy,vmin=0,vmax=2,rasterized=True)
hp3=ax[1].pcolormesh(distmat/1000,CTD_data["depth"],DOmol,cmap=cmocean.cm.oxy,vmin=0,vmax=60,rasterized=True)
ax[1].plot(dist_bathy/1000,depth_bathy,'k-')
ax[1].plot(distsort/1000,depth_oxycline[indprof[ind_sort_dist]],'r.-')
ax[1].set_ylabel('Depth (m)')
ax[1].set_xlabel('Transect distance (km)')
cb2=fig_grid.colorbar(hp3)
cb2.set_label("DO ($\mu$mol l$^{-1}$)")

ax[0].set_yticks(ticks=np.arange(0,200,50))
ax[0].invert_yaxis()

fig_grid.set_tight_layout(True)

# Export the figures
if savefig:
    fig_grid.savefig('transect_grid.png',dpi=400)
    fig_grid.savefig('transect_grid.svg',dpi=2000)

#%% Plot grid as a function of profile number
dx=0.5 # Band width [m]
# indprof=np.arange(len(CTD_data["time"]))
# indprof=np.arange(12)# transect 1
# indprof=np.concatenate((np.arange(12,15),np.arange(-4,0))) # transect 2
indprof=np.concatenate((np.arange(4),np.arange(6,15),np.arange(-4,0))) # transect 1+2 (without 4)
#indprof=np.arange(15,21) # perpendicular transect

ind_sort_dist=np.argsort(CTD_data["disttrans"][indprof])
distsort=CTD_data["disttrans"][indprof][ind_sort_dist]
tempmat=CTD_data["temp"][:,indprof][:,ind_sort_dist]
Ltmat=CTD_data["Lt"][:,indprof][:,ind_sort_dist]
DOmat=CTD_data["DO"][:,indprof][:,ind_sort_dist]
turbmat=CTD_data["turb"][:,indprof][:,ind_sort_dist]
indmat,tempmat=array_bands(np.arange(len(distsort)),tempmat,dx)
_,Ltmat=array_bands(np.arange(len(distsort)),Ltmat,dx)
_,DOmat=array_bands(np.arange(len(distsort)),DOmat,dx)
_,turbmat=array_bands(np.arange(len(distsort)),turbmat,dx)

fig,ax=plt.subplots(3,1,sharey=True,sharex=True)

hp1=ax[0].pcolormesh(indmat,CTD_data["depth"],tempmat)
fig.colorbar(hp1)

hp2=ax[1].pcolormesh(indmat,CTD_data["depth"],np.log10(Ltmat))
fig.colorbar(hp2)

hp3=ax[2].pcolormesh(indmat,CTD_data["depth"],DOmat,vmin=0,vmax=2)
fig.colorbar(hp3)

ax[0].invert_yaxis()

#%% Plot shifted profiles of profile number
dT=0.1 # Temperature shift [°C]
dlog_Lt=3 # log10(Lt) shift 
dDO=1.5 # DO shift [mg/l]
dDOmol=45 # DO shift [umol/l]
drho=0.05# density shift [kg/m3]
dturb=0.5
dmin=115 # min depth [m]
# dmin=0 # min depth [m]
ind_depth=np.where(CTD_data["depth"]>dmin)[0]
# colprof=plt.get_cmap('tab10',len(ind_sort_dist))
colprof=get_cmap_discrete('tab10',len(ind_sort_dist))

# indprof=np.arange(len(CTD_data["time"]))
# indprof=np.arange(12)# transect 1
# indprof=np.concatenate((np.arange(12,15),np.arange(-4,0))) # transect 2
# indprof=np.concatenate((np.arange(4),np.arange(6,15),np.arange(-4,0))) # transect 1+2 (without 4)
#indprof=np.arange(15,21) # perpendicular transect
indprof=np.concatenate((np.arange(6,11),np.arange(14,15),np.arange(-4,0))) # keep only the deep profiles

ind_sort_dist=np.argsort(CTD_data["disttrans"][indprof])
distsort=CTD_data["disttrans"][indprof][ind_sort_dist]
tempmat=CTD_data["temp"][:,indprof][:,ind_sort_dist]
tempsortmat=CTD_data["temp_sorted"][:,indprof][:,ind_sort_dist]
overturnmat=CTD_data["overturn_bool"][:,indprof][:,ind_sort_dist]
Ltmat=CTD_data["Lt"][:,indprof][:,ind_sort_dist]
DOmat=CTD_data["DO"][:,indprof][:,ind_sort_dist]
rhomat=CTD_data["rho"][:,indprof][:,ind_sort_dist]
turbmat=CTD_data["turb"][:,indprof][:,ind_sort_dist]
DOmol=DOmat/32*1000

fig_prof,ax=plt.subplots(2,1,figsize=(18/2.54,13/2.54),sharey=True)

for k in range(len(ind_sort_dist)):

    ax[0].plot(tempmat[ind_depth,k]+dT*k,CTD_data["depth"][ind_depth],'-',color=colprof(k))
    ax[0].plot(tempsortmat[ind_depth,k]+dT*k,CTD_data["depth"][ind_depth],'k:')
    #ax[0].plot(tempmat[ind_depth[overturnmat[ind_depth,k]==1],k]+dT*k,CTD_data["depth"][ind_depth[overturnmat[ind_depth,k]==1]],'-r')
    ax[0].set_ylabel('Depth (m)')
    ax[0].set_xlabel('Temperature (°C)')
    
    
    #ax[1].plot(np.log10(Ltmat[ind_depth,k])+dlog_Lt*k,CTD_data["depth"][ind_depth],'-',color=colprof(k)) 
    #ax[1].plot(rhomat[bool_depth,k]+drho*k,CTD_data["depth"][bool_depth],'-k')
    # ax[1].plot(DOmat[ind_depth,k]+dDO*k,CTD_data["depth"][ind_depth],'-',color=colprof(k))
    ax[1].plot(DOmol[ind_depth,k]+dDOmol*k,CTD_data["depth"][ind_depth],'-',color=colprof(k))
    ax[1].plot([dDOmol*k,dDOmol*k],CTD_data["depth"][ind_depth[np.array([0,-1])]],'--',color=colprof(k))
    ax[1].set_ylabel('Depth (m)')
    # ax[1].set_xlabel('DO (mg l$^{-1}$)')
    ax[1].set_xlabel("DO ($\mu$mol l$^{-1}$)")
    ax[1].set_xticks(ticks=np.arange(0,450,50))

ax[0].set_ylim((dmin,175))
ax[0].invert_yaxis()
fig_prof.set_tight_layout(True)

# Export the figure
if savefig:
    fig_prof.savefig('transect_prof.png',dpi=400)
    fig_prof.savefig('transect_prof.svg',dpi=400)

