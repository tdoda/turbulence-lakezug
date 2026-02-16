# Computes Thorpe scale from CTD profiles (note that Thorpe scale is already computed in ctd but with potential temperature)
# and interpolate the transect data to grid

# -*- coding: utf-8 -*-
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), r'..\..\functions\ctd'))
from functions_ctd import thorpe_scale
sys.path.append(os.path.join(os.path.dirname(__file__), r'..\..\functions\general_functions'))
from general_functions import export_netCDF, dist_transect
import netCDF4
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt

#%% Specify field campaign here:
plt.close('all')
date_campaign='20260113' # date of the campaign (format YYYYMMDD)
ctd_data_folder=os.path.join('../../../data/CTD',date_campaign)
data_folder=os.path.join(ctd_data_folder,'Level2')
filenames=os.listdir(data_folder)

cut_top=0.8 # distance to remove from upper part of the profiles [m]
dz_grid=0.01 # [m]

#%% Import the CTD data
CTD_dict=[None]*len(filenames)
Thorpe_temp_dict=[None]*len(filenames)
profnames=[None]*len(filenames)
xcoord=[None]*len(filenames)
ycoord=[None]*len(filenames)
disttrans=[None]*len(filenames)
for k,file in enumerate(filenames):
    nc=netCDF4.Dataset(os.path.join(data_folder,file))
    profnames[k]=getattr(nc,"Profile name")
    depthval=nc.variables["depth"][:].data[nc.variables["depth"][:]>cut_top]
    indsort=np.argsort(depthval)
    xcoord[k]=float(getattr(nc,"X Coordinate (CH1903)"))
    ycoord[k]=float(getattr(nc,"Y Coordinate (CH1903)"))
    disttrans[k]=round(dist_transect(xcoord[k],ycoord[k],xtrans=np.nan,ytrans=np.nan,monotonic_dir='y',method='projected'))
    CTD_dict[k]={"timeall":nc.variables["time"][:].data[nc.variables["depth"][:]>cut_top][indsort],
                         "depth":depthval[indsort],
                         "temp":nc.variables["Temp"][:].data[nc.variables["depth"][:]>cut_top][indsort],
                         "cond":nc.variables["Cond"][:].data[nc.variables["depth"][:]>cut_top][indsort],
                         "rho":nc.variables["rho"][:].data[nc.variables["depth"][:]>cut_top][indsort],
                         "DO":nc.variables["DO_mg"][:].data[nc.variables["depth"][:]>cut_top][indsort],
                         "DOsat":nc.variables["sat"][:].data[nc.variables["depth"][:]>cut_top][indsort],
                         "turb":nc.variables["Turb"][:].data[nc.variables["depth"][:]>cut_top][indsort],
                         "chl":nc.variables["Chl_A"][:].data[nc.variables["depth"][:]>cut_top][indsort]}
    
    Lt, thorpe_disp, temp_sorted, _, idx_overturns, idx_sorted=thorpe_scale(CTD_dict[k]["depth"],CTD_dict[k]["temp"],"decreasing")
    Thorpe_temp_dict[k]={"temp_sorted":temp_sorted,"idx_overturns":idx_overturns,"thorpe_disp":thorpe_disp,"Lt":Lt}
    nc.close()
    

#%% Grid interpolation

depthgrid=np.array([np.arange(0,200,dz_grid),]).transpose()
CTDdata={"time":np.full(len(CTD_dict),np.nan),
         "depth":depthgrid}
varnames=["temp","cond","rho","DO","DOsat","turb","chl","temp_sorted","thorpe_disp","Lt","overturn_bool"]
for var in varnames:
    CTDdata[var]=np.full((len(depthgrid),len(CTD_dict)),np.nan)
CTDdata["xcoord"]=np.full(len(CTD_dict),np.nan)
CTDdata["ycoord"]=np.full(len(CTD_dict),np.nan)
CTDdata["disttrans"]=np.full(len(CTD_dict),np.nan)

for k in range(len(CTD_dict)):
    CTDdata["time"][k]=CTD_dict[k]["timeall"][0]
    for var in varnames:
        if var in CTD_dict[k].keys():
            CTDdata[var][:,k]=np.interp(CTDdata["depth"][:,0],CTD_dict[k]["depth"],CTD_dict[k][var],left=np.nan, right=np.nan)
        if var in Thorpe_temp_dict[k].keys():
            CTDdata[var][:,k]=np.interp(CTDdata["depth"][:,0],CTD_dict[k]["depth"],Thorpe_temp_dict[k][var],left=np.nan, right=np.nan)
    
    overturn_bool=np.zeros((len(depthgrid),))
    for ko in range(len(Thorpe_temp_dict[k]["idx_overturns"])):
        ind_start=np.where(depthgrid>=CTD_dict[k]["depth"][Thorpe_temp_dict[k]["idx_overturns"][ko][0]])[0][0]                          
        ind_end=np.where(depthgrid<=CTD_dict[k]["depth"][Thorpe_temp_dict[k]["idx_overturns"][ko][1]])[0][-1]
        overturn_bool[ind_start:ind_end+1]=1 
    CTDdata["overturn_bool"][:,k]=overturn_bool 
    
    CTDdata["xcoord"][k]=xcoord[k]
    CTDdata["ycoord"][k]=ycoord[k]
    CTDdata["disttrans"][k]=disttrans[k]
    
# Sort in chronological order
ind_sort_chrono=np.argsort(CTDdata["time"])
CTDdata["time"]=CTDdata["time"][ind_sort_chrono]
for var in varnames:
    CTDdata[var]=CTDdata[var][:,ind_sort_chrono]
CTDdata["xcoord"]=CTDdata["xcoord"][ind_sort_chrono]
CTDdata["ycoord"]=CTDdata["ycoord"][ind_sort_chrono]
CTDdata["disttrans"]=CTDdata["disttrans"][ind_sort_chrono]
CTDdata["profnames"]=[profnames[k] for k in ind_sort_chrono]
    
#%% Plot all profiles
fig,ax=plt.subplots(1,2,sharey=True)

for k in range(len(filenames)):
    ax[0].plot(CTD_dict[k]["temp"],CTD_dict[k]["depth"])
    
    ax[1].plot(CTD_dict[k]["rho"],CTD_dict[k]["depth"])
    
ax[0].invert_yaxis()
ax[0].set_xlabel("Temp [°C]")
ax[0].set_ylabel("Depth [m]")
ax[1].set_xlabel(r"$\rho$ [kg/m$^3$]")

#%% Plot Thorpe for a given profile
#Raw data
k_time=0
k_file=ind_sort_chrono[k_time]

fig,ax=plt.subplots(1,2,sharey=True)
hp1,=ax[0].plot(CTD_dict[k_file]["temp"],CTD_dict[k_file]["depth"])
hp2,=ax[0].plot(Thorpe_temp_dict[k_file]["temp_sorted"],CTD_dict[k_file]["depth"],'-k')
for ko in range(len(Thorpe_temp_dict[k_file]["idx_overturns"])):
    ind_overturn=np.arange(Thorpe_temp_dict[k_file]["idx_overturns"][ko][0],Thorpe_temp_dict[k_file]["idx_overturns"][ko][1]+1)
    hp3,=ax[0].plot(CTD_dict[k_file]["temp"][ind_overturn],CTD_dict[k_file]["depth"][ind_overturn],'-r')
ax[0].set_xlabel("Temp [°C]")
ax[0].set_ylabel("Depth [m]")
ax_top = ax[0].twiny()
hp4,=ax_top.plot(CTD_dict[k_file]["DO"],CTD_dict[k_file]["depth"])
ax_top.set_xlabel("DO [mg/l]")
ax[0].legend([hp2,hp3,hp4],("Sorted Temp","Overturns","DO"))

ax[1].plot(np.abs(Thorpe_temp_dict[k_file]["thorpe_disp"]),CTD_dict[k_file]["depth"],'-k')
ax[1].plot(Thorpe_temp_dict[k_file]["Lt"],CTD_dict[k_file]["depth"],'-r')
#ax[1].set_xscale('log')
ax[1].set_xlabel("Displacements [m]")
ax[1].legend(("Thorpe displ","Thorpe scale"))
    
ax[0].invert_yaxis()
fig.suptitle("Prof {} raw (dist {} m)".format(profnames[k_file],disttrans[k_file]))

# Interpolated data

fig,ax=plt.subplots(1,2,sharey=True)
hp1,=ax[0].plot(CTDdata["temp"][:,k_time],CTDdata["depth"])
hp2,=ax[0].plot(CTDdata["temp_sorted"][:,k_time],CTDdata["depth"],'-k')
show_overturn=np.copy(CTDdata["temp"][:,k_time])
show_overturn[CTDdata["overturn_bool"][:,k_time]==0]=np.nan
hp3,=ax[0].plot(show_overturn,CTDdata["depth"],'-r')
ax[0].set_xlabel("Temp [°C]")
ax[0].set_ylabel("Depth [m]")
ax_top = ax[0].twiny()
hp4,=ax_top.plot(CTDdata["DO"][:,k_time],CTDdata["depth"])
ax_top.set_xlabel("DO [mg/l]")
ax[0].legend([hp2,hp3,hp4],("Sorted Temp","Overturns","DO"))

ax[1].plot(np.abs(CTDdata["thorpe_disp"][:,k_time]),CTDdata["depth"],'-k')
ax[1].plot(CTDdata["Lt"][:,k_time],CTDdata["depth"],'-r')
#ax[1].set_xscale('log')
ax[1].set_xlabel("Displacements [m]")
ax[1].legend(("Thorpe displ","Thorpe scale"))

fig.suptitle("Prof {} interpolated (dist {} m)".format(CTDdata["profnames"][k_time],CTDdata["disttrans"][k_time]))
    
ax[0].invert_yaxis()


#%% Export the data

CTD_attributes={"title":"CTD profiles along transects"}
CTD_dimensions={"time":{"dim_name": "time", "dim_size":None},"depth":{"dim_name": "depth", "dim_size":None}}
CTD_variables={"time": {"var_name": "time", "dim": ("time"),'unit': "seconds since 01.01.1970", 'longname': "time"},
               "depth":{"var_name": "depth", "dim": ("depth"),'unit': "m", 'longname': "depth (positive)"},
               "temp": {"var_name": "temp", "dim": ("depth","time"),'unit': "degC", 'longname': "water temperature"},
               "cond":{"var_name": "cond", "dim": ("depth","time"),'unit': "mS/cm", 'longname': "water conductivity"},
               "rho":{"var_name": "rho", "dim": ("depth","time"),'unit': "kg/m3", 'longname': "water density"},
               "DO":{"var_name": "DO", "dim": ("depth","time"),'unit': "mg/l", 'longname': "dissolved oxygen concentration"},
               "DOsat":{"var_name": "DOsat", "dim": ("depth","time"),'unit': "%", 'longname': "percentage of oxygen saturation"},
               "turb":{"var_name": "turb", "dim": ("depth","time"),'unit': "FTU", 'longname': "water turbidity"},
               "temp_sorted":{"var_name": "temp_sorted", "dim": ("depth","time"),'unit': "degC", 'longname': "sorted water temperature"},
               "thorpe_disp":{"var_name": "thorpe_disp", "dim": ("depth","time"),'unit': "m", 'longname': "Thorpe displacements"},
               "Lt":{"var_name": "Lt", "dim": ("depth","time"),'unit': "m", 'longname': "Thorpe length scale"},
               "overturn_bool":{"var_name": "overturn_bool", "dim": ("depth","time"),'unit': "-", 'longname': "Presence (1) or absence (0) of an overturn"},
               "xcoord":{"var_name": "xcoord", "dim": ("time"),'unit': "m", 'longname': "X coordinates (CH1903)"},
               "ycoord":{"var_name": "ycoord", "dim": ("time"),'unit': "m", 'longname': "Y coordinates (CH1903)"},
               "disttrans":{"var_name": "disttrans", "dim": ("time"),'unit': "m", 'longname': "Distance along the transect"},
               "profnames":{"var_name": "profnames", "dim": ("time"),'unit': "-", 'longname': "Profile name","var_type": "str"}}

export_netCDF(os.path.join(ctd_data_folder,"Level3","CTD_transect2_{}.nc".format(date_campaign)),CTD_attributes,CTD_dimensions,CTD_variables,CTDdata)

