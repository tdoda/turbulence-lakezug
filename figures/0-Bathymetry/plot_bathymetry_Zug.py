# Plot bathymetry Lake Zug

# -*- coding: utf-8 -*-
import os
import sys
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import cmocean
from skimage import measure
from mpl_toolkits.mplot3d import Axes3D  # needed for 3D plotting
from matplotlib.colors import LightSource
sys.path.append(os.path.join(os.path.dirname(__file__), r'..\..\functions\general_functions'))
from general_functions import read_bathy, read_netCDF_xr
from pathlib import Path
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 12
plt.rcParams['font.family'] = 'Arial'  # or any installed font
plt.close('all')

# Set working directory to the script's location
os.chdir(Path(__file__).resolve().parent)

#%% Load interpolated bathymetry data from netCDF
depth_interp = read_netCDF_xr(os.path.join('..','..','data','Bathymetry','Level2','bathymetry_Zug_10m.nc'))

#%% Load HR data:
depth_xr = read_netCDF_xr(os.path.join('..','..','data','Bathymetry','Level2','bathymetry_Zug_1m.nc'))

#%% Figure 

# Create a figure with specified size in centimeters
fig,ax = plt.subplots(1,1,figsize=(9 / 2.54, 10 / 2.54))  # Convert from cm to inches

hp=ax.pcolormesh(depth_interp.x, depth_interp.y,depth_interp.depth,vmin=0,vmax=np.nanmax(depth_interp.depth),cmap="cmo.deep")
ax.plot(depth_interp.x_contour, depth_interp.y_contour, color='k', linewidth=1)
ax.set_aspect("equal")
cb = fig.colorbar(hp, ax=ax)
cb.set_label('Depth (m)')

#%% Figure 3D
xzoom_lim=[2679500,2680600]
yzoom_lim=[1218400,1218800]
indzoom_x=np.where((depth_xr.x.values>xzoom_lim[0])&(depth_xr.x.values<xzoom_lim[1]))[0]
indzoom_y=np.where((depth_xr.y.values>yzoom_lim[0])&(depth_xr.y.values<yzoom_lim[1]))[0]
contour_offset=5 # [m]

# Create a figure with specified size in centimeters
fig = plt.figure(figsize=(18 / 2.54, 10 / 2.54))
ax = fig.add_subplot(111, projection='3d')
X, Y = np.meshgrid(depth_xr.x, depth_xr.y)

hp = ax.plot_surface(X[np.ix_(indzoom_y,indzoom_x)], Y[np.ix_(indzoom_y,indzoom_x)],-depth_xr.depth.values[np.ix_(indzoom_y,indzoom_x)]*10, cmap="cmo.deep", antialiased=False,rstride=2, cstride=2) 
#ax.contour(X[np.ix_(indzoom_y,indzoom_x)], Y[np.ix_(indzoom_y,indzoom_x)],(-depth_xr.values[np.ix_(indzoom_y,indzoom_x)]+contour_offset)*10,levels=np.arange(-200*10, 0, 10*10),colors='k', linestyles='-', linewidths=1)
#ax.contourf(X[np.ix_(indzoom_y,indzoom_x)], Y[np.ix_(indzoom_y,indzoom_x)],-depth_xr.values[np.ix_(indzoom_y,indzoom_x)]*10, 50)
ax.set_aspect("equal")
cbar = fig.colorbar(hp, ax=ax)


#%% Figure 3D with shading

Z=depth_xr.depth.values
Z[np.isnan(Z)]=0

# Create a figure with specified size in centimeters
fig = plt.figure(figsize=(18 / 2.54, 10 / 2.54))
ax = fig.add_subplot(111, projection='3d')
X, Y = np.meshgrid(depth_xr.x, depth_xr.y)

# Create a LightSource object
ls = LightSource(azdeg=270, altdeg=45)  # light direction
# Shade the data
rgb = ls.shade(-Z[np.ix_(indzoom_y,indzoom_x)]*10, cmap=cmocean.cm.deep_r, vert_exag=2, blend_mode='soft')




hp = ax.plot_surface(X[np.ix_(indzoom_y,indzoom_x)], Y[np.ix_(indzoom_y,indzoom_x)],-Z[np.ix_(indzoom_y,indzoom_x)]*10, facecolors=rgb, antialiased=False,rstride=2, cstride=2,shade=False,linewidth=0)  

# ax.set_aspect("equal")
#cbar = fig.colorbar(hp, ax=ax)

#%% Display figures
plt.show()

