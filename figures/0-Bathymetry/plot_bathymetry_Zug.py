# Plot VMP or microCTD profiles from a specific campaign

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
from general_functions import read_bathy
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 12
plt.rcParams['font.family'] = 'Arial'  # or any installed font
plt.close('all')

#%% Parameters

datafolder=r'C:\Users\tdoda\Big_datasets_lakes\Zug'
dx_m=10 # [m]
dy_m=10 # [m]
lake_elev=414 # [m]

#%% Load HR bathymetry data

datafiles=os.listdir(datafolder)
depth_xr = read_bathy(datafolder,datafiles,lake_elev)

#%% Interpolation
dx=int(dx_m/(depth_xr.x.values[1]-depth_xr.x.values[0]))
dy=int(dy_m/(depth_xr.y.values[1]-depth_xr.y.values[0]))

x_interp=np.arange(np.min(depth_xr.x.values),np.max(depth_xr.x.values),dx)
y_interp=np.arange(np.min(depth_xr.y.values),np.max(depth_xr.y.values),dy)

depth_interp = depth_xr.interp(y=y_interp, x=x_interp)

#%% Lake boundary

# Mask of valid data
mask=~np.isnan(depth_interp.values)

# Find contours at level 0.5 (boundary between False and True)
contours = measure.find_contours(mask.astype(float), level=0.5)

# Select the **longest contour** (assume it's the outer lake)
outer_contour = max(contours, key=len)

# The contour coordinates are (row, col) indices, convert to X/Y
x_contour = depth_interp['x'].values[outer_contour[:, 1].astype(int)]
y_contour = depth_interp['y'].values[outer_contour[:, 0].astype(int)]

#%% Figure 

# Create a figure with specified size in centimeters
fig,ax = plt.subplots(1,1,figsize=(9 / 2.54, 10 / 2.54))  # Convert from cm to inches

hp=ax.pcolormesh(depth_interp.x, depth_interp.y,depth_interp,vmin=0,vmax=np.nanmax(depth_interp),cmap="cmo.deep")
#ax.contour(depth_interp.x, depth_interp.y,mask , levels=[0.5], colors='k', linewidths=1.5)
ax.plot(x_contour, y_contour, color='k', linewidth=1)
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

hp = ax.plot_surface(X[np.ix_(indzoom_y,indzoom_x)], Y[np.ix_(indzoom_y,indzoom_x)],-depth_xr.values[np.ix_(indzoom_y,indzoom_x)]*10, cmap="cmo.deep", antialiased=False,rstride=2, cstride=2) 
#ax.contour(X[np.ix_(indzoom_y,indzoom_x)], Y[np.ix_(indzoom_y,indzoom_x)],(-depth_xr.values[np.ix_(indzoom_y,indzoom_x)]+contour_offset)*10,levels=np.arange(-200*10, 0, 10*10),colors='k', linestyles='-', linewidths=1)
#ax.contourf(X[np.ix_(indzoom_y,indzoom_x)], Y[np.ix_(indzoom_y,indzoom_x)],-depth_xr.values[np.ix_(indzoom_y,indzoom_x)]*10, 50)
ax.set_aspect("equal")
cbar = fig.colorbar(hp, ax=ax)


#%% Figure 3D with shading

Z=depth_xr.values
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