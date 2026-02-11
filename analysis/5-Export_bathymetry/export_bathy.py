# Plot bathymetry Lake Zug

# -*- coding: utf-8 -*-
import os
import sys
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from skimage import measure
sys.path.append(os.path.join(os.path.dirname(__file__), r'..\..\functions\general_functions'))
from general_functions import read_bathy
from pathlib import Path
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 12
plt.rcParams['font.family'] = 'Arial'  # or any installed font
plt.close('all')

# Set working directory to the script's location
os.chdir(Path(__file__).resolve().parent)

#%% Parameters

datafolder=r'C:\Users\tdoda\Big_datasets_lakes\Zug'
dx_m=10 # [m]
dy_m=10 # [m]
lake_elev=414 # [m]

export_txt=True # To export the interpolated bathymetry to a txt file (tab-separated, with x/y coordinates as row/column names)

#%% Load HR bathymetry data

datafiles = [f for f in os.listdir(datafolder) if f.endswith(".asc")]
depth_xr = read_bathy(datafolder,datafiles,lake_elev)

print("Data loaded, interpolate it")
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

#%% Export to netCDF

print("Saving data") 

# Convert to Dataset
ds = depth_interp.to_dataset(name="depth")

# Add x and y as variables
ds["x"] = ("x", ds.x.values)
ds["y"] = ("y", ds.y.values)

# Add contour coordinates as variables
ds = ds.assign_coords(
    contour_points=np.arange(len(x_contour))
)
ds["x_contour"] = ("contour_points", x_contour)
ds["y_contour"] = ("contour_points", y_contour)

# Add attributes
ds["depth"].attrs = {
    "long_name": "Lake depth",
    "units": "m",
    "positive": "down"
}

ds["x"].attrs = {
    "long_name": "East coordinates (CH1903+/LV95)",
    "units": "m"
}

ds["y"].attrs = {
    "long_name": "North coordinates (CH1903+/LV95)",
    "units": "m"
}

ds["x_contour"].attrs = {
    "long_name": "Lake boundary east coordinates (CH1903+/LV95)",
    "units": "m"
}

ds["y_contour"].attrs = {
    "long_name": "Lake boundary north coordinates (CH1903+/LV95)",
    "units": "m"
}

# Save in the same folder as the script:
# ds.to_netcdf('bathymetry_Zug_'+str(dx_m)+'m.nc', mode="w")

# Save in the same folder as the original data (if permission error):
ds.to_netcdf(os.path.join(datafolder, 'bathymetry_Zug_'+str(dx_m)+'m.nc'), mode="w")

print("Data saved!")

#%% Export to txt file

if export_txt:
    # Convert to DataFrame (rows=y, columns=x)
    df = pd.DataFrame(
        depth_interp.values,
        index=depth_interp.y.values,
        columns=depth_interp.x.values
    )

    # Optional: give index/column names
    df.index.name = "y"
    df.columns.name = "x"

    # Export to txt
    df.to_csv(os.path.join(datafolder, "bathymetry_Zug_"+str(dx_m)+"m.txt"), sep="\t", float_format="%.3f",na_rep="NaN")
