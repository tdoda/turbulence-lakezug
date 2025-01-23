# -*- coding: utf-8 -*-
"""
Create KML file with GPS coordinates 
Created on Thu Jan 16 15:33:33 2025

@author: tdoda
"""
import os
import sys
import pandas as pd
import numpy as np
from pyproj import Proj, transform
import simplekml
import netCDF4
sys.path.append(os.path.join(os.path.dirname(__file__), r'..\..\functions\general_functions'))
from general_functions import read_netCDF_xr

#%% Extract the coordinates from netCDF files (already WGS84 format)
date_campaign='20250115'
folder = os.path.join("../../data/CTD/",date_campaign,"Level1")
files = [os.path.join(folder, f) for f in os.listdir(folder) if f.endswith(".nc")]
files.sort()

nameprof=[]
long=np.full((len(files),),np.nan)
lat=np.full((len(files),),np.nan)

for kf in range(len(files)):
   data_nc=read_netCDF_xr(files[kf]) # Read netCDF
   long[kf]=data_nc.attrs["longitude"]
   lat[kf]=data_nc.attrs["latitude"]
   nameprof.append(data_nc.attrs["Profile name"])

#%% Convert coordinates
# # Define the coordinate systems
# ch1903 = Proj(init='epsg:21781')  # CH1903 (Swiss Grid)
# wgs84 = Proj(init='epsg:4326')    # WGS84 (Lat/Lon)


# # Convert each point from CH1903 (East, North) to WGS84 (Lat, Lon)
# for index, row in df.iterrows():
#     # Convert from CH1903 (East, North) to WGS84 (Longitude, Latitude)
#     lon, lat = transform(ch1903, wgs84, row['Longitude_CH1903'], row['Latitude_CH1903'])
#     latitudes.append(lat)
#     longitudes.append(lon)

# # Add the converted coordinates as new columns to the DataFrame
# df['Latitude'] = latitudes
# df['Longitude'] = longitudes

#%% Create a KML object
kml = simplekml.Kml()

# Add each point to the KML file as a placemark
for krow in range(len(nameprof)):
    kml.newpoint(name=nameprof[krow], coords=[(long[krow], lat[krow])])

# Save the KML file
kml.save("coordinates.kml")

print('KML file saved')