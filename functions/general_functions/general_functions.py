# -*- coding: utf-8 -*-
import pandas as pd
import os 
import sys
import numpy as np
from datetime import datetime, timezone, timedelta, UTC
import netCDF4 as nc
import xarray as xr
import shutil
import copy


def read_Vemco(filename):
    """ Function read_Vemco
    
    Read temperature data from Vemco loggers.

    Inputs:
    ----------
    filename (string): path to the file (*.csv)

    Outputs:
    -------
    tnum (numpy array (n,) of floats): timestamps [number of seconds since 01.01.1970]
    tdate (numpy array (n,) of datetime) : datetime values
    temp (numpy array (n,) of floats): temperature values [°C]

    """
    #datatemp=pd.read_csv(filename,skiprows=7,sep=",",parse_dates=[["Date(yyyy-mm-dd)","Time(hh:mm:ss)"]],encoding='latin1')
    datatemp=pd.read_csv(filename,skiprows=7,sep=",",dtype={'Date(yyyy-mm-dd)': str, 'Time(hh:mm:ss)': str},encoding='latin1')
    datatemp['datetime'] = pd.to_datetime(datatemp.pop('Date(yyyy-mm-dd)')+' '+datatemp.pop('Time(hh:mm:ss)'), format="%Y-%m-%d %H:%M:%S")
    datatemp=datatemp[[list(datatemp.columns)[i] for i in [1,0]]]
    tdate=datatemp['datetime'].to_numpy().astype('datetime64[s]').astype(datetime)
    tnum=datatemp.iloc[:,0].to_numpy().astype('datetime64[s]').astype(np.int64)
    temp=datatemp.iloc[:,1].to_numpy()
    
    return tnum, tdate, temp

def read_Tinytag(filename):
    """ Function read_Tinytag
    
    Read temperature data from Tinytag loggers.

    Inputs:
    ----------
    filename (string): path to the file (*.xlsx)

    Outputs:
    -------
    tnum (numpy array (n,) of floats): timestamps [number of seconds since 01.01.1970]
    tdate (numpy array (n,) of datetime) : datetime values
    temp (numpy array (n,) of floats): temperature values [°C]

    """
    datatemp=pd.read_excel(filename,skiprows=4,usecols=[1,2],names=["Date","Temp"])
    tdate=datatemp['Date'].to_numpy().astype('datetime64[s]').astype(datetime)
    tnum=datatemp['Date'].to_numpy().astype('datetime64[s]').astype(np.int64)
    temp=datatemp['Temp'].to_numpy()
    
    return tnum, tdate, temp
    
    
def calculate_hML(z,T,grad_tresh=0.05):
    """ Function calculate_hML
    
    Calculates mixed layer depth and thermocline depth from a temperature profile,
    using the T gradient method.
   
    Inputs:
    ----------
    z (numpy array (n,) of floats): z values (negative values, increasing upward) [m]
    T (numpy array (n,) of floats): temperature values [°C]
    grad_tresh (float): threshold for the temperature gradient [C°/m]
    
    Outputs:
    -------
    hML (float): mixed layer depth [m]
    htherm (float): thermocline depth [m]
    
    """
    
    if len(z)!=len(T):
        raise Exception('z et T do not have the same dimension')

    if len([z[z<0]])==0:
        raise Exception('z values must be negative')

    indsort=np.argsort(z)[-1::-1] # From surface to bottom
    z=z[indsort]
    T=T[indsort]
    
    gradT=np.concatenate((np.array([np.nan]),np.diff(T)/np.diff(z)));
    zgradT=np.concatenate((np.array([np.nan]),z[:-1]+np.diff(z)/2));
    indmaxgrad=np.nanargmax(gradT)
    maxgrad=gradT[indmaxgrad]
    htherm=-zgradT[indmaxgrad]

    if maxgrad<grad_tresh or indmaxgrad==len(gradT): # Fully mixed
        hML=np.nanmax(-z); # Take maximum depth
    else:
        indz_all=np.where(gradT[:indmaxgrad+1]<grad_tresh)[0]
        if not list(indz_all): # Entirely stratified
            hML=np.nan;
        else:
            indz=indz_all[-1]
            dz=zgradT[indz]-zgradT[indz+1]; # Positive
            dgrad=gradT[indz+1]-gradT[indz]; # Positive
            hML=-(zgradT[indz]-(grad_tresh-gradT[indz])/dgrad*dz); # Interpolation
  
    
    return hML,htherm


def calculate_N2(z,rho,window_smooth=1,rho0=np.nan,g=9.81):
    """ Function calculate_N2
    
    Calculates the squared buoyancy frequency profile from a density profile.
   
    Inputs:
    ----------
    z (numpy array (n,) of floats): z values (negative values, increasing upward) [m]
    rho (numpy array (n,) of floats): water density values [kg.m-3]
    window_smooth (int): size of the smoothing window
    rho0 (float): reference density [kg.m-3]. If nan, use the average density from the profile.
    g (float): gravitational acceleration [m.s-2]
    
    Outputs:
    -------
    zN2 (numpy array (n,) of floats): z values for which N2 is computed [m]
    N2 (numpy array (n,) of floats): squared buoyancy frequency as a function of depth [s-2]
    
    """
    N2=np.full(z.shape,np.nan)
    
    if np.sum(z<0)==0:
        raise Exception("Depth values must be negative")
    
    if window_smooth<=0:
        raise Exception("The window size must be >0")
        
    ind_start=int(np.floor(window_smooth/2))
    if window_smooth%2==0: # Even window size
        ind_end=-ind_start
        zN2=z
    else: # Odd window size
        ind_end=-ind_start-1
        zN2=np.concatenate((np.mean(np.array([z[1:],z[:-1]]),axis=0),np.array([np.nan])))
        
    if np.isnan(rho0):
        rho0=np.nanmean(rho)

    N2[ind_start:ind_end]=-g/rho0*(rho[window_smooth:]-rho[:-window_smooth])/(z[window_smooth:]-z[:-window_smooth])
    
    return zN2, N2

def calculate_St(z,rho,hypso_z,hypso_A,ztop=np.nan,zbot=np.nan,g=9.81):
    """ Function calculate_St
    
    Calculates the Schmidt stability [J.m-2] from a density profile.
   
    Inputs:
    ----------
    z (numpy array (n,) of floats): z values (negative values, increasing upward) [m]
    rho (numpy array (n,) of floats): water density values [kg.m-3]
    hypso_z (numpy array (n,) of floats): z values (negative values, increasing upward) for hypsometry data [m]
    hypso_A (numpy array (n,) of floats): lake area as a function of depth [m2]
    ztop (float): upper depth of the layer where St is computed [m]. If nan, use the top value of z.
    zbot (float): lower depth of the layer where St is computed [m]. If nan, use the bottom value of z.
    g (float): gravitational acceleration [m.s-2]
    
    Outputs:
    ------- 
    St_cv (float): area-specific Schmidt stability based on center of volume [J.m-2] (e.g., Read et al., 2011; Sahoo et al., 2016)
    St_cg_rel (float): area-specific Schmidt stability based on center of gravity and avg density [J.m-2] (e.g., Sahoo et al., 2016)
    St_cv_rel (float): area-specific Schmidt stability based on center of volume and avg density [J.m-2] (e.g., Wetzel et al., 2024)
    Ps (float): potential energy of stratification with respect to avg density [J.m-2] (e.g., Simpson et al., 1978; Ward et al., 1990; Kastev et al., 2010)
    z (numpy array (m,) of floats): z values used to compute St (negative values, increasing upward) [m]
    rho (numpy array (m,) of floats): water density values used to compute St [kg.m-3]
    Az (numpy array (m,) of floats): depth-varying lake area used to compute St [m2]
    """

    if ~np.isnan(zbot):
        rho=rho[z>=zbot]
        z=z[z>=zbot]
    
    if ~np.isnan(ztop):
        rho=rho[z<=ztop]
        z=z[z<=ztop]
        # Note: we can also set make z relative to ztop so that z starts at 0m but this will not affect St computed with (z-zv) or (z-zg)
    
    # Order the data from bottom to top 
    indz=np.argsort(z)
    z=z[indz]
    rho=rho[indz]
    ind_hypso=np.argsort(hypso_z)
    hypso_z=hypso_z[ind_hypso]
    hypso_A=hypso_A[ind_hypso]


    Az=np.interp(z,hypso_z,hypso_A,left=np.nan,right=np.nan)
    A0=np.max(hypso_A)
    
    ind0=np.where(np.logical_and(np.logical_and(~np.isnan(z), ~np.isnan(rho)),~np.isnan(Az)))[0][0]
    indf=np.where(np.logical_and(np.logical_and(~np.isnan(z), ~np.isnan(rho)),~np.isnan(Az)))[0][-1]
    if ind0>0 or indf<(len(z)-1):
        print("Only a part of the water column is used to compute St: {:.2f}-{:.2f} m".format(-z[indf],-z[ind0]))
    z=z[ind0:indf+1]
    rho=rho[ind0:indf+1]
    Az=Az[ind0:indf+1]
    
    # rhomean=np.nanmean(rho)
    rhomean=np.trapz(Az*rho,z)/np.trapz(Az,z)

    zv=np.trapz(z*Az,z)/np.trapz(Az,z)
    zg=np.trapz(z*Az*rho,z)/np.trapz(Az*rho,z)
    
    St_cv=g/A0*np.trapz((zv-z)*rho*Az,z) # [J/m2]
    St_cg_rel=g/A0*np.trapz((zg-z)*(rho-rhomean)*Az,z) # [J/m2]
    St_cv_rel=g/A0*np.trapz((zv-z)*(rho-rhomean)*Az,z) # [J/m2]
    Ps=g/A0*np.trapz(z*(rhomean-rho)*Az,z) # [J/m2]
    
    
    return St_cv, St_cg_rel, St_cv_rel, Ps, z, rho, Az  

def movmean(X,windowsize,axis=0):
    """Function movmean

    Computes the moving average of an array centered at the given index.

    Inputs:
    ----------
    X (numpy array (m,n) of floats): array to average
    windowsize (int): size of the averaging window
    axis (int): index of the axis along which the averaging is applied
        
    Outputs:
    ----------
    X_smooth (numpy array (m,n) of floats): smoothed array
    
    """

    if len(X.shape)==1:
        X=np.expand_dims(X,axis=1)
    if axis==1:
        X=X.transpose()
    X_smooth=np.full(X.shape,np.nan)
    for k in range(X.shape[1]): 
        df=pd.DataFrame({'val':X[:,k]})
        X_smooth[:,k]=df.rolling(windowsize,center=True).mean().values[:,0]
    if axis==1:
        X_smooth=X_smooth.transpose()
    return X_smooth

def grad_order(x,f,order=1):
    """ To fill

    """
    
    dfdx=np.full(x.shape,np.nan)
    dfdx[order:]=(f[order:]-f[:-order])/(x[order:]-x[:-order])

    
    return dfdx

def grad_order(x,f,order=1):
    """ To fill

    """
    
    dfdx=np.full(x.shape,np.nan)
    dfdx[order:]=(f[order:]-f[:-order])/(x[order:]-x[:-order])

    
    return dfdx

def export_netCDF(filename,general_attributes,dimensions,variables,data):
    """Function export_netCDF

    Export data to a netCDF file

    Inputs:
    ----------
    filename (string): netCDF filename with path included
    general_attributes (dictionary): file attributes with the format {"name_attribute1": "attribute1_value",…}
    dimensions (dictionary): dimensions with the format {'name_dimension1': {'dim_name': 'name_dimension1', 'dim_size': ...},…}
    variables (dictionary): variable names with the format {'name_variable1': {'var_name': 'name_variable1', 'dim': ('name_dim1',’name_dim2’,…),'unit': “name_units”, 'longname': 'long_name_variable1', 'var_type':'float' or 'str'},…}
    data (dictionary): data to export with the format {“name_variable1”:numpy_array,…}
    vartype 
    
        
    Outputs: None
    
    """

    nc_file = nc.Dataset(filename, mode='w', format='NETCDF4')
    #nc_file = create_netCDF(filename)

    for key in general_attributes:
        setattr(nc_file, key, general_attributes[key])
        
    for key, values in dimensions.items():
         nc_file.createDimension(values['dim_name'], values['dim_size'])


    for key, values in variables.items(): 
        if 'var_type' not in values.keys() or values["var_type"]=="float":
            var = nc_file.createVariable(values["var_name"], np.float64, values["dim"], fill_value=np.nan)
        elif values["var_type"]=="str":
            var = nc_file.createVariable(values["var_name"], str, values["dim"])
        else:
            raise Exception("Variable type is unknown")
        
        var.units = values["unit"]
        var.long_name = values["longname"]
        if key in data.keys():
            if 'var_type' not in values.keys() or values["var_type"]=="float":
                var[:] = data[key]
            else:
                for k in range(len(data[key])):
                    var[k]=data[key][k]
        else:
            raise Exception("Data is missing for {}".format(key))
    nc_file.close()
    print("Data exported to netCDF!")
    
def read_netCDF_xr(pathname,rootpath='C:/Users/tdoda'):
    """Function read_netCDF_xr

    Read a netCDF file as an xarray, working even if the path contains non-ASCII characters

    Inputs:
    ----------
    pathname (string): netCDF filename with path included
    rootpath (string): basic path without accent (typically C:/Users/username)
    
        
    Outputs:
    ----------
    data_xr (xarray dataset): netCDF data as an xarray
    
    """
    
    # Get the file name
    ind_slash=pathname.rfind("/")
    ind_backslash=pathname.rfind("\\")
    ind_filename=max(ind_slash,ind_backslash)
    if ind_filename>=0:
        filename=pathname[ind_filename+1:]
    else:
        filename=pathname
    
    
    # current_path=os.getcwd()
    # source_path=os.path.join(current_path, pathname)
    
    destination_path = os.path.join(rootpath, filename)  # Full path for the destination file

    # Copy the file 
    try:
        shutil.copy(pathname, destination_path)
    except Exception as e:
        print(f"Error copying file: {e}")
        
    # Open file
    data_xr=xr.open_dataset(destination_path)
    data_xr.close() # Close the file
    data_xr=data_xr.load() # Load data into memory to be independent of teh file
    
    # Delete file
    os.remove(destination_path)
    
    return data_xr

def create_netCDF(pathname,mode_name='w', format_name='NETCDF4', rootpath='C:/Users/tdoda'):
    """Function create_netCDF_xr

    Create a netCDF file even if the path contains non-ASCII characters.

    Inputs:
    ----------
    pathname (string): netCDF filename with path included
    mode_name
    format_name
    rootpath (string): basic path without accent (typically C:/Users/username)
    
        
    Outputs:
    ----------
    nc_file (netCDF object): netCDF dataset
    
    """
    
    # Get the file name
    ind_slash=pathname.rfind("/")
    ind_backslash=pathname.rfind("\\")
    ind_filename=max(ind_slash,ind_backslash)
    if ind_filename>=0:
        filename=pathname[ind_filename+1:]
    else:
        filename=pathname
    
    
    # current_path=os.getcwd()
    # source_path=os.path.join(current_path, pathname)
    
    destination_path = os.path.join(rootpath, filename)  # Full path for the destination file
    
    # Create netCDF file there
    nc_file = nc.Dataset(destination_path, mode=mode_name, format=format_name)
    
    return nc_file

def close_netCDF(nc_file,pathname,rootpath=r'C:/Users/tdoda'):
    """Function create_netCDF_xr

    Closed an open netCDF file created with create_netCDF() and save it in the path containing non-ASCII characters.

    Inputs:
    ----------
    nc_file (netCDF object): netCDF dataset
    pathname (string): netCDF filename with path included (where to save it)
    rootpath (string): basic path without accent (typically C:/Users/username)
    
        
    Outputs:
    ----------
    None
    """
    # Close the file
    nc_file.close()
    
    # Get the file name
    ind_slash=pathname.rfind("/")
    ind_backslash=pathname.rfind("\\")
    ind_filename=max(ind_slash,ind_backslash)
    if ind_filename>=0:
        filename=pathname[ind_filename+1:]
    else:
        filename=pathname
    
    current_path = os.path.join(rootpath, filename)  # Full path where the file has been created
    
    if not os.path.exists(current_path):
        raise Exception("File does not exist")
    
    # Copy the file 
    try:
        shutil.copy(current_path,pathname)
    except Exception as e:
        print(f"Error copying file: {e}")
        
    
    # Delete file
    os.remove(current_path)

def netCDF2dict(nc):
    """Function netCDF2dict

    Converts a netCDF object to a dictionary
    
    """
    
    nc_data=dict()
    nc_genatt=dict()
    nc_varatt=dict()
    
    for key,value in nc.variables.items():
        if value.dtype==np.float64:
            nc_data[key]=value[:].data
        else:
            nc_data[key]=value[:]
        nc_varatt[key]=dict()
        for att in value.ncattrs():
            nc_varatt[key][att]=getattr(value,att)

    for att in nc.ncattrs():
        nc_genatt[att]=getattr(nc,att)

                
    return nc_data, nc_genatt, nc_varatt

def dist_transect(xpt,ypt,xtrans=np.nan,ytrans=np.nan,monotonic_dir='y',method='projected'):
    """Function dist_transect

    Computes distance along transect
    
    method='projected' on transect or 'sum' of final distance
    
    """

    if (np.sum(~np.isnan(xtrans))==0) | (np.sum(~np.isnan(ytrans))==0): # Use default points
        xtrans=np.array([679815.,679803., 679803., 679663., 679415., 679505., 680030., 680150.,
               680223., 680409., 680814., 681350., 681909.])
        ytrans=np.array([225711.,224100., 222690., 221650., 220433., 219544., 218872., 218462.,
               218126., 217437., 216239., 214957., 213760.])
    
    dist_trans=np.concatenate((np.array([0]),np.sqrt(np.diff(xtrans)**2+np.diff(ytrans)**2)))
    
    if monotonic_dir=='y':
        transcoord=ytrans
        ptcoord=ypt
    else:  # x
        transcoord=xtrans
        ptcoord=xpt
    
    if transcoord[1]>transcoord[0]: # ascending order
        ind_point=np.where(transcoord<=ptcoord)[0]
    elif transcoord[1]<transcoord[0]: # descending order
        ind_point=np.where(transcoord>=ptcoord)[0]
    else:
        raise Exception("Transect coordinates should be monotonous in the selected direction")
        
    if np.size(ind_point)==0 or (len(ind_point)==len(transcoord) and transcoord[-1]!=ptcoord): 
        raise Exception("Point outside of the transect")
    if len(ind_point)<len(transcoord):
        if method=='projected':
            dist_pt=np.cumsum(dist_trans)[ind_point[-1]]+(ptcoord-transcoord[ind_point[-1]])/(transcoord[ind_point[-1]+1]-transcoord[ind_point[-1]])*dist_trans[ind_point[-1]+1]
        elif method=='sum':
            dist_pt=np.cumsum(dist_trans)[ind_point[-1]]+np.sqrt((xpt-xtrans[ind_point[-1]])**2+(ypt-ytrans[ind_point[-1]])**2)
        return dist_pt
    else:
        return np.sum(dist_trans)
    
    
def bathy_transect(x,y,D,xtrans=np.nan,ytrans=np.nan,monotonic_dir='y',dx=1):
    """Function bathy_transect

    Extract bathymetry along transect
    
    """
    if D.shape[0]!=len(y) or D.shape[1]!=len(x):
        raise Exception("Wrong dimensions")

    if (np.sum(~np.isnan(xtrans))==0) | (np.sum(~np.isnan(ytrans))==0): # Use default points
        xtrans=np.array([679815.,679803., 679803., 679663., 679415., 679505., 680030., 680150.,
               680223., 680409., 680814., 681350., 681909.])
        ytrans=np.array([225711.,224100., 222690., 221650., 220433., 219544., 218872., 218462.,
               218126., 217437., 216239., 214957., 213760.])
    
    if monotonic_dir=='y': 
        ybathy=np.arange(np.min(ytrans),np.max(ytrans),dx)
        ind_sort=np.argsort(ytrans)
        xbathy=np.interp(ybathy,ytrans[ind_sort],xtrans[ind_sort],left=np.nan,right=np.nan)
    else:  # x
        xbathy=np.arange(np.min(xtrans),np.max(xtrans),dx)
        ind_sort=np.argsort(xtrans)
        ybathy=np.interp(xbathy,xtrans[ind_sort],ytrans[ind_sort],left=np.nan,right=np.nan)
        
    dist_bathy=np.full(xbathy.shape,np.nan)
    depth_bathy=np.full(xbathy.shape,np.nan)
    for k in range(len(dist_bathy)):
        dist_bathy[k]=dist_transect(xbathy[k],ybathy[k],xtrans=xtrans,ytrans=ytrans,monotonic_dir=monotonic_dir,method="projected")
        try:
            depth_bathy[k]=D[np.where(y>=ybathy[k])[0][0],np.where(x>=xbathy[k])[0][0]]
        except:
            breakpoint()
    
    ind_sort_dist=np.argsort(dist_bathy)
    xbathy=xbathy[ind_sort_dist]
    ybathy=ybathy[ind_sort_dist]
    dist_bathy=dist_bathy[ind_sort_dist]
    depth_bathy=depth_bathy[ind_sort_dist]
    
    
    return xbathy,ybathy,dist_bathy,depth_bathy


def matstruct2dict(mat_struct):
    """Function matstruct2dict

    Converts Matlab structures into dictionary
    
    """
    key_names=list(dir(mat_struct))
    D=dict()
    for key in key_names:
        if key[0]!="_":
            D[key]=getattr(mat_struct,key)

    return D
    


