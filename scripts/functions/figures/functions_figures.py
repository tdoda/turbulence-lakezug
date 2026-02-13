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
import matplotlib.pyplot as plt
import matplotlib.colors as mcol


def array_bands(x,v,dx):
    """Function array_bands

    Creates a grid with nan columns to plot it as bands with pcolormesh
    
    """

    xmat=np.copy(x)
    vmat=np.copy(v)

    for kd in range(len(x)-1):
        if (x[kd+1]-x[kd])>dx: # Add 2 NaN columns
            indmat=np.where(xmat==x[kd+1])[0][0]
            xmat=np.insert(xmat, indmat,x[kd+1]-dx/2)
            xmat=np.insert(xmat, indmat,x[kd]+dx/2)
            vmat=np.insert(vmat, indmat, np.full(vmat[:,0].shape,np.nan), axis=1)
            vmat=np.insert(vmat, indmat, np.full(vmat[:,0].shape,np.nan), axis=1)
    indmat=len(xmat)
    xmat=np.insert(xmat,indmat ,x[-1]+dx/2)
    vmat=np.insert(vmat,indmat, np.full(vmat[:,0].shape,np.nan), axis=1)
            
    return xmat,vmat

def get_cmap_discrete(cmap_name,n):
    """Function get_cmap_discrete

   Returns the discrete colors in a repeated cycle like in Matlab.
    
    """

    
    base_colors = plt.get_cmap(cmap_name).colors
    colprof = [base_colors[i % len(base_colors)] for i in range(n)]  
            
    return mcol.ListedColormap(colprof)

def compute_logmedian_dict(DATA,cell_size,varname,ind_file=np.nan,fieldname="BINNED"):
    """Function compute_logmedian_dict

    Computes the median of the log of a variable in depth cells of size cell_size for several profiles contained in a dictionary. The profiles can be either in a list (fieldname) or directly in the dictionary (fieldname=None).
    
    Inputs: 
        DATA (dictionary): dictionary containing the profiles data
        cell_size (float): size of the depth cells in meters
        varname (string): name of the variable to compute the log median for
        ind_file (list) [optional]: list of indices of the profiles to use in the dictionary. Default: np.nan (all profiles are used)
        fieldname (string) [optional]: name of the field in the dictionary containing the profiles. Default: None (the profiles are directly in the dictionary). The profiles can be either in a list (fieldname) or directly in the dictionary (fieldname=None). 

    Outputs:
        depth_allprof (array): all depths of the profiles combined
        var_allprof (array): all values of the variable combined
        depth_cells (array): depths of the cells
        logmed_cells (array): median of the log of the variable between cell depth (i-1) and cell depth i
   
    
     """

    if np.sum(np.isnan(ind_file))>0: # At least one NaN value
        ind_file=np.arange(len(DATA))

    # Combine the data
    var_allprof=np.array([])
    depth_allprof=np.array([])
    
    for k in range(len((ind_file))):
        if type(DATA[ind_file[k]][fieldname])==np.ndarray: # Several profiles
            for kprof in range(len(DATA[ind_file[k]][fieldname])): 
                datavar=DATA[ind_file[k]][fieldname][kprof][varname] 
                var_allprof=np.concatenate((var_allprof,datavar[~np.isnan(datavar)]))
                depth_allprof=np.concatenate((depth_allprof,DATA[ind_file[k]][fieldname][kprof]["depth"][~np.isnan(datavar)])) 
        else:
            datavar=DATA[ind_file[k]][fieldname][varname] 
            var_allprof=np.concatenate((var_allprof,datavar[~np.isnan(datavar)]))
            depth_allprof=np.concatenate((depth_allprof,DATA[ind_file[k]][fieldname]["depth"][~np.isnan(datavar)])) 
            
    #Compute the median of the log
    depth_start=np.arange(0,np.nanmax(depth_allprof)-cell_size ,cell_size)
    depth_end=np.arange(cell_size,np.nanmax(depth_allprof),cell_size)
    depth_cells=np.sort(np.concatenate((depth_start,depth_end))) 
    logmed_cells=np.full(depth_cells.shape,np.nan) 
    for k in range(len(depth_start)):
        indcell=np.where(depth_cells==depth_start[k])[0][-1]
        data_cell=var_allprof[(depth_allprof>=depth_start[k])&(depth_allprof<depth_end[k])]
        data_cell[data_cell<=0]=np.nan # Remove non positive values for log calculation
        if len(data_cell[~np.isnan(data_cell)])>0:
            logmed_cells[indcell:indcell+2]=np.nanmedian(np.log10(data_cell)) 

    return depth_allprof,var_allprof,depth_cells,logmed_cells

def compute_logmedian_list(depth_list,var_list,cell_size):
    """Function compute_logmedian_list
    
    var_list: list containing already the variable data
    cell_size in meters

    Returns median of the log of the variable between cell depth (i-1) and cell
    depth i
   
    
     """
 
    # Combine the data
    var_allprof=np.array([])
    depth_allprof=np.array([])
    

    for kprof in range(len(var_list)): 
        datavar=var_list[kprof] 
        var_allprof=np.concatenate((var_allprof,datavar[~np.isnan(datavar)]))
        depth_allprof=np.concatenate((depth_allprof,depth_list[kprof][~np.isnan(datavar)])) 
        
            
    #Compute the median of the log
    depth_start=np.arange(0,np.nanmax(depth_allprof)-cell_size ,cell_size)
    depth_end=np.arange(cell_size,np.nanmax(depth_allprof),cell_size)
    depth_cells=np.sort(np.concatenate((depth_start,depth_end))) 
    logmed_cells=np.full(depth_cells.shape,np.nan) 

    for k in range(len(depth_start)):
        indcell=np.where(depth_cells==depth_start[k])[0][-1]
        data_cell=var_allprof[(depth_allprof>=depth_start[k])&(depth_allprof<depth_end[k])]
        data_cell[data_cell<=0]=np.nan # Remove non positive values for log calculation
        if len(data_cell[~np.isnan(data_cell)])>0:
            logmed_cells[indcell:indcell+2]=np.nanmedian(np.log10(data_cell)) 

    return depth_allprof,var_allprof,depth_cells,logmed_cells

