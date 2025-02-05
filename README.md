# Turbulence_LakeZug

## Introduction

This repository contains the data and scripts related to turbulence measurements in Lake Zug that are part of the Interbasin Exchanges in Stratified Lakes (IBEX) project.

Link to the remote repository: https://gitlab.renkulab.io/tomy.doda/turbulence-lakezug.git

## Working with the project on the Renku platform (online)

The simplest way to start your project is right from the Renku
platform - just click on the `Sessions` tab and start a new session.
This will start an interactive environment right in your browser.

## Installation on a local repository 

To work with the project anywhere outside the Renku platform,
click the `Settings` tab where you will find the
git repo URLs - use `git` to clone the project on whichever machine you want.

**You need to have [git](https://git-scm.com/downloads) and [git-lfs](https://git-lfs.github.com/) installed in order to successfully clone the repository.**

- Clone the repository to your local machine using the command:

<span style="color:red">REPLACE WITH LINK TO YOUR REPOSITORY</span>

 `git clone 
 
 Note that the repository will be copied to your current working directory.

- Use Python 3 and install the requirements with:

 `pip install -r requirements.txt`

 The python version can be checked by running the command `python --version`. In case python is not installed or only an older version of it, it is recommend to install python through the anaconda distribution which can be downloaded [here](https://www.anaconda.com/products/individual). 

## Sensors
### VMP profiler

<font color='red'>*Add details about the VMP here.*</font>

### CTD profiler

<font color='red'>*Add details about the CTD here.*</font>

## Geospatial information

<font color='red'>*Add information about the coordinates of the profiles.*</font>

## Temporal coverage

<font color='red'>*Add information about the field campaigns.*</font>

## Structure of the repository

### Folder `data`
VMP and CTD data organized in levels. 

<font color='red'>*Add details about the data here.*</font>

<br />

### Folder `analysis`
Scripts performing data analysis.
#### Folder `1-Analyze_turbulence`
<font color='red'>*Add details about the turbulence analysis here.*</font>

#### Folder `1-Analyze_ctd`

- Script `read_ctd.py`: read and export the CTD data to netCDF files.

<br />

### Folder `functions`
Functions used by the analysis scripts.

<font color='red'>*Add more information here.*</font>

<br />

### Folder `figures`
Figures showing the results of the analyses.

<font color='red'>*Add more information here.*</font>

<br />

### Folder `notebooks`
Jupyter notebooks to visualize the data:
- CTD data: `plot_ctd.ipynb`
- VMP data: <font color='red'>*not available yet*</font>

## Future updates

- Export of the data as netcdf files (including level 1)
- Jupyter notebooks to visualize the data

## Collaborators

- **Data collection**: Tomy Doda, Damien Bouffard, Jemima Rama, Michael Plüss
- **Scripts and analysis**: Tomy Doda, Oscar Sepúlveda Steiner, Bieito Fernández Castro, Sebastiano Piccolroaz, Hugo N. Ulloa, Damien Bouffard

## Contact
- [Tomy Doda](mailto:tomy.doda@unil.ch)
- [Damien Bouffard](mailto:damien.bouffard@eawag.ch)
