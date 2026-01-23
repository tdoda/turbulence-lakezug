# Turbulence_LakeZug

## Introduction

This repository contains the data and scripts related to turbulence measurements in Lake Zug that are part of the Interbasin Exchanges in Stratified Lakes (IBEX) project.

Link to the remote repository: https://github.com/tdoda/turbulence-lakezug.git

This is a test.

## Sensors
### VMP profiler

<font color='red'>*Add details about the VMP here.*</font>

### CTD profiler

<font color='red'>*Add details about the CTD here.*</font>

## Geospatial information

<font color='red'>*Add information about the coordinates of the profiles.*</font>

## Temporal coverage

<font color='red'>*Add information about the field campaigns.*</font>

## Use the repository 

### Install the repository locally

To work with the project on your local machine, you need to have [git](https://git-scm.com/downloads) and [git-lfs](https://git-lfs.github.com/) installed in order to successfully clone the repository.

- Clone the repository to your local machine using the command:

 `git clone https://github.com/tdoda/turbulence-lakezug.git`
 
 Note that the repository will be copied to your current working directory.

- Use Python 3 and install the requirements with:

 `pip install -r requirements.txt`

 The python version can be checked by running the command `python --version`. In case python is not installed or only an older version of it, it is recommend to install python through the anaconda distribution which can be downloaded [here](https://www.anaconda.com/products/individual). 

### Process new VMP data

Here are the steps to process new data:
1. Create a new folder in `data\VMP` with the date of the campaign
2. Add the raw data files (`*.P`) in a `Level0` subfolder, with the `*.cfg` file and a `readme.txt` file
3. If default parameters are used for the turbulence analysis (quick method), go to the next step directly.\
Otherwise, modify the file `\functions\microstructure\load_parameters_Zug.m` by adding the campaign-specific parameters in the function `load_campaign_Zug()`.
4. In `analysis\1-Export_turbulence\analyze_profiles.m`, update the parameters of the campaign (`date_campaign`,`cfg_file` if not already specified in step 3, etc.)
5. Run `analysis\1-Export_turbulence\analyze_profiles.m`.


## Working with the project on the Renku platform (online)

The Renku project can be found here: https://renkulab.io/p/tomy.doda/turbulence-lakezug 

To run the scripts via Renku, click on the `Sessions` tab and start a new session.
This will start an interactive environment right in your browser.

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

- **Data collection**: Tomy Doda, Haydn Herrema, Damien Bouffard, Jemima Rama, Michael Plüss
- **Scripts and analysis**: Tomy Doda, Haydn Herrema, Oscar Sepúlveda Steiner, Bieito Fernández Castro, Sebastiano Piccolroaz, Hugo N. Ulloa, Damien Bouffard

## Contact
- [Tomy Doda](mailto:tomy.doda@unil.ch)
- [Haydn Herrema](mailto:haydn.herrema@eawag.ch)
- [Damien Bouffard](mailto:damien.bouffard@eawag.ch)
