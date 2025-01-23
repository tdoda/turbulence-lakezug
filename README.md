# Turbulence_LakeZug

## Introduction

This repository contains the data and scripts related to turbulence measurements in Lake Zug that are part of the Interbasin Exchanges in Stratified Lakes (IBEX) project.

Link to the remote repository: https://gitlab.renkulab.io/tomy.doda/turbulence-lakezug.git

## Working with the project

The simplest way to start your project is right from the Renku
platform - just click on the `Sessions` tab and start a new session.
This will start an interactive environment right in your browser.

To work with the project anywhere outside the Renku platform,
click the `Settings` tab where you will find the
git repo URLs - use `git` to clone the project on whichever machine you want.

Initially we install a very minimal set of packages to keep the images small.
However, you can add python and conda packages in `requirements.txt` and
`environment.yml` to your heart's content. If you need more fine-grained
control over your environment, please see [the documentation](https://renku.readthedocs.io/en/stable/topic-guides/customizing-sessions.html).

## Sensors
### VMP profiler

<font color='red'>*Add details about the VMP here...*</font>

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
