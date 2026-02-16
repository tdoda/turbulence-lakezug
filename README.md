# Turbulence_LakeZug

## Introduction

This repository contains the data and scripts related to turbulence measurements in Lake Zug that are part of the Interbasin Exchanges in Stratified Lakes (IBEX) project.

Link to the remote repository: https://github.com/tdoda/turbulence-lakezug.git

## Sensors
### VMP profiler

<font color='red'>*Add details about the VMP here.*</font>

### CTD profiler

<font color='red'>*Add details about the CTD here.*</font>

## Geospatial information

<font color='red'>*Add information about the coordinates of the profiles.*</font>

## Temporal coverage

<font color='red'>*Add information about the field campaigns.*</font>

## Installation

### 1. Matlab installation

Matlab is required to run some of the scripts.

### 2. Python installation

Python 3 is required to run some of the scripts. Three types of installation are possible:
- Recommended option: download [Miniforge](https://github.com/conda-forge/miniforge). 
- User-friendly option: download the [Anaconda distribution](https://www.anaconda.com/products/individual).
- Classic option: download Python from the [official website](https://www.python.org/downloads/).

### 3. Repository installation

- If using GIT, clone the repository to your local machine using the command in Git Bash: 

    ``` 
    git clone https://github.com/tdoda/turbulence-lakezug.git
    ```
 
    Note that the repository will be copied to your current working directory.
- Without GIT, just download the entire ZIP folder from https://github.com/tdoda/turbulence-lakezug.git ("Code" > "Download ZIP") and extract it.

### 4. Python packages installation

<font color='red'>*Need to update environment and requirements files.*
</font>

1. Open the terminal (e.g., Anaconda Prompt), and move to the `turbulence_methods` repository.
2. Create a new environment *turbulence-methods* and install the packages as follows:
    - If using conda (Anaconda or Miniforge installation):
        ```
        conda env create -f environment.yml
        conda activate turbulence-methods 
        ```
        It is also possible to install the packages from `requirements.txt` with pip instead:
        ```
        conda create -n turbulence-methods python=3.11
        conda activate turbulence-methods
        pip install -r requirements.txt
        ```
    - If using mamba (Anaconda or Miniforge installation):
        ```
        mamba env create -f environment.yml
        mamba activate turbulence-methods 
        ```
    - If using pip (classic Python installation):
        ```
        python -m venv turbulence-methods       
        source turbulence-methods /bin/activate  # For Linux/macOS
        turbulence-methods\Scripts\activate     # For Windows
        pip install -r requirements.txt
        ```

## Working with the project online without any installation (Renku)

To run **Python** scripts and **Jupyter Notebooks** on the browser (without any installation), you can use the Renku
platform. Note that this approach **does not** work for Matlab scripts, that can only be ran [locally](#installation-on-a-local-repository) with Matlab installed.

Two options with Renku:
- Directly launch the created session by clicking here: [![launch - renku](https://renkulab.io/renku-badge.svg)](https://renkulab.io/p/tomy.doda/turbulence-lakezug/sessions/01KHBHMKZX1PCQAWQNA3S3GSHN/start)
- Access the Renku project [here](https://renkulab.io/p/tomy.doda/turbulence-lakezug) and create a new session in the `Sessions` tab. This will start an interactive environment right in your browser.

## Usage 

### Data processing

The raw (Level0) and processed (>Level1) data are stored in the `data` folder. The raw data is processed using the scripts in `scripts/data_processing` as explained below.

#### Turbulence data (microstructure)

The scripts and procedure are explained on the [turbulence-methods](https://github.com/tdoda/turbulence_methods.git) repository.

#### CTD

<font color='red'>*Add procedure.*</font>

#### Bathymetry

<font color='red'>*Add procedure.*</font>

### Data visualization
The processed data can be quickly visualized using the notebooks in `notebooks/visualization`.

### Data analysis
The processed data is further analyzed using the scripts in `scripts/data_analysis` and saved in the same folder.

### Data exploration
Data is explored (e.g., for a specific field campaign) using reports as notebooks in `notebooks/reports`.

### Figures
Specific figures from analyzed data are generated in `figures`.


## Structure of the repository

### Folder `data`
Data organized in processing levels. 

<font color='red'>*Add details about the data here.*</font>



### Folder `scripts`
Scripts processing the data performing further analyses.
#### `data_processing`
- `Turbulence`
<font color='red'>*Add details about the turbulence analysis here.*</font>

- `CTD`
    - Script `read_ctd.py`: read and export the CTD data to netCDF files.

#### `data_analysis`


#### `functions`
Functions used by the analysis scripts.

<font color='red'>*Add more information here.*</font>

### Folder `figures`
Figures showing the results of the analyses.

<font color='red'>*Add more information here.*</font>

<br />

### Folder `notebooks`
#### `visualization`
Jupyter notebooks to visualize the data:
- `plot_ctd.ipynb`: CTD data
- `plot_turbulence.ipynb`: turbulence data (including spectra)

#### `reports`
Jupyter notebooks exploring the datasets.

## Future updates

- ...

## Collaborators

- **Data collection**: Tomy Doda, Haydn Herrema, Damien Bouffard, Jemima Rama, Michael Plüss
- **Scripts and analysis**: Tomy Doda, Haydn Herrema, Oscar Sepúlveda Steiner, Bieito Fernández Castro, Sebastiano Piccolroaz, Hugo N. Ulloa, Damien Bouffard

## Contact
- [Tomy Doda](mailto:tomy.doda@unil.ch)

