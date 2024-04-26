# The passing probabilities of earthquake gates and their effect in surface rupture length
A set of scripts to estimate the geometry of earthquake gates and passing probabilities as a function of geometry

<!-- ABOUT THE PROJECT -->
## About The Project
Propagating earthquakes must overcome geometrical complexity on fault networks to grow into large, surface rupturing events. We map step-overs, bends, gaps, splays, and strands of length scales ~100-500 meters from the surface ruptures of 31 strike-slip earthquakes, recording whether ruptures propagated past or halted at the feature. This set of scripts measures the geometry of each feature and classifies passing probabilities as a function of geometry for those features with distinct breached and unbreached populations. The scripts also enable characterizing and assessing the relationship between event likelihood (refer to manuscript) and surface rupture length given the distribution of geometrical complexity on a fault.

<!-- GETTING STARTED -->
## Data access

This repository contains the scripts required to reproduce the results in Rodriguez Padilla et al. 202X, and to measure the geometry and passing probabilities of different types of earthquake gates. The data required to run the scripts is provided in the following [repository](https://caltech-my.sharepoint.com/:f:/r/personal/alba_caltech_edu/Documents/GRL%20EQ%20gates%20materials/Data%20and%20code?csf=1&web=1&e=J9RgIm) and listed below:

- [ ] Data repository directories and content
    - [ ] primary_EQgate_shapefiles_v1
          Shapefiles of the earthquake gates mapped for each event. 
    - [ ] Regional_maps
          Shapefiles with the regional fault maps for each event.
    - [ ] FDHI_data
          Event information and displacement data for each event. Also available from the FDHI database appendix. 
    - [ ] event_kmz
          kmz files and primary rupture shapefiles for each event in the FDHI database. Kmz files sourced from the FDHI database appendix. Shapefiles generated using the kzm2shp script. 
          
Refer to the readme file in the repository for additional details. 

### Prerequisites for running the scripts

The subset of the scripts that measure the geometry of earthquake gates from shapefiles are in Matlab and require the Matlab Mapping Toolbox. Some of the scripts rely on functions downloable from Mathworks, and are provided as part of this repository in the source_code directory. The specific dependencies for each Matlab script to run are listed at the beginning of the corresponding script. The scripts for estimating passing probabilities and event likelihood are available as Python Jupyter Notebooks, with the functions stored in the .utils file in the directory.

<!-- CODE ROADMAP -->
### Measuring earthquake gate geometry, estimating passing probabilities, and estimating event likelihood

- [ ] To measure the geometry of earthquake gates in a shapefile
    - [ ] Run the "measure_EQgates.m" Matlab script (must run in directory containing the earthquake gate shapefiles)
    - [ ] This will output a csv file with the characterized gates
    - [ ] To measure the spacing between earthquake gate, run the "gatespacing.m" script. This script produces a pdf output fitting log-normal, Weibull, and exponential CDFs to the ECDF of the gate spacings (supplemental Figure S5).

- [ ] To estimate passing probabilities and event likelihood
    - [ ] Run the "analysis_EQgates_probabilities.ipynb" Jupyter Notebook. Requires the csv containing the gate geometries generated from the Matlab code in the previous step. This file ("aEQgate_geometries.csv") is also provided as part of this repository for users lacking access to Matlab or interested in accessing the geometry measurements directly without downloading the shapefiles and running the Matlab code. 
    - [ ] This script estimates passing probability as a function of geometry using logistic models.
    - [ ] This script also estimates the event likelihood based on the probabilities.
    - [ ] All figures in the main body of the manuscript except for Figure 1 and Figure 4b can be reproduced by running this code. All supplemental figures except figure S5 can be reproduced using this code too.

- [ ] To reproduce the rupture maps in the appendix with the earthquake gates plotted over them
    - [ ] Run the "map_maker.ipynb" script (requires the rupture maps and regional fault maps to be in shapefile format, provided in the data repository)

<!-- CONTACT -->
## Contact

Please report suggestions and issues:

Email: alba@caltech.edu

Project Link: [https://github.com/absrp/passing_probabilities_EQgates_strikeslip](https://github.com/absrp/passing_probabilities_EQgates_strikeslip)

Data link: [https://drive.google.com/drive/folders/1ecbHHmdSKvSZC_7zf8Pam6WqLQ-UTaZv?usp=share_link](https://drive.google.com/drive/folders/1ecbHHmdSKvSZC_7zf8Pam6WqLQ-UTaZv?usp=share_link)

Manuscript Link: * under review, stay tuned * 
[Pre-print link](https://essopenarchive.org/users/700353/articles/687281-the-influence-of-earthquake-gates-on-surface-rupture-length)


<p align="right">(<a href="#readme-top">back to top</a>)</p>
