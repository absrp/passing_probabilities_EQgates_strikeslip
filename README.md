# The passing probabilities of earthquake gates and their effect in surface rupture length
A set of scripts to estimate the geometry of earthquake gates and passing probabilities as a function of geometry

<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->
<a name="readme-top"></a>
<!--
*** Thanks for checking out the Best-README-Template. If you have a suggestion
*** that would make this better, please fork the repo and create a pull request
*** or simply open an issue with the tag "enhancement".
*** Don't forget to give the project a star!
*** Thanks again! Now go create something AMAZING! :D
-->


<!-- ABOUT THE PROJECT -->
## About The Project
Earthquake magnitude is controlled by the rupture area of the fault network hosting the event. For surface-rupturing large strike-slip earthquakes (~MW6+), ruptures must overcome zones of geometrical complexity along fault networks. These zones, or earthquake gates, act as barriers to rupture propagation. We map step-overs, bends, gaps, splays, and strands from the surface ruptures of 31 strike-slip earthquakes, classifying each population into breached and unbreached groups. We develop a statistical model for passing probability as a function of geometry for each group. Step-overs, and single bends are more predictable earthquake gates than double bends and gaps, and ~20% of ruptures terminate on straight segments. Based on our modeled probabilities, we estimate event likelihood as the joint passing probabilities of breached gates and straight segments along a rupture. Event likelihood decreases inversely with  rupture length squared. Our findings support a barrier model as a factor in limiting large earthquake size.

<!-- GETTING STARTED -->
## Data access

This repository contains the scripts required to reproduce the results in Rodriguez Padilla et al. 202X, and to measure the geometry and passing probabilities of different types of earthquake gates. The data required to run the scripts is provided in the following [repository](https://drive.google.com/drive/folders/1ecbHHmdSKvSZC_7zf8Pam6WqLQ-UTaZv?usp=share_link) and listed below:

- [ ] Data repository directories and content
    - [ ] primary_EQgate_shapefiles_v1
          Shapefiles of the earthquake gates mapped for each event.
    - [ ] Regional_maps
          Shapefiles with the regional fault maps for each event.
    - [ ] FDHI_data
          Event information and displacement data for each event. Also available from the FDHI database appendix. 
    - [ ] event_kmz
          kmz files and primary rupture shapefiles for each event in the FDHI database. Kmz files sourced from the FDHI database appendix.
          
Refer to the readme file in the repository for additional details. 

### Prerequisites for running the scripts

The subset of the scripts that measure the geometry of earthquake gates from shapefiles are in Matlab and require the Matlab Mapping Toolbox. Some of the scripts rely on functions downloable from Mathworks, and are provided as part of this repository in the source_code directory. The specific dependencies for each Matlab script to run are listed at the beginning of the corresponding script. The scripts for estimating passing probabilities and event likelihood are available as Python Jupyter Notebooks, with the functions stored in the .utils file in the directory.

<!-- CODE ROADMAP -->
### Measuring earthquake gate geometry, estimating passing probabilities, and estimating event likelihood

- [ ] To measure the geometry of earthquake gates in a shapefile
    - [ ] Run the "measure_EQgates.m" Matlab script (must run in directory containing shapefiles)
    - [ ] This will output a csv file with the characterized gates
    - [ ] To measure the spacing between earthquake gate, run the "gatespacing.m" script. This script produces a pdf output fitting log-normal and exponential CDFs to the ECDF of the gate spacings.


- [ ] To estimate passing probabilities and event likelihood
    - [ ] Run the "analysis_EQgates_probabilities.ipynb" Jupyter Notebook. Requires the csv containing the gate geometries generated from the Matlab code in the previous step. This file is also provided as part of this repository for users lacking access to Matlab or interested in accessing the geometry measurements directly without downloading the shapefiles and running the Matlab code.
    - [ ] This script estimates passing probability as a function of geometry using logistic models.
    - [ ] This script also estimates the event likelihood based on the probabilities.
    - [ ] All figures in the main body of the manuscript except for Figure 1 can be reproduced by running this code. 

- [ ] To reproduce the rupture maps in the appendix with the earthquake gates plotted over them
    - [ ] Run the "map_maker.ipynb" script (requires the rupture maps to be in shapefile format, provided in the data repository)

<!-- CONTACT -->
## Contact

Please report suggestions and issues:

Email: alba@caltech.edu

Project Link: [https://github.com/absrp/passing_probabilities_EQgates](https://github.com/absrp/passing_probabilities_EQgates)

<p align="right">(<a href="#readme-top">back to top</a>)</p>

Manuscript Link: * under review, stay tuned *



