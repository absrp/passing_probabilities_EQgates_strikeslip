# The passing probabilities of zones of geometrical complexity and their effect in surface rupture length
A set of scripts to characterize different zones of geometrical complexity along surface ruptures and passing probabilities as a function of geometry

<!-- ABOUT THE PROJECT -->
## About The Project
Propagating earthquakes must overcome geometrical complexity on fault networks to grow into large, surface rupturing events. We map step-overs, bends, gaps, splays, and strands of length scales ~100-500 meters from the surface ruptures of 31 strike-slip earthquakes, recording whether ruptures propagated past the feature. We find that step-overs and bends can arrest rupture and develop a statistical model for passing probability as a function of geometry for each group. Step-overs wider than 1.2 km, single bends larger than 32°, and double bends larger than 38° are breached by rupture half of the time. ~20% of the ruptures terminate on straight segments. We examine how the distribution of geometrical complexity influences surface rupture length, inferring an exponential relationship between rupture length and event probability for a given fault. Our findings support that geometrical complexity limits the size of large events and help discriminate between different proposed models of rupture propagation.

<!-- GETTING STARTED -->
## Data access

This repository contains the scripts required to reproduce the results in Rodriguez Padilla et al. 202X, and to measure the geometry and passing probabilities of different types of earthquake gates. The data required to run the scripts is provided in the following [repository](https://caltech-my.sharepoint.com/:f:/r/personal/alba_caltech_edu/Documents/GRL%20EQ%20gates%20materials/Data%20and%20code?csf=1&web=1&e=J9RgIm) and listed below:

- [ ] Data repository directories and content
    - [ ] primary_EQgate_shapefiles_v1
          Shapefiles of the zones of geometrical complexity mapped for each event. 
    - [ ] Regional_maps
          Shapefiles with the regional fault maps for each event.
    - [ ] FDHI_data
          Event information and displacement data for each event. Also available from the FDHI database appendix. 
    - [ ] event_kmz
          kmz files and primary rupture shapefiles for each event in the FDHI database. Kmz files sourced from the FDHI database appendix. Shapefiles generated using the kzm2shp script. 
          
Refer to the readme file in the repository for additional details. 

### Prerequisites for running the scripts

The subset of the scripts that measure the geometry of different features from shapefiles are in Matlab and require the Matlab Mapping Toolbox. Some of the scripts rely on functions downloable from Mathworks, and are provided as part of this repository in the source_code directory. The specific dependencies for each Matlab script to run are listed at the beginning of the corresponding script. The scripts for estimating passing probabilities and event likelihood are available as Python Jupyter Notebooks, with the functions stored in the .utils file in the directory.

<!-- CODE ROADMAP -->
### Measuring fault geometry, estimating passing probabilities, and estimating event likelihood

- [ ] To measure the geometry offeatures in a shapefile
    - [ ] Run the "measure_EQgates.m" Matlab script (must run in directory containing the earthquake gate shapefiles)
    - [ ] This will output a csv file with the characterized features
    - [ ] To measure the spacing between zones of geometrical complexity along a rupture, run the "gatespacing.m" script. This script produces a pdf output fitting log-normal, Weibull, and exponential CDFs to the ECDF of the feature spacings (supplemental Figure S5).

- [ ] To estimate passing probabilities and event likelihood
    - [ ] Run the "analysis_EQgates_probabilities.ipynb" Jupyter Notebook. Requires the csv containing the gate geometries generated from the Matlab code in the previous step. This file ("aEQgate_geometries.csv") is also provided as part of this repository for users lacking access to Matlab or interested in accessing the geometry measurements directly without downloading the shapefiles and running the Matlab code. 
    - [ ] This script estimates passing probability as a function of geometry using logistic models.
    - [ ] This script also estimates the event likelihood based on the probabilities.
    - [ ] All figures in the main body of the manuscript except for Figure 1 and Figure 4b can be reproduced by running this code. All supplemental figures except figure S5 can be reproduced using this code too.

- [ ] To reproduce the rupture maps in the appendix with the earthquake gates plotted over them
    - [ ] Run the "map_maker.ipynb" script (requires the rupture maps and regional fault maps to be in shapefile format, provided in the data repository). The rupture maps in kzm format may be transformed in to shapefile format using the script kmz2shp.m

<!-- CONTACT -->
## Contact

Please report suggestions and issues:

Email: alba@caltech.edu

Project Link: [https://github.com/absrp/passing_probabilities_EQgates_strikeslip](https://github.com/absrp/passing_probabilities_EQgates_strikeslip)

Data link: [https://zenodo.org/records/11095762](https://zenodo.org/records/11095762)

Manuscript Link: * under review, stay tuned * 
[Pre-print link](https://essopenarchive.org/users/700353/articles/687281-the-influence-of-earthquake-gates-on-surface-rupture-length)


<p align="right">(<a href="#readme-top">back to top</a>)</p>
