# california-inlabru-forecasts

This folder contains all neccessary data and code to reproduce the forecasts and figures in "Pseudo-prospective testing of 5-year earthquake forecasts for California using inlabru", Bayliss, Naylor, Kamranzad and Main (
https://doi.org/10.5194/nhess-2021-403).  
UCERF3 data downloaded from https://pubs.usgs.gov/of/2013/1165/ in June 2019.  
The strain rate model used in this work is available at https://platform.openquake.org/maps/82 and was downloaded in April 2021.   
The exact versions are included in the 'Data' folder to ensure future reproducibility.
The Comcat_catalogues folder contains the catalogues used to test the forecasts. The included code downloads these directly, but they are saved here in case of any changes to the catalogue.

California_forecast_inlabru_testing.Rmd is an R-markdown file that runs all steps for model generation and then tests the models with the `pyCSEP` python package using the R package `reticulate`.  
inlabru_import_functions.R contains functions to generate the two types of forecast as well as process some of the input data.  
pycsep_testing.py runs only the pyCSEP testing using the inlabru forecasts, which are included in the 'Forecasts' folder.  

-----------------------------------------------------------------------------------

System information:
R version 4.1.0 (2021-05-18)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19041)
AMD Ryzen 5 3600 processor, 64GB RAM

R packages (forecast model assembly, fitting and simulation):  
reticulate_1.20  
RColorBrewer_1.1-2  
maptools_1.1-1      
dplyr_1.0.6       
raster_3.4-13
future.apply_1.8.1
future_1.22.1       
sf_1.0-2
INLA_21.06.11      
foreach_1.5.1     
Matrix_1.3-3      
fields_12.5       
viridis_0.6.1      
viridisLite_0.4.0   
spam_2.7-0         
dotCall64_1.0-1    
rgeos_0.5-5        
rgdal_1.5-23       
inlabru_2.3.1      
sp_1.4-5           
ggplot2_3.3.5   


Most packages can be downloaded directly from CRAN, with the exception of INLA which is using the current (as of December 2021) testing version,
downloaded with `install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)`  
The latest version of inlabru should also be installed using `remotes::install_github("inlabru-org/inlabru", ref="devel")` 

The forecast models take approx 2 hours each to run fully, including the simulation step and grid-based projections. 

------------------------------------------------------------------------------

Python packages (pycsep testing + plotting):  
Python 3.8.5  
pycsep 0.5.1  
numpy 1.20.3  
matplotlib  3.4.3  
cartopy 0.18.0  
pandas  1.3.4  
seaborn 0.11.1  

All python packages installed with conda, pyCSEP installed with `conda install --channel conda-forge pycsep`  
File pycsep_testing.py runs these steps only.

For further details on installation of the INLA package see https://www.r-inla.org/download-install  
For further details on installation of inlabru see https://github.com/inlabru-org/inlabru  
For further details on installation of pyCSEP, see https://docs.cseptesting.org/getting_started/installing.html  

All websites last accessed 16/12/2021
