
# Small area estimation of labour force indicators under unit-level multinomial mixed models

### Bugallo M., Esteban M.D., Hozba T., Morales D., Perez A.

### University Miguel Hernandez, Spain

## Objetive

This repository contains the R scripts of the *Application to real data* from the paper entlited "Small area estimation of labour force indicators under unit-level multinomial mixed models".

To cite this publication, please put the following paragraph: "Bugallo M., Esteban M.D., Hozba T., Morales D., Perez A. Small area estimation of labour force indicators under unit-level multinomial mixed models. Journal of the Royal Statistical Society: Series A."

## Installation and requirements

The R code available in this repository was tested on a 3.3 GHz Intel Xeon W  Mac computer with 32 GB RAM. For data manipulation, task processing and graphical results, the following R packages are required:

-   library(car)
-   library(nnet)
-   library(mefa)
-   library(lme4)
-   library(mvtnorm)
-   library(cubature)
-   library(nloptr)
-   library(fastDummies)
-   library(numDeriv)
-   library(KernSmooth)
-   library(maptools)
-   library(RColorBrewer)

In order to facilitate the manipulation of the code, we recommend the GUI, "Graphical User Interface", [R Studio](https://posit.co/downloads/). Once all the libraries are installed, simply download the project sources and run the scripts.

## Guide of use

To obtain a complete set of results of the *Application to real data* section of the paper, you need to follow the steps:

1.  *Obtain direct estimators*. Run the R script `Direct_Estimator2021T1.R` which returns the output result file `EPA_2021T1_Hajek.csv`. The data set required to run the script is `EPA_2021T1.csv`, which must be allocated in the same folder as the script.

2.  *Fit a unit-level multinomial mixed model*. Run the R script `Umlogit.R` which returns the output files `UNEMP.RATE.txt`, `INACTIVE.txt`, `EMPLOYED.txt` and `UNEMPLOYED.txt`, and some plots of the estimators of the multinomial model with random effects by province-sex. The data sets needed to run this script are `EPA_2021T1.csv`, `Censo_2021T1.csv` and `EPA_2021T1_Hajek.csv` (results of the previous step), which must be in the same folder as the script.

3.  *Obtain Bootstrap mse estimators*. Run the R script `Umlogit_BootRMSE.R` which returns the output files `mse_UNEMPRate_500.txt`, `mse_lab1_500.txt`, `mse_lab2_500.txt` and `mse_lab3_500.txt`. The data sets required to run this script are `EPA_2021T1.csv` and `Censo_2021T1.csv`, which must be in the same folder as the script. Please note that the bootstrap procedure uses a lot of computing time. Our script has been created to perform 500 bootstrap iterations. We recommend decreasing this number if you want to test the code.

4.  *Plot graphs of the results*. Using the maptools package, the `Maps2021T1.R` script provides several map plots of the above results. Please note that you need to have the files corresponding to the location you want to plot. In our case we have these files in the `spainmap' folder.


## Advertisement

The maptools package is now archived. If you experience problems with this package, maybe you can fix it with the following instance:

``` r
install.packages("maptools", repos = "https://packagemanager.posit.co/cran/2023-10-13")
```
