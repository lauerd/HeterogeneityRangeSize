# Introduction

This repository contains the R code and the associated data files that were used to execute the entire workflow for the following project, which is now a published manuscript (cited in APA format):

Lauer, D. A., Shipley, B. R., & McGuire, J. L. (2023). Habitat and not topographic heterogeneity constrains the range sizes of African mammals. _Journal of Biogeography_, _50_(5), 846-857.

# Usage Instructions

By running the "code_for_project.R" file, a user can reproduce all of the analytical outputs and figures that are included in the published manuscript.

To run the code file, first download the file and the "data_files" folder. Ensure that both the file and folder are in the same directory, because the folder contains the data files that the code processes and analyses. Then, download and install the [RStudio Desktop](https://posit.co/download/rstudio-desktop/) integrated development environment for R. Check out the [RStudio User Guide](https://docs.posit.co/ide/user/) for more information on how to use it. Then, open RStudio, and use it to install the following R packages (note that many of these packages may have more updated versions at this point):

* ape (version 5.3)
* car (version 3.0.3)
* caret (version 6.0.84)
* cluster (version 2.1.0)
* DescTools (version 0.99.34)
* exactextractr (version 0.2.1)
* FSA (version 0.8.26)
* ggplot2 (version 3.3.5)
* geosphere (version 1.5.10)
* missForest (version 1.4)
* MuMIn (version 1.43.17)
* nlme (version 3.1.140)
* phytools (version 0.7.47)
* psych (version 1.9.12.31)
* raster (version 3.4.5)
* rgeos (version 0.5.5)
* sf (version 1.0.5)
* sp (version 1.4.5)
* spatialRF (version 1.1.2)
* tibble (version 3.1.6)

With the packages installed, open the code file in RStudio and run the code. Note that the code was written using R version 3.6.1.
