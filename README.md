# sprom

The R package *sprom* is the companion package for the paper Castillo-Mateo et al. "Spatio-temporal modeling for record-breaking temperature events in Spain". The package includes scripts and functions to reproduce all the results from the paper. 

## Installation

You can install version 0.0.1 from the provided folder
```s
if (!require("remotes")) install.packages("remotes")
remotes::install_local("reproducibility_materials.zip")
```

You can install the **development** version from
[GitHub](https://github.com/JorgeCastilloMateo/sprom)
```s
if (!require("remotes")) install.packages("remotes")
remotes::install_github("JorgeCastilloMateo/sprom")
```

## Workflow of the paper

* The data that contains basic information for the weather stations and that will be updated in subsequent scripts is in `data\stations.rda`
* The scripts to reproduce the results from the paper are in `inst\scripts\`:
  + `01_1_download_Tmax_data.R` downloads the daily maximum temperature data from ECA
  + `01_2_grid_dist.R` builds a 10 x 10 km grid for peninsular Spain and obtains the distance to the coast for each grid cell and weather station
  + `01_3_figure_1_map.R` reproduces the map in Figure 1 of the Manuscript
  + `01_4_simulations_missing_ties.R` includes simulations to check the impact of missing data on the results (see Section 2.1 of the Manuscript) and explores tied records (see Section 2.2 of the Manuscript)
  + `01_5_EDA.R` includes code to reproduce some figures for the exploratory data analysis (see Section 2.3 of the Manuscript)
  + `02_1_data.R` builds the data frame that will be used to fit the models
  + `02_2_KFCV.R` obtains the metrics using 10-fold cross-validation for each model (see Table 2 of the Manuscript and Figures 7 and 8 in Section 5.1 of the Supplementary material)
  + `03_1_fitting_full.R` implements model fitting for the full model, and assesses its convergence (see Section 4 of the Supplementary material)
  + `03_2_fitting_others.R` implements model fitting for all the other models
  + `03_3_DIC_checks.R` obtains DIC for all models and does the posterior predictive checks (adequacy) for the full model (see last column in Table 2 of the Manuscript and Table 4 and Figure 9 in Section 5.1 of the Supplementary material)
  + `04_1_results_parameters.R` uses the full model parameters to do posterior inference (see _Posterior distribution of the model parameters_ in Section 4.2 of the Manuscript and Section 5.2 of the Supplementary material)
  + `04_2_results_inference` obtains model-based tools using posterior predictive samples from the full model (see _Model-based analysis over peninsular Spain_ in Section 4.2 of the Manuscript and Section 5.3 of the Supplementary material)
   + `05_1_simulation_study` obtains the results of the illustrative simulation study that compares a binary model approach against a mean model approach (see Table 9 in Section 6 of the Supplementary material)
