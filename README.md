# INLA_groupCV

This repository contains the [**R code**](https://github.com/spatialstatisticsupna/INLA_groupCV/tree/main/R/) to reproduce the analyses of the paper entitled *"Automatics cross-validation in structured models: Is it time to leave out leave-one-out?"* [(Adin et al., 2024)](https://doi.org/10.1016/j.spasta.2024.100843).

## Table of contents

-   [Section 3: Evaluating cross-validation measures in structured models](#section-3-evaluating-cross-validation-measures-in-structured-models)
-   [Section 4: Applications](#section-4-applications)
-   [Supplementary Material](#supplementary-material)
-   [References](#references)

# Section 3: Evaluating cross-validation measures in structured models

In the following examples, we are interested in analyzing the extrapolation performance of structured models, in which the unobserved data is often assumed less dependent on the observed data. For this end, we compare the behaviour of different models in terms of (i) Bayesian model selection criteria (such as DIC or WAIC), and (ii) predictive performance measures based on leave-one-out cross-validation (LOOCV) and the recently proposed automatic group construction procedure for leave-group-out cross-validation (LGOCV) for latent Gaussian models using the well-known INLA estimation technique [(Liu and Rue, 2022)](https://doi.org/10.48550/arXiv.2210.04482).

## Example 1: temporal models

In this example, we want to reproduce a real situation where the objective is to perform long-term forecasting using temporal models.

-   [**Example1_TemporalData.R**](./R/Section3/Example1_TemporalData.R)

    R code to fit with R-INLA two different temporally structured models (linear and first-order autoregressive model) for simulated data of imperfect observations of shortwave radition (X) and surface temperature (Y) for t=1,...,2000 data points.

-   [**Example1_TemporalData_Predictions.R**](./R/Section3/Example1_TemporalData_Predictions.R)

    R code to validate the predictive performance measures based on LOOCV/LGOCV for temporally structured models through a simulation study for long-term forecasting of surface temperature.

## Example 2: spatial models

The aim of this example is to validate the use of automatic groups construction of LGOCV under an extrapolation prediction task in disease mapping models, that is, when the objective is to predict disease risks or rates in unobserved areas using spatially structured models.

-   [**Example2_SpatialData.R**](./R/Section3/Example2_SpatialData.R)

    R code to fit with R-INLA different spatially structured models based on conditional autoregresive (CAR) priors and two-dimensional spatial P-splines for the estimation of dowry deaths relative risks in the 70 districts of Uttar Pradesh (the most populous state in India).

-   [**Example2_SpatialData_Predictions.R**](./R/Section3/Example2_SpatialData_Predictions.R)

    R code to validate the predictive performance measures based on LOOCV/LGOCV for spatial structured models through a simulation study where k-order neighbouring areas are removed to compare the estimated predictive distributions of dowry death counts in district against true values.

# Section 4: Applications

## Joint modelling of pancreatic cancer mortality and incidence data

In this section, we show how to use the automatic group construction procedure for LGOCV when the objective is to jointly modelling incidence and mortality data for a highly lethal disease such as pancreatic cancer using different spatio-temporal disease mapping models.

-   [**1_DiseaseMapping_FitModels.R**](./R/Section4/DiseaseMapping/1_DiseaseMapping_FitModels.R)

    R code to fit with R-INLA different univariate and multivariate (M-models and shared-component models) spatio-temporal models for pancreatic cancer mortality and incidence data in the 105 Clinical Commissioning Groups of England male population during the period 2001-2020. Official data provided by the [National Health Service in England](https://www.cancerdata.nhs.uk/incidence_and_mortality).

-   [**2_DiseaseMapping_Results.R**](./R/Section4/DiseaseMapping/2_DiseaseMapping_Results.R)

    R code to reproduce the results (figures and tables) from Section 4.1.

## Spatial Compositional Data: the case of *Arabidopsis thaliana*

In this section, we illustrate the utility of an automatic group construction procedure for LGOCV in a spatial Compositional Data problem: the geographical study of the plant *Arabidopsis thaliana* in the Iberian Peninsula.

We conducted our study using a dataset consisting of 301 accessions of the annual plant *Arabidopsis thaliana*. For each accession, the genetic cluster (GC) membership proportions for the four genetic clusters (GC1, GC2, GC3, GC4) are available. The entire dataset is available at <http://doi.org/10.5281/zenodo.2552025>.

-   [**Spatial_Coda.R**](./R/Section4/CoDa/Spatial_Coda.R)

    R code to fit the Logistic Normal Dirichlet models for spatial compositional data using R-INLA and reproduce the results (figures and tables) from Section 4.2.


## Space-time models for the United Kingdom wind speed data

Here, we illustrate how to use the LGOCV procedure to perform cross-validation on space-time models for wind speed data in the United Kingdom (UK).

-   [**1_get_wind_data.R**](./R/Section4/Spacetime/1_get_wind_data.R): R code to download, select and plot the wind speed data.

-   [**2_prepare_fit.R**](./R/Section4/Spacetime/2_prepare_fit.R): R code to build mesh, priors and likelihood model setup.
    
-   [**3_fit_models.R**](./R/Section4/Spacetime/3_fit_models.R): R code to run the four models and save selected results.

-   [**4_groups_plot.R**](./R/Section4/Spacetime/4_groups_plot.R): R code to make the groups figure.

-   [**5_results_models.R**](./R/Section4/Spacetime/5_results_models.R): R code to R code to extract some results.


# Supplementary Material

The purpose of the study presented in the Supplementary Material of the paper is to evaluate the ability of the automatic groups generated from different model candidate models in the context of model selection. All the scripts to reproduce the study are available at this [folder](./R/Supplement/).


# References
[Adin, A., Krainski, E., Lenzi, A., Liu, Z., Martínez-Minaya, J., and Rue, H. (2024). Automatic cross-validation in structured models: Is it time to leave out leave-one-out?. *Spatial Statistics*, 62, 100843.](https://doi.org/10.1016/j.spasta.2024.100843)

[Liu, Z., and Rue, H. (2022). Leave-group-out cross validation for latent Gaussian models. *arXiv preprint*.](https://doi.org/10.48550/arXiv.2210.04482)
