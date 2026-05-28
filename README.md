# Scales_of_effect-SoE-

## Overview
This repository contains derived data, scripts for the SoE project

## Data
- Observations: `Data_derived/observations.Rds`
- Covariates: `Data_derived/covariates.Rds`
- Checklists: `Data_derived/checklists.Rds`

## Workflow
1. Prepare data
2. Run `SoE_BLISS()`
3. Save outputs per species–BCR

## Usage
```r
SoE_results <- SoE_BLISS(species = "Baltimore Oriole",bcr = "Atlantic Northern Forest",
                          niterations = 160000,nchains = 1,burnin = 10000, nthin = 5)
