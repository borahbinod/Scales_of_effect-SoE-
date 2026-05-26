# Scales_of_effect-SoE-
This repository contains derived data, scripts for the SoE project
# SoE BLISS Modeling

Short description of what this repo does.

## Overview
This repository contains code to run BLISS occupancy models for multiple species and BCR regions.

## Data
- Observations: `Data_derived/observations`
- Covariates: `Data_derived/covariates`

## Workflow
1. Prepare data
2. Run `SoE_BLISS()`
3. Save outputs per species–BCR

## Usage
```r
SoE_results <- SoE_BLISS(
  species = "Baltimore Oriole",
  bcr = "Atlantic Northern Forest"
)
