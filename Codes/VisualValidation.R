# to create caterpillar plots to check them
library(MCMCvis)
### reading individual model output files
y<-list.files(path="Data_derived/SoE_results",pattern="\\.Rds$",full.names = T) # change path accordingly
rds_list <- lapply(y, readRDS)
for (i in 1:length(rds_list)){
SoE_results<- rds_list[[i]] # reading out species output 

cov.pars<- c("beta.0","beta.developed_low_intensity","beta.developed_low_intensity_sq","beta.developed_high_intensity","beta.developed_high_intensity_sq",
             "beta.barren_land","beta.barren_land_sq","beta.deciduous_forest","beta.deciduous_forest_sq","beta.evergreen_forest","beta.evergreen_forest_sq","beta.mixed_forest","beta.mixed_forest_sq",
             "beta.shrub_scrub","beta.shrub_scrub_sq","beta.grassland_herbaceous","beta.grassland_herbaceous_sq",
             "beta.pasture_hay","beta.pasture_hay_sq","beta.cultivated_crop","beta.cultivated_crop_sq",
             "beta.woody_wetland","beta.woody_wetland_sq","beta.emergent_wetland","beta.emergent_wetland_sq",
             "beta.ocean","beta.ocean_sq","beta.lakes","beta.lakes_sq","beta.rivers","beta.rivers_sq",
             "beta.elevation","beta.elevation_sq", "beta.aspect","beta.aspect_sq", "beta.slope","beta.slope_sq",
             "beta.precipitation","beta.precipitation_sq","beta.max_temperature","beta.max_temperature_sq", 
             "beta.effort_hours", "beta.number_observers",
             "beta.time_of_day")
scale.pars<- c("scale.1", "scale.2", "scale.3", "scale.4","scale.5", 
               "scale.6", "scale.7", "scale.8","scale.9", "scale.10", 
               "scale.11", "scale.12","scale.13", "scale.14", 
               "scale.15", "scale.16", "scale.17", 
               "scale.18", "scale.19", "scale.20")
message(paste0("Model outputs for species: ", SoE_results$scales_results$species[1]))
message(paste0("presences: ", length(which(SoE_results$data$z==1))," absences: ",length(which(SoE_results$data$z==0)) ))
SoE_results$start_time; SoE_results$finish_time
SoE_results$niterations; SoE_results$nthin
MCMCvis::MCMCtrace(SoE_results$samples, 
                   params = c(cov.pars,scale.pars), # checking three coefficients at a time
                   ISB = FALSE, iter= 30000,
                   exact = TRUE,
                   pdf = TRUE,filename = paste0("Data_derived/Trace_plots/",SoE_results$scales_results$species[1],"_",SoE_results$niterations,"iters","_",SoE_results$nthin,"thin")) # change path accordingly, for BCR you will need to paste the BCR name also

}

# For BCR the code is:
y<-list.files(path="Data_derived/BCR_results",pattern="\\.Rds$",full.names = T)
rds_list <- lapply(y, readRDS)
for (i in seq_along(y)) {
  
  SoE_results <- rds_list[[i]]
  
  # Get filename only
  fname <- basename(y[i])
  
  # Extract species and region
  species <- strsplit(fname, "_")[[1]][1]
  region  <- strsplit(fname, "_")[[1]][2]
  
  # Clean filename components
  species_file <- gsub("[^[:alnum:] ]", "", species)
  region_file  <- gsub("[^[:alnum:] ]", "", region)
  
  outfile <- file.path(
    "Data_derived",
    "Trace_plots_BCR",
    paste0(
      species_file, "_",
      region_file, "_",
      SoE_results$niterations, "iters_",
      SoE_results$nthin, "thin"
    )
  )
  
  cat("Writing:", outfile, "\n")
  
  message("Model outputs for species: ", SoE_results$scales_results$species[1])
  message("presences: ", sum(SoE_results$data$z == 1),
          " absences: ", sum(SoE_results$data$z == 0))
  
  MCMCvis::MCMCtrace(
    SoE_results$samples,
    params = c(cov.pars, scale.pars),
    ISB = FALSE,
    iter = SoE_results$niterations,
    exact = TRUE,
    pdf = TRUE,
    filename = outfile
  )
} # ignore the warnings() 
