# to create caterpillar plots to check them
## analyzing some results
cov.pars<- c("beta.0","beta.developed_low_intensity","beta.developed_high_intensity",                   # include the quadratic terms if you want to check them too
             "beta.barren_land","beta.deciduous_forest","beta.evergreen_forest","beta.mixed_forest","beta.shrub_scrub","beta.grassland_herbaceous",
             "beta.pasture_hay","beta.cultivated_crop","beta.woody_wetland","beta.emergent_wetland",
             "beta.ocean","beta.lakes","beta.rivers",
             "beta.elevation", "beta.aspect", "beta.slope",
             "beta.precipitation","beta.max_temperature", "beta.snow",
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
MCMCvis::MCMCtrace(SoE_results$samples, 
                   params = cov.pars, 
                   ISB = FALSE, iter=30000,
                   exact = TRUE,
                   pdf = TRUE)

stats::acf(SoE_results$samples[,"beta.max_temperature"])
