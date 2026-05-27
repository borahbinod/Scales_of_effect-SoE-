library(lubridate)
library(dplyr)
library(tidyverse)
library (auk)
library (tidyverse)
library (nimble)
library (sf)
library (MCMCvis)
library (terra)

### Data preparation for BLISS
SoE_BLISS <- function(species, bcr, niterations,nchains, burnin, nthin){
  message(paste0("running SoE_BLISS for species:",species," for niterations =", niterations))
  start_time<- lubridate::now()
  # read checklists and observations
  observations<- readRDS("C:/Users/A02291907/Documents/SoE_Project/SoE_Main/Data_files/Data_derived/observations/observations.Rds")
  checklists<- readRDS("C:/Users/A02291907/Documents/SoE_Project/SoE_Main/Data_files/Data_derived/checklists/checklists.Rds")
  
  # read covariates data
  covar<- readRDS("C:/Users/A02291907/Documents/SoE_Project/SoE_Main/Data_files/Data_derived/covariates/covariates.Rds")
  
  # attach covariates and expand to checklists
  checklists <- inner_join(checklists,covar,by="locality_id")
  
  # remove observations without matching checklists
  observations_selected <- semi_join(observations, checklists, by = "checklist_id")
  
  # select a species
  observations_species<- observations_selected %>% filter(common_name %in% species)
  nrow(observations_species)
  
  # read sepecies range map
  x.range<- st_read(paste0("C:/Users/A02291907/Documents/SoE_Project/SoE_Main/GIS_files/Data_derived/bcr_shapefiles/",bcr,".shp"))
  #plot(st_geometry(xy.range))
  layer_mask<- st_as_sf(st_geometry(x.range))
  
  
  # mask observations by species range
  observations_df = st_as_sf(observations_species, coords = c("longitude", "latitude"), crs = 4326)
  observations_df<- st_transform(observations_df, st_crs(layer_mask))
  shape_observations <- st_filter(observations_df, layer_mask)
  nrow(shape_observations)
  #plot(st_geometry(layer_mask))
  #points(observations_df, col = "black")
  #points(shape_observations, col = "green")
  
  # attach and expand covarites to observation data
  observations_covar<- inner_join(shape_observations,covar, by="locality_id") %>% st_drop_geometry()
  
  # mask checklists by species ranged
  checklists_df = st_as_sf(checklists, coords = c("longitude", "latitude"), crs = 4326)
  checklists_df<-st_transform(checklists_df, st_crs(layer_mask))
  shape_checklists<- st_filter(checklists_df , layer_mask)
  #plot(st_geometry(layer_mask))
  #points(checklists_df, col = "black")
  #points(shape_checklists, col = "green")
  shape_checklists<- shape_checklists %>% filter(!is.na(all_species_reported))%>% st_drop_geometry()
  nrow(shape_checklists)
  
  # zerofill
  ebd_sp<- auk::auk_zerofill(observations_covar,sampling_events = shape_checklists)
  ebd_df<- collapse_zerofill(ebd_sp) ## combining ebd and sampling data
  
  ### convert time to decimal
  source("C:/Users/A02291907/Documents/SoE_Project/SoE_Main/Data_files/Codes/time_to_decimal.R")
  
  # clean up variables
  ebd_df <- ebd_df %>%
    mutate(
      effort_hours = duration_minutes/60,
      # convert time to decimal hours since midnight
      hours_of_day = time_to_decimal(time_observations_started),
      # split date into year and day of year
      year = year(observation_date),
      day_of_year = yday(observation_date)
    )
  
  # restricting checklists to  stationary counts less than 5 hours in duration and with fewer than 10 observers.
  # additional filtering
  species_df <- ebd_df %>%
    filter(protocol_type %in% c("Stationary"),
           effort_hours < 5,
           number_observers < 10)
  
  species_df<- species_df[!duplicated(species_df[ , c("locality_id")]),]
  x.p<-length(which(species_df$species_observed==T));x.a<-length(which(species_df$species_observed==F))
  message(paste0 ("#species:", species, " present=", x.p," and"," #species:", species," absent=", x.a))
  
  ##covariate data for Nimble
  vars<- c("open_water","developed_open_space","developed_low_intensity","developed_medium_intensity","developed_high_intensity",
           "barren_land","deciduous_forest","evergreen_forest","mixed_forest","shrub_scrub","grassland_herbaceous","pasture_hay","cultivated_crops",
           "woody_wetlands","emergent_wetlands","land","ocean","river","lake","elevation","slope","aspect",
           "precipitation","temperature_avg","temperature_max","temperature_min")
  vars_sq<- paste(vars,"sq",sep="_")
  df.list<- list()
  ddf.list<- list()
  for (v in 1:length(vars)){
    df<- species_df[, grepl(vars[v], names(species_df))]
    dfx<- df[,c(1:11,seq(12,50,by=2))]
    ldx<- data.frame(matrix(NA,ncol = ncol(dfx),nrow = nrow(dfx)))
    pdx<- data.frame(matrix(NA,ncol = ncol(dfx),nrow = nrow(dfx)))
    for(d in 1:ncol(dfx)){
      ifelse(length(unique(dfx[,d,drop=TRUE]))>2,ldx[,d]<- poly(dfx[,d,drop=TRUE],degree = 2)[,1],ldx[,d]<- (dfx[,d,drop=TRUE] - mean(dfx[,d,drop=TRUE])) / sd(dfx[,d,drop=TRUE]) ^ as.logical(sd(dfx[,d,drop=TRUE])))
      ifelse(length(unique(dfx[,d,drop=TRUE]))>2,pdx[,d]<- poly(dfx[,d,drop=TRUE],degree = 2)[,2],pdx[,d]<- (dfx[,d,drop=TRUE] - mean(dfx[,d,drop=TRUE])) / sd(dfx[,d,drop=TRUE]) ^ as.logical(sd(dfx[,d,drop=TRUE])))
    }
    ldx <- ldx[,colSums(is.na(ldx))<nrow(ldx)]
    pdx <- pdx[,colSums(is.na(pdx))<nrow(pdx)]
    df.list[[v]]<- ldx
    ddf.list[[v]]<- pdx
    assign(vars[v],df.list[[v]])
    assign(vars_sq[v],ddf.list[[v]])
  }
  effort_hours<- as.numeric(scale(species_df$effort_hours)); number_observers<- as.integer(scale(species_df$number_observers)); hours_of_day<- as.numeric(scale(species_df$hours_of_day))
  
  ## response data for Nimble
  sp_observed<- as.integer(as.logical(species_df$species_observed))
  
  ## Nimble model
  Nimble_BLISS<- nimbleCode( {
    # prior for SD of the covariate parameters
    SD_covar <- 1/100
    
    # priors for linear covariates of ecological model
    beta.0 ~ dnorm (0.0, SD_covar) # model intercept 
    beta.developed_low_intensity ~ dnorm (0.0, SD_covar) # model coefficients representing relationships between  covariates and species observation
    beta.developed_high_intensity ~ dnorm (0.0, SD_covar)
    beta.barren_land ~ dnorm (0.0, SD_covar)
    beta.deciduous_forest ~ dnorm (0.0, SD_covar)
    beta.evergreen_forest ~ dnorm (0.0, SD_covar)
    beta.mixed_forest ~ dnorm (0.0, SD_covar)
    beta.shrub_scrub ~ dnorm (0.0, SD_covar)
    beta.grassland_herbaceous ~ dnorm (0.0, SD_covar)
    beta.pasture_hay ~ dnorm (0.0, SD_covar)
    beta.cultivated_crop ~ dnorm (0.0, SD_covar)
    beta.woody_wetland ~ dnorm (0.0, SD_covar)
    beta.emergent_wetland ~ dnorm (0.0, SD_covar)
    beta.ocean ~ dnorm (0.0, SD_covar)
    beta.lakes ~ dnorm (0.0, SD_covar)
    beta.rivers ~ dnorm (0.0, SD_covar)
    beta.elevation ~ dnorm (0.0, SD_covar)
    beta.aspect ~ dnorm (0.0, SD_covar)
    beta.slope ~ dnorm (0.0, SD_covar)
    beta.precipitation ~ dnorm (0.0, SD_covar)
    beta.max_temperature ~ dnorm (0.0, SD_covar)
    beta.effort_hours ~ dnorm (0.0, SD_covar)
    beta.number_observers ~ dnorm (0.0, SD_covar)
    beta.time_of_day ~ dnorm (0.0, SD_covar)
    
    # priors for quadratic terms
    beta.developed_low_intensity_sq ~ dnorm (0.0, SD_covar) # model coefficients representing relationships between  covariates and species observation
    beta.developed_high_intensity_sq ~ dnorm (0.0, SD_covar)
    beta.barren_land_sq ~ dnorm (0.0, SD_covar)
    beta.deciduous_forest_sq ~ dnorm (0.0, SD_covar)
    beta.evergreen_forest_sq ~ dnorm (0.0, SD_covar)
    beta.mixed_forest_sq ~ dnorm (0.0, SD_covar)
    beta.shrub_scrub_sq ~ dnorm (0.0, SD_covar)
    beta.grassland_herbaceous_sq ~ dnorm (0.0, SD_covar)
    beta.pasture_hay_sq ~ dnorm (0.0, SD_covar)
    beta.cultivated_crop_sq ~ dnorm (0.0, SD_covar)
    beta.woody_wetland_sq ~ dnorm (0.0, SD_covar)
    beta.emergent_wetland_sq ~ dnorm (0.0, SD_covar)
    beta.ocean_sq ~ dnorm (0.0, SD_covar)
    beta.lakes_sq ~ dnorm (0.0, SD_covar)
    beta.rivers_sq ~ dnorm (0.0, SD_covar)
    beta.elevation_sq ~ dnorm (0.0, SD_covar)
    beta.aspect_sq ~ dnorm (0.0, SD_covar)
    beta.slope_sq ~ dnorm (0.0, SD_covar)
    beta.precipitation_sq ~ dnorm (0.0, SD_covar)
    beta.max_temperature_sq ~ dnorm (0.0, SD_covar)
    
    
    # scale indicator for each of the covariate that will perform scale selection for each covariate
    scale.1 ~ dcat (priors[1:31]) 
    scale.2 ~ dcat (priors[1:31]) 
    scale.3 ~ dcat (priors[1:31]) 
    scale.4 ~ dcat (priors[1:31]) 
    scale.5 ~ dcat (priors[1:31]) 
    scale.6 ~ dcat (priors[1:31]) 
    scale.7 ~ dcat (priors[1:31]) 
    scale.8 ~ dcat (priors[1:31]) 
    scale.9 ~ dcat (priors[1:31]) 
    scale.10 ~ dcat (priors[1:31]) 
    scale.11 ~ dcat (priors[1:31]) 
    scale.12 ~ dcat (priors[1:31]) 
    scale.13 ~ dcat (priors[1:31]) 
    scale.14 ~ dcat (priors[1:31]) 
    scale.15 ~ dcat (priors[1:31]) 
    scale.16 ~ dcat (priors[1:31]) 
    scale.17 ~ dcat (priors[1:31]) 
    scale.18 ~ dcat (priors[1:31]) 
    scale.19 ~ dcat (priors[1:31]) 
    scale.20 ~ dcat (priors[1:31]) 
    
    ## priors on scales
    for (sc in 1:31){ 
      priors[sc] <- 1/31
    }
    
    
    # Model
    for (i in 1:nlocalities){
      logit_psi[i] <- beta.0 + beta.developed_low_intensity * developed_low_intensity [i, scale.1]+ beta.developed_low_intensity_sq * developed_low_intensity_sq [i, scale.1]  + beta.developed_high_intensity * developed_high_intensity [i, scale.2]+ beta.developed_high_intensity_sq * developed_high_intensity_sq [i, scale.2]+ 
        beta.barren_land * barren_land [i, scale.3]+ beta.barren_land_sq * barren_land_sq [i, scale.3] + beta.deciduous_forest * deciduous_forest [i,scale.4]+ beta.deciduous_forest_sq * deciduous_forest_sq [i,scale.4] + beta.evergreen_forest * evergreen_forest [i, scale.5]+ beta.evergreen_forest_sq * evergreen_forest_sq [i, scale.5] + beta.mixed_forest * mixed_forest [i, scale.6]+ beta.mixed_forest_sq * mixed_forest_sq [i, scale.6]+ 
        beta.shrub_scrub * shrub_scrub [i, scale.7]+ beta.shrub_scrub_sq * shrub_scrub_sq [i, scale.7] + beta.grassland_herbaceous * grassland_herbaceous [i, scale.8]+ beta.grassland_herbaceous_sq * grassland_herbaceous_sq [i, scale.8] + beta.pasture_hay * pasture_hay [i, scale.9]+ beta.pasture_hay_sq * pasture_hay_sq [i, scale.9] + beta.cultivated_crop * cultivated_crop [i, scale.10]+ beta.cultivated_crop_sq * cultivated_crop_sq [i, scale.10] + beta.woody_wetland * woody_wetland [i,scale.11]+ beta.woody_wetland_sq * woody_wetland_sq [i,scale.11] + beta.emergent_wetland * emergent_wetland [i, scale.12]+ beta.emergent_wetland_sq * emergent_wetland_sq [i, scale.12]+
        beta.ocean * ocean [i, scale.13]+ beta.ocean_sq * ocean_sq [i, scale.13] + beta.lakes * lakes [i, scale.14]+ beta.lakes_sq * lakes_sq [i, scale.14] + beta.rivers * rivers [i, scale.15]+ beta.rivers_sq * rivers_sq [i, scale.15]+
        beta.elevation * elevation [i, scale.16]+ beta.elevation_sq * elevation_sq [i, scale.16] + beta.aspect * aspect [i, scale.17]+ beta.aspect_sq * aspect_sq [i, scale.17] + beta.slope * slope [i, scale.18]+ beta.slope_sq * slope_sq [i, scale.18]+ 
        beta.precipitation * precipitation [i, scale.19]+ beta.precipitation_sq * precipitation_sq [i, scale.19] + beta.max_temperature * max_temperature [i, scale.20]+ beta.max_temperature_sq * max_temperature_sq [i, scale.20]+
        beta.effort_hours * effort_hours [i] + beta.number_observers * number_observers [i] + beta.time_of_day * time_of_day [i]
      
      psi[i] <- 1 / (1 + exp(-logit_psi[i]))  
      z[i] ~ dbern (psi[i]) # species observed/ not observed
    }
  })
  
  
  ## Nimble constants
  "nlocalities"= nrow(species_df)
  my_constants <- list(
    nlocalities= nlocalities
  )
  
  ## Nimble data
  observation_data= list (
    "z"= sp_observed,
    "developed_low_intensity"= developed_low_intensity,
    "developed_low_intensity_sq"= developed_low_intensity_sq,
    "developed_high_intensity"= developed_high_intensity,
    "developed_high_intensity_sq"= developed_high_intensity_sq,
    "barren_land"= barren_land,
    "barren_land_sq"= barren_land_sq,
    "deciduous_forest"= deciduous_forest,
    "deciduous_forest_sq"= deciduous_forest_sq,
    "evergreen_forest"= evergreen_forest,
    "evergreen_forest_sq"= evergreen_forest_sq,
    "mixed_forest"= mixed_forest,
    "mixed_forest_sq"= mixed_forest_sq,
    "shrub_scrub"= shrub_scrub,
    "shrub_scrub_sq"= shrub_scrub_sq,
    "grassland_herbaceous"= grassland_herbaceous,
    "grassland_herbaceous_sq"= grassland_herbaceous_sq,
    "pasture_hay"= pasture_hay,
    "pasture_hay_sq"= pasture_hay_sq,
    "cultivated_crop"= cultivated_crops,
    "cultivated_crop_sq"= cultivated_crops_sq,
    "woody_wetland"= woody_wetlands,
    "woody_wetland_sq"= woody_wetlands_sq,
    "emergent_wetland"= emergent_wetlands,
    "emergent_wetland_sq"= emergent_wetlands_sq,
    "ocean"= ocean,
    "ocean_sq"= ocean_sq,
    "rivers"= river,
    "rivers_sq"= river_sq,
    "lakes"= lake,
    "lakes_sq"= lake_sq,
    "elevation"= elevation,
    "elevation_sq"= elevation_sq,
    "slope"= slope,
    "slope_sq"= slope_sq,
    "aspect"= aspect,
    "aspect_sq"= aspect_sq,
    "precipitation"= precipitation,
    "precipitation_sq"= precipitation_sq,
    "max_temperature"= temperature_max,
    "max_temperature_sq"= temperature_max_sq,
    "effort_hours"= effort_hours,
    "number_observers"= number_observers,
    "time_of_day"= hours_of_day
  )
  
  # check data for NA
  if (any(sapply(observation_data, function(x) { any(is.na(x))}))) {
    NA_check = sapply(observation_data, function(x) { any(is.na(x))})
    
    stop(paste0("Some input data are NA:\n", paste(names(NA_check[NA_check]), collapse = "\n")))
  }
  inits <-function(){list(
    beta.0 = dnorm(0,2),
    beta.developed_low_intensity = dnorm(0,2),
    beta.developed_low_intensity_sq = dnorm(0,2),
    beta.developed_high_intensity = dnorm(0,2),
    beta.developed_high_intensity_sq= dnorm(0,2),
    beta.barren_land = dnorm(0,2),
    beta.barren_land_sq = dnorm(0,2),
    beta.deciduous_forest = dnorm(0,2),
    beta.deciduous_forest_sq = dnorm(0,2),
    beta.evergreen_forest = dnorm(0,2),
    beta.evergreen_forest_sq = dnorm(0,2),
    beta.mixed_forest = dnorm(0,2),
    beta.mixed_forest_sq = dnorm(0,2),
    beta.shrub_scrub = dnorm(0,2),
    beta.shrub_scrub_sq = dnorm(0,2),
    beta.grassland_herbaceous = dnorm(0,2),
    beta.grassland_herbaceous_sq = dnorm(0,2),
    beta.pasture_hay = dnorm(0,2),
    beta.pasture_hay_sq = dnorm(0,2),
    beta.cultivated_crop = dnorm(0,2),
    beta.cultivated_crop_sq = dnorm(0,2),
    beta.woody_wetland = dnorm(0,2),
    beta.woody_wetland_sq = dnorm(0,2),
    beta.emergent_wetland = dnorm(0,2),
    beta.emergent_wetland_sq = dnorm(0,2),
    beta.lakes = dnorm(0,2),
    beta.lakes_sq = dnorm(0,2),
    beta.ocean= dnorm(0,2),
    beta.ocean_sq= dnorm(0,2),
    beta.rivers = dnorm(0,2),
    beta.rivers_sq = dnorm(0,2),
    beta.elevation = dnorm(0,2),
    beta.elevation_sq = dnorm(0,2),
    beta.aspect = dnorm(0,2),
    beta.aspect_sq = dnorm(0,2),
    beta.slope = dnorm(0,2),
    beta.slope_sq = dnorm(0,2),
    beta.precipitation = dnorm(0,2),
    beta.precipitation_sq = dnorm(0,2),
    beta.max_temperature = dnorm(0,2),
    beta.max_temperature_sq = dnorm(0,2),
    beta.effort_hours = dnorm(0,2),
    beta.number_observers = dnorm(0,2),
    beta.time_of_day = dnorm(0,2),
    scale.1 = sample.int(31,1), # possible scale values
    scale.2 = sample.int(31,1),
    scale.3 = sample.int(31,1),
    scale.4 = sample.int(31,1), 
    scale.5 = sample.int(31,1),
    scale.6 = sample.int(31,1),
    scale.7 = sample.int(31,1),
    scale.8 = sample.int(31,1),
    scale.9 = sample.int(31,1),
    scale.10 = sample.int(31,1),
    scale.11 = sample.int(31,1),
    scale.12 = sample.int(31,1),
    scale.13 = sample.int(31,1),
    scale.14 = sample.int(31,1),
    scale.15 = sample.int(31,1),
    scale.16 = sample.int(31,1),
    scale.17 = sample.int(31,1),
    scale.18 = sample.int(31,1),
    scale.19 = sample.int(31,1),
    scale.20 = sample.int(31,1)
  )} 
  timestamp()
  samples <- nimbleMCMC(
    code = Nimble_BLISS,
    constants = my_constants, 
    data = observation_data, 
    inits = inits,
    monitors = c(
      "beta.0","beta.developed_low_intensity","beta.developed_low_intensity_sq","beta.developed_high_intensity","beta.developed_high_intensity_sq",
      "beta.barren_land","beta.barren_land_sq","beta.deciduous_forest","beta.deciduous_forest_sq","beta.evergreen_forest","beta.evergreen_forest_sq","beta.mixed_forest","beta.mixed_forest_sq",
      "beta.shrub_scrub","beta.shrub_scrub_sq","beta.grassland_herbaceous","beta.grassland_herbaceous_sq",
      "beta.pasture_hay","beta.pasture_hay_sq","beta.cultivated_crop","beta.cultivated_crop_sq",
      "beta.woody_wetland","beta.woody_wetland_sq","beta.emergent_wetland","beta.emergent_wetland_sq",
      "beta.ocean","beta.ocean_sq","beta.lakes","beta.lakes_sq","beta.rivers","beta.rivers_sq",
      "beta.elevation","beta.elevation_sq", "beta.aspect","beta.aspect_sq","beta.slope","beta.slope_sq",
      "beta.precipitation","beta.precipitation_sq","beta.max_temperature", "beta.max_temperature_sq",
      "beta.effort_hours", "beta.number_observers",
      "beta.time_of_day",
      "scale.1", "scale.2", "scale.3", "scale.4","scale.5", 
      "scale.6", "scale.7", "scale.8","scale.9", "scale.10", 
      "scale.11", "scale.12","scale.13", "scale.14", 
      "scale.15", "scale.16", "scale.17", 
      "scale.18", "scale.19", "scale.20"
    ),
    niter = niterations, ## total saved iterations will be niter - nburnin
    nburnin = burnin,
    nchains = nchains,
    thin = nthin
  )
  gc()
  model_results <- MCMCvis::MCMCsummary(
    samples,
    digits= 2
  ) 
  results_df<- samples
  ## get selected scales
  scale_vars<- read.csv("C:/Users/A02291907/Documents/SoE_Project/SoE_Main/Data_files/Data_derived/covariates/scale_indicators.csv",header=T)
  all_variables <- colnames(results_df)
  scales_df <- results_df[,grepl("scale", all_variables)]
  selected_scales<- vector()
  for (s in 1:ncol(scales_df)){
    selected_scales[s]<- as.integer(names(which.max(table(scales_df[,s]))))
  }
  
  selected_scales_prob<- vector()
  for (s in 1:ncol(scales_df)){
    selected_scales_prob[s] <- proportions(table(scales_df[,s]))[as.integer(names(which.max(table(scales_df[,s]))))] 
  }
  scales_results<- data.frame("scale_indicator"=colnames(scales_df),"selected_scale"=selected_scales, "selected_scale_prob"=selected_scales_prob)
  scales_results<- cbind(merge(x = scales_results,y = scale_vars), species=species)# 
  
  finish_time<- lubridate::now()  
  # return result
  list(
    data= observation_data,
    samples= samples,
    niterations = niterations,
    nchains = nchains,
    burnin = burnin,
    nthin= nthin,
    model_results = model_results,
    scales_results = scales_results,
    range=x.range$name_en,
    start_time = start_time,
    finish_time = finish_time
  )
  
}

SoE_results<- SoE_BLISS(species = "Baltimore Oriole",bcr = "Atlantic Northern Forest",niterations = 160000,nchains = 1,burnin = 10000, nthin = 5)

saveRDS(SoE_results,paste0("Data_derived/BCR_results/", species,"_",bcr,"_Results_", (160000-10000)/5,"iterations",".Rds"))

bcr_sample<- read.csv("Data_derived/bcr_samples.csv",header=T,check.names = FALSE)
names(bcr_sample) <- gsub("\\.", " ", names(bcr_sample))


# columns to ignore
ignore_cols <- c("species", "common_name", "model_counts", "status")

# get BCR columns
bcr_cols <- setdiff(colnames(bcr_sample), ignore_cols)

# loop over each species
for (i in 320:330) {
  
  species_i <- bcr_sample$common_name[i]
  
  row_i <- bcr_sample[i, bcr_cols, drop = FALSE]
  valid_bcrs <- names(row_i)[which(row_i >= 20)]
  
  if (length(valid_bcrs) == 0) {
    message(paste("Skipping", species_i, "- no BCR >= 20"))
    next
  }
  
  for (bcr_i in valid_bcrs) {
    
    message(paste("Running:", species_i, "| BCR:", bcr_i))
    
    # wrap EVERYTHING in tryCatch
    tryCatch({
      
      SoE_results <- SoE_BLISS(
        species = species_i,
        bcr = bcr_i,
        niterations = 160000,
        nchains = 1,
        burnin = 10000,
        nthin = 5
      )
      
      species_clean <- gsub("[^A-Za-z0-9]", "_", species_i)
      bcr_clean <- gsub("[^A-Za-z0-9]", "_", bcr_i)
      
      outfile <- paste0(
        "Data_derived/BCR_results/",
        species_clean, "_",
        bcr_clean,
        "_Results_",
        (160000 - 10000) / 5,
        "iterations.Rds"
      )
      
      saveRDS(SoE_results, outfile)
      
      rm(SoE_results)
      gc()
      
    }, error = function(e) {
      
      message(paste(
        "ERROR for species:", species_i,
        "| BCR:", bcr_i,
        "| Message:", e$message
      ))
      
      # optional: log errors to file
      write(
        paste(Sys.time(), species_i, bcr_i, e$message),
        file = "error_log.txt",
        append = TRUE
      )
      
      # continue loop automatically
      NULL
    })
  }
}




