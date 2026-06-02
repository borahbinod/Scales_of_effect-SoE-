library(randomForest)
library(tidymodels)
library(caret)
library(tidyverse)
library(lightgbm)
library(mice)

# read the selected scales output
#df_scale<- read.csv("Data_derived/csv_outputs/scale_results.csv",header=T)
df_scale<- read.csv("Traits/NA_Birds_Estimated_Breeding_SoE_radius.csv",header=T)

# read the traits list of those species
traits<- read.csv("Traits/species_traits.csv",header = T)# whole list of traits
traitx<- read.csv("Traits/BIRDBASE.csv",header=T)

traits<-traits %>%
  left_join(traitx, by = c("common_name" = "species"))
traits <- traits %>%
  dplyr::select(where(~ sum(is.na(.)) <= 10))
traits<- traits %>% dplyr::select(species=common_name,habitat,primary_nesting,sociality,beak.length_Culmen,beak.length_nares,beak_width,beak_depth,tarsus_length,
                                  wing_length,kipps_distance,secondary1,hand.wing_index,tail_Length,mass,habitat_density,migration,trophic_level,trophic_niche,primary_lifestyle,
                                  ESI,Clutch_avg,movement_type)

traits<- traits %>% dplyr::select(species=common_name,beak.length_Culmen,beak.length_nares,beak_width,beak_depth,tarsus_length,
                                  wing_length,kipps_distance,secondary1,hand.wing_index,tail_Length,mass,habitat_density,migration,trophic_level,trophic_niche,primary_lifestyle
                                 )

###Imputation
set.seed(100)
traits_id <- traits$species        # keep species separately
traits_data <- traits %>% dplyr::select(-species)   # only trait columns
# Run MICE
imp <- mice(traits_data, m = 5, method = "pmm", seed = 123)
traits_imputed <- complete(imp, 1)
final_traits <- cbind(species = traits_id, traits_imputed)

source("Codes/combine_traits.R")
source("Codes/filter_by_probability.R")
source("Codes/summary_values.R")

new_df <- combine_traits(df_scale, final_traits, scale_cols = c("selected_max_temperature"))
new_df<- filter_by_probability(df_scale = df_scale,final_traits = final_traits,habitat ="developed_high_intensity",prob_threshold = 0.2,direction = "greater" )
new_df<- summary_values(df_scale = df_scale,cols = c(selected_elevation,selected_slope,selected_aspect),final_traits = final_traits,stat = "min")



sp_one<- new_df
sp_one <- sp_one %>%
  mutate(
    across(
      where(is.numeric) & !all_of("value"),
      ~ as.numeric(scale(.))
    )
  )

train_index <- createDataPartition(sp_one$value, p = 0.8, list = FALSE)
train_data <- sp_one[train_index, ]
test_data  <- sp_one[-train_index, ]

X_train <- train_data %>% dplyr::select(-value)
y_train <- train_data$value

X_test  <- test_data %>% dplyr::select(-value)
y_test  <- test_data$value

factor_cols_train <- names(X_train)[sapply(X_train, is.factor)]
factor_cols_test  <- names(X_test)[sapply(X_test, is.factor)]

# Ensure the same encoding in train and test
X_train[factor_cols_train] <- lapply(X_train[factor_cols_train], function(x) as.integer(as.factor(x)))
X_test[factor_cols_test]   <- lapply(X_test[factor_cols_test], function(x) as.integer(as.factor(x)))

y_train_log <- log(y_train)
y_test_log  <- log(y_test)

X_train_mat <- data.matrix(X_train)
X_test_mat  <- data.matrix(X_test)

dtrain <- lgb.Dataset(data = X_train_mat, label = y_train_log)
dtest  <- lgb.Dataset(data = X_test_mat, label = y_test_log)
# 3. Define parameter grid
# -------------------------------
param_grid <- expand.grid(
  num_leaves = c(15, 31, 63),
  learning_rate = c(0.01, 0.05, 0.1),
  max_depth = c(-1, 5, 10),
  min_data_in_leaf = c(10, 20, 50),
  feature_fraction = c(0.8, 1.0)
)
# 4. Iteratively test parameters
# -------------------------------
param_grid[] <- lapply(param_grid, as.numeric)  # ensure numeric

results <- list()

for (i in 1:nrow(param_grid)) {
  params <- as.list(param_grid[i, ])
  params$objective <- "regression"
  params$metric <- "rmse"
  params$verbosity <- -1
  params$feature_pre_filter <- FALSE
  
  cat("Testing parameter set", i, "of", nrow(param_grid), "\n")
  
  model <- lgb.train(
    params = params,
    data = dtrain,
    nrounds = 200,
    valids = list(test = dtest),
    early_stopping_rounds = 20,
    verbose = -1
  )
  
  results[[i]] <- cbind(param_grid[i, ],
                        best_iter = model$best_iter,
                        best_score = model$best_score)
}

# Combine results
results_df <- do.call(rbind, results)

# Best parameters
best_params <- results_df[which.min(results_df$best_score), ]
print(best_params)

# Prepare final parameters for LightGBM
final_params <- as.list(best_params[ , names(best_params) %in% names(param_grid)])
final_params$objective <- "regression"
final_params$metric <- "rmse"
final_params$verbosity <- -1

# Alternative: gradual learning 
final_params$learning_rate <- 0.01
final_params$min_data_in_leaf <- 5
final_params$max_depth <- 10
final_params$num_leaves <- 63
final_params$objective <- "regression"

# Train final model
final_model <- lgb.train(
  params = final_params,
  data = dtrain,
  nrounds = 500    ##best_params$best_iter
)

# Predictions
pred_log <- predict(final_model, X_test_mat)
pred <- exp(pred_log)  # back-transform

# Evaluate
rmse <- sqrt(mean((pred - y_test)^2))
cat("Test RMSE:", rmse, "\n")





##########################
# combining with BIRDBASE
traitx<- read.csv("Traits/BIRDBASE.csv",header = T)
traitx <- traitx %>%
  rename(scientific_name = `eBird.Clements..V2024b.`)

##########################
# combining with Bioshifts
bioshift<- read.csv("Traits/BioShifts.csv",header=T)
subcopy<- bioshift %>% dplyr::select(Species,ShiftR,Unit)
subcopy <- subcopy %>%
  mutate(
    ShiftR = ifelse(Unit == "km/year", ShiftR * 1000, ShiftR),
    Unit = "m/year"  # now all units are meters per year
  )%>% dplyr::select(-Unit)
subcopy <- subcopy %>% distinct(Species, .keep_all = TRUE)

traits_combined <-  traits %>%
  left_join(subcopy, by = c("scientific_name" = "Species"))

bioshifts_exposure<- read.csv("Traits/Bioshifts_merge_Exposure_all.csv",header=T)
##############################################################################

### appending AVONICHE
avoniche<- read.csv("Traits/Avoniche2_eBird2021.csv",header=T)
avoniche<- avoniche %>% rename("scientific_name"="Species2")
avoniche<- avoniche %>% dplyr::select(-c(Family2,Order2,Avibase_id,Inferred,MainDiet,MainNiche))
traits_combined<- traits_combined %>% left_join(avoniche,by="scientific_name")
traits_combined <- traits_combined %>%
  distinct()
traits_combined_clean <- traits_combined %>%
  # remove all columns ending with .y
  dplyr::select(-ends_with(".y")) %>%
  # rename columns ending with .x by removing the .x
  rename_with(~ gsub("\\.x$", "", .x), ends_with(".x"))
#########################################################################
## combining brain data
brain<- read_excel("Traits/Brain_data.xls")

traitx_unique <- traitx %>%
  distinct(scientific_name, .keep_all = TRUE)

traits_combined <- traits %>%
  left_join(traitx_unique, by = "scientific_name")
write.csv(traits_combined,"Traits/species_traits.csv",row.names = F)
