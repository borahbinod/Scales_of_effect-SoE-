lapply(c("randomForest", "tidymodels", "tidyverse", "caret", "lightgbm","mice), library, character.only = TRUE)

df_scale<- read.csv("Data_derived/NA_Birds_Estimated_Breeding_SoE_radius.csv",header=T) # read the selected scales output
traits<- read.csv("Traits/species_traits.csv",header = T)# # read the traits list of those species
traits <- traits %>%
  dplyr::select(where(~ sum(is.na(.)) <= 10)) # select traits with 10 or less missing values
set.seed(100)
traits_id <- traits$species        # keep species separately
traits_data <- traits %>% dplyr::select(-species)   # only trait columns
imp <- mice(traits_data, m = 5, method = "pmm", seed = 123) # imputing for missing values using multiple imputations by chained equations (MICE)
traits_imputed <- complete(imp, 1) # selecting the first set of imputed values, average of all sets can also be selected
final_traits <- cbind(species = traits_id, traits_imputed) # binding imputed traits with corresponding species

source("Codes/combine_traits.R") # this code combines the species traits with the scale of effect of an environmental variable
new_df <- combine_traits(df_scale, final_traits, scale_cols = c("selected_max_temperature")) # this is an example
source("Codes/filter_by_probability.R") # this code filters selected scale of effect values by their posterior probability
new_df<- filter_by_probability(df_scale = df_scale,final_traits = final_traits,habitat ="developed_high_intensity",prob_threshold = 0.2,direction = "greater" ) #this is an example
source("Codes/summary_values.R") # this code summarizes scale of effect values among two or more environmental variables
new_df<- summary_values(df_scale = df_scale,cols = c(selected_elevation,selected_slope,selected_aspect),final_traits = final_traits,stat = "min") # this is an example
# predictive model
sp_one<- new_df
sp_one <- sp_one %>%
  mutate(
    across(
      where(is.numeric) & !all_of("value"),
      ~ as.numeric(scale(.))
    )
  )

train_index <- createDataPartition(sp_one$value, p = 0.8, list = FALSE) # data is split into 80:20 train and test data
train_data <- sp_one[train_index, ]
test_data  <- sp_one[-train_index, ]

X_train <- train_data %>% dplyr::select(-value) # predictors in train data
y_train <- train_data$value # responses in train data

X_test  <- test_data %>% dplyr::select(-value) # predictors in train data
y_test  <- test_data$value # responses in train data

factor_cols_train <- names(X_train)[sapply(X_train, is.factor)]
factor_cols_test  <- names(X_test)[sapply(X_test, is.factor)]

# Ensure the same encoding in train and test
X_train[factor_cols_train] <- lapply(X_train[factor_cols_train], function(x) as.integer(as.factor(x))) # factors converted into integer codes
X_test[factor_cols_test]   <- lapply(X_test[factor_cols_test], function(x) as.integer(as.factor(x)))

y_train_log <- log(y_train) # log transforamtion of responses
y_test_log  <- log(y_test)

X_train_mat <- data.matrix(X_train)
X_test_mat  <- data.matrix(X_test)

dtrain <- lgb.Dataset(data = X_train_mat, label = y_train_log) # creating lightgbm datasets
dtest  <- lgb.Dataset(data = X_test_mat, label = y_test_log)

param_grid <- expand.grid(  #  Define parameter grid
  num_leaves = c(15, 31, 63),
  learning_rate = c(0.01, 0.05, 0.1),
  max_depth = c(-1, 5, 10),
  min_data_in_leaf = c(10, 20, 50),
  feature_fraction = c(0.8, 1.0)
)
param_grid[] <- lapply(param_grid, as.numeric)  #Iteratively test parameters, ensure numeric
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
results_df <- do.call(rbind, results) # Combine results

best_params <- results_df[which.min(results_df$best_score), ] # Best parameters
print(best_params)
final_params <- as.list(best_params[ , names(best_params) %in% names(param_grid)]) # Prepare final parameters for LightGBM
final_params$objective <- "regression"
final_params$metric <- "rmse"
final_params$verbosity <- -1

# Alternative: gradual learning , this values were tested separately instead of the best parameter values derived from iterations
final_params$learning_rate <- 0.01
final_params$min_data_in_leaf <- 5
final_params$max_depth <- 10
final_params$num_leaves <- 63
final_params$objective <- "regression"

final_model <- lgb.train( # Train final model
  params = final_params,
  data = dtrain,
  nrounds = 500    ##best_params$best_iter
)
pred_log <- predict(final_model, X_test_mat) # Predictions
pred <- exp(pred_log)  # back-transform
rmse <- sqrt(mean((pred - y_test)^2)) # Evaluate
cat("Test RMSE:", rmse, "\n")
