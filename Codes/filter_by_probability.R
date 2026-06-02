library(dplyr)
library(stringr)

filter_by_probability<- function(df_scale, final_traits, 
                            habitat, 
                            prob_threshold = 0.1,
                            direction = "greater") {
  
  # Build column names dynamically
  selected_col    <- paste0("selected_", habitat)
  probability_col <- paste0("probability_", habitat)
  
  # Check columns exist
  if (!all(c(selected_col, probability_col) %in% names(df_scale))) {
    stop("Habitat columns not found in df_scale.")
  }
  
  # Apply filtering condition
  if (direction == "greater") {
    df_filtered <- df_scale %>%
      filter(.data[[probability_col]] > prob_threshold)
  } else if (direction == "less") {
    df_filtered <- df_scale %>%
      filter(.data[[probability_col]] < prob_threshold)
  } else {
    stop("direction must be 'greater' or 'less'")
  }
  
  # Keep only species + selected column
  df_filtered <- df_filtered %>%
    dplyr::select(species, !!selected_col)
  
  # Join with traits
  df_joined <- df_filtered %>%
    left_join(final_traits, by = "species")
  # Rename selected_* to "value" and remove species
  df_joined <- df_joined %>%
    rename(value = !!selected_col) %>%
    dplyr::select(-species)
  return(df_joined)
}