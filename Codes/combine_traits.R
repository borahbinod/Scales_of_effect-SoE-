combine_traits <- function(df_scale, final_traits, scale_cols) {
  library(dplyr)
  
  # Check that requested scale columns exist
  missing_cols <- setdiff(scale_cols, names(df_scale))
  if(length(missing_cols) > 0) {
    stop("These columns are missing in df_scale: ", paste(missing_cols, collapse = ", "))
  }
  
  # Select scale columns + species for joining
  df_scale_sel <- df_scale %>% dplyr::select(species, all_of(scale_cols))
  
  # Select all traits except species
  traits_sel <- final_traits %>% dplyr::select(-species)
  
  # Join by species
  result <- df_scale_sel %>% bind_cols(traits_sel[match(df_scale_sel$species, final_traits$species), ])
  result<- result[,-1]
  names(result)[1] <- "value"
  return(result)
}

