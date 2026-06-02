library(dplyr)


summary_values <- function(df_scale, cols, final_traits, stat = "max"){
  
  stat_fun <- switch(stat,
                     max = max,
                     min = min,
                     median = median,
                     sd = sd,
                     stop("stat must be one of: max, min, median, sd"))
  
  df_summary <- df_scale %>%
    select(species, {{cols}}) %>%
    rowwise() %>%
    mutate(value = stat_fun(c_across({{cols}}), na.rm = TRUE)) %>%
    ungroup() %>%
    dplyr::select(species, value)
  
  final_traits %>%
    left_join(df_summary, by = "species") %>%
    select(-species)
}