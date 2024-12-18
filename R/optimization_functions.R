# Optimization functions supporting {targets} analysis workflow


# Compute the optimized POCP threshold allowing the maximum MCC score
mcc_optim <- Vectorize(function(data,th){
  data %>%
    mutate(prediction = forcats::as_factor(.data[["pocp"]] > th),
           same_genus_truth = forcats::as_factor(.data[["same_genus_truth"]])) %>% 
    yardstick::mcc(
      truth = all_of("same_genus_truth"),
      estimate = all_of("prediction")
    ) %>% pull(all_of(".estimate"))
}, vectorize.args = "th")

find_pocp_threshold <- function(df){
  ranges <- df %>% pull(all_of("pocp")) %>% range()
  optimise(mcc_optim, interval = ranges, data = df, maximum = TRUE) %>%
    as_tibble() %>% 
    rename("optimized_threshold"="maximum", "maximum_mcc"="objective")
}

# Getting the POCP threshold and associated MCC max value,
# merged together with the MCC with current threshold at 50%
get_optimized_pocp_threshold <- function(df, mcc_df_per_family, family_label){
  df %>% group_by(Family) %>%
    nest() %>%
    mutate(optimisation = map(data, find_pocp_threshold)) %>%
    unnest(optimisation) %>%
    ungroup() %>% select(-data) %>%
    full_join(mcc_df_per_family, by="Family") %>% # with current MCC
    left_join(family_label, by="Family") %>% # with Phylum information
    select(-c(label, mcc_random)) %>%
    mutate(
      threshold_change = optimized_threshold-50,
      mcc_change = round(maximum_mcc-mcc, digits = 2),
      improved_classification = abs(mcc_change) > 0.1,
    ) %>% 
    select(Phylum, Family,mcc_change, threshold_change, improved_classification, optimized_threshold, mcc, maximum_mcc)
  
}
