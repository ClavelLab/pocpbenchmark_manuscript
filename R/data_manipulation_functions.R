# Data manipulation functions for the analysis {targets} workflow

# Ensure the tool variable has the correct levels
format_factor_tool <- function(tool_table){
  tool_table %>%
    mutate(
      tool= factor(tool, levels = c(
        "BLAST_BLASTP", "BLAST_BLASTPDB",
        "DIAMOND_FAST", "DIAMOND_SENSITIVE",
        "DIAMOND_VERYSENSITIVE", "DIAMOND_ULTRASENSITIVE",
        "MMSEQS2_S1DOT0","MMSEQS2_S2DOT5",
        "MMSEQS2_S6DOT0", "MMSEQS2_S7DOT5"))
    )
}

# Pivot POCP values to a longer form to compare with BLASTP
pivot_pocp <- function(pocp_values,family_metadata, type = c("POCP", "POCPu")){
  families_to_keep <- dplyr::filter(family_metadata, benchmark_type == "full")
  dplyr::filter(pocp_values, type == {{ type }}) %>% 
    dplyr::semi_join(families_to_keep, by = "Family") %>% 
    dplyr::select(type, tool, is_recommended_tool, pocp, query, subject, Family) %>% 
    # Widen the data with all tools as column and prepend pocp or pocpu to the tool
    pivot_wider(names_from = tool, values_from = pocp, id_cols = c(type, Family,query, subject)) %>%
    # Select only the columns of diamond and mmseqs
    pivot_longer(cols = -c(query, subject,type,Family, BLAST_BLASTP), names_to = "tool", values_to = "pocp") %>%
    separate(tool, into = c("aligner", "parameter"), sep = "_", remove = FALSE) %>% 
    format_factor_tool()
}

# Helpers to replace values or get variable in order to parse computing metrics
# 
replace_ms <- function(chr_time){
  dplyr::if_else(stringr::str_detect(chr_time, "ms"),
                 true = stringr::str_remove(chr_time, "ms") %>%
                   strtoi() %>% sapply(.,function(x)x*0.001)%>% paste0("s"),
                 false = chr_time
  )
}
replace_percent <- function(chr_percent) {
  as.numeric(gsub("%", "", chr_percent))
}
get_comparison_id <- function(id){
  dplyr::if_else(stringr::str_detect(id, "-"),
                 true = stringr::str_split(id, "-") %>%
                   sapply(.,function(x)str_sort(x) %>% paste(collapse = "-")),
                 false = id
  )
}

# Parse the nextflow execution trace
parse_computing_metrics <- function(computing_metrics){
  computing_metrics %>% mutate(
    comparison_id =get_comparison_id(dataset_id),
    time = replace_ms(realtime) %>% 
      parse_date_time(orders =  c("S","MS")) %>%
      difftime(parse_date_time("0", orders = "S")) %>% as.duration(),
    memory = rlang::parse_bytes(peak_vmem) %>% as.numeric(),
    cpu = replace_percent(`%cpu`),
    io = rlang::parse_bytes(rchar) %>% as.numeric() +
      rlang::parse_bytes(wchar) %>% as.numeric()
  ) %>% select(Family,category, tool, comparison_id, dataset_id, time, memory, cpu, io)
}

# From the figure <https://lpsn.dsmz.de/statistics/figure/30>
# the data used can be guessed to be at <https://lpsn.dsmz.de/statistics/data/30>
# downloaded manually and fix json to be parse correctly for genera stats
# Cannot use tarchetypes::tar_download as the LPSN does not provide a
# ETag or Last-Modified for url
# see: https://github.com/ropensci/tarchetypes/discussions/114
parse_lpsn_stats <- function(lpsn_json){
  readLines(lpsn_json) %>% 
    paste0(collapse = " ") %>% # fix broken strings to get valid json to be parsed
    jsonlite::fromJSON() %>%
    select(amount_is,grouping,subgroup,amount) %>% 
    arrange(grouping) %>% 
    filter(subgroup == "genus") %>% 
    mutate(
      n = cumsum(amount),
      grouping=as.factor(grouping)
    ) 
}

# Get Matthews Correlation Coefficient for classification and random classification
get_mcc <- function(all_pocpu, per_family = TRUE){
  all_pocpu %>% 
  select(type, tool, starts_with("same_"), Family) %>% 
    mutate(across(starts_with("same_"), ~forcats::as_factor(.x))) %>%
    {if(per_family) group_by(., Family) else . } %>%
    summarise(.,
      mcc = mcc_vec(truth = same_genus_truth, estimate = same_genus),
      mcc_random = mcc_vec(truth = same_genus_truth, estimate = same_genus_random)
    )
}

# Helper function to format the count of the confusion matrix for a given Family
get_family_confusion_matrix <- function(df, Family){
  filter(df, Family == {{ Family }}) %>%
    arrange(desc(class)) %>%
    mutate(n = prettyNum(n,big.mark = ",")) %>%
    glue::glue_data("{class} = {n}") %>%
    glue::glue_collapse(sep = ", ", last = " and ")
}
