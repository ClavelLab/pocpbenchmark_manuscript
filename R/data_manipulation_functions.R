# Data manipulation functions for the analysis {targets} workflow


# Pivot POCP values to a longer form to compare with BLASTP
pivot_pocp <- function(pocp_values,family_metadata, type = c("POCP", "POCPu")){
  families_to_keep <- dplyr::filter(family_metadata, benchmark_type == "full")
  dplyr::filter(pocp_values, type == {{ type }}) %>% 
    dplyr::semi_join(families_to_keep, by = "Family") %>% 
    dplyr::select(type, tool, is_recommended_tool, pocp, query, subject, Family) %>% 
    # Widen the data with all tools as column and prepend pocp or pocpu to the tool
    pivot_wider(names_from = tool, values_from = pocp, id_cols = c(type, Family,query, subject)) %>%
    # Select only the columns of diamond and mmseqs
    pivot_longer(cols = -c(query, subject,type,Family, blast_blastp), names_to = "tool", values_to = "pocp") %>%
    separate(tool, into = c("aligner", "parameter"), sep = "_", remove = FALSE)
}

# Helpers to replace values or get variable in order to parse computing metrics
# 
replace_ms <- function(chr_time){
  require(magrittr)
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
  require(magrittr)
  dplyr::if_else(stringr::str_detect(id, "-"),
                 true = stringr::str_split(id, "-") %>%
                   sapply(.,function(x)str_sort(x) %>% paste(collapse = "-")),
                 false = id
  )
}

# Parse the nextflow execution trace
parse_computing_metrics <- function(computing_metrics){
  require(magrittr)
  computing_metrics %>% mutate(
    comparison_id =get_comparison_id(dataset_id),
    time = replace_ms(realtime) %>% 
      parse_date_time(orders =  c("S","MS")) %>%
      difftime(parse_date_time("0", orders = "S")) %>% as.duration(),
    memory = rlang::parse_bytes(peak_vmem) %>% as.numeric(),
    cpu = replace_percent(`%cpu`),
    io = rlang::parse_bytes(rchar) %>% as.numeric() +
      rlang::parse_bytes(wchar) %>% as.numeric()
  ) %>% select(family,category, tool, comparison_id, dataset_id, time, memory, cpu, io) 
}
