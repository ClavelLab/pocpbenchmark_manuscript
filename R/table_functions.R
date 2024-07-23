# Table functions supporting {targets} analysis workflow

# Extract only the computing metrics:
# - for the DB creation steps
# - for the many-vs-many alignment steps
# 
get_db_parsed_stats <- function(computing_metrics){
  computing_metrics %>%
    parse_computing_metrics() %>% 
    filter(tool %in% c("MMSEQS2_CREATEDB", "BLAST_MAKEBLASTDB", "DIAMOND_MAKEDB")) %>% 
    select(-dataset_id)
}
get_tool_parsed_stats <- function(computing_metrics){
  computing_metrics %>%
    parse_computing_metrics() %>% 
    filter(!tool %in% c("MMSEQS2_CREATEDB", "BLAST_MAKEBLASTDB", "DIAMOND_MAKEDB")) %>% 
    mutate(
      query = str_remove(dataset_id, "-.*")
    )
}

# Aggregate the metrics into meaningful indexes for both tool and db tables
generate_table_tool <- function(computing_metrics_tool){
  computing_metrics_tool %>% group_by(dataset_id) %>%
    mutate(
      time_fold = time/time[tool == "BLAST_BLASTP"],
      memory_fold = memory/memory[tool == "BLAST_BLASTP"],
      cpu_fold = cpu/cpu[tool == "BLAST_BLASTP"],
      io_fold = io/io[tool == "BLAST_BLASTP"]
    ) %>% ungroup() %>% 
    group_by(tool) %>% 
    summarise(
      n = n(),
      time_fold = median(time_fold),
      memory_fold = median(memory_fold),
      cpu_fold = median(cpu_fold),
      io_fold = median(io_fold))
}
generate_table_db <- function(computing_metrics_db){
  computing_metrics_db %>%
    group_by(tool) %>%
    summarise(
      n = n(),
      time = as.double(time) %>% median(),
      memory = median(memory),
      cpu = median(cpu),
      io = median(io)
    )
  # mutate(time=as.duration(time), memory=rlang::as_bytes(memory), cpu=(cpu/100) %>% scales::percent(), io=rlang::as_bytes(io)
}

# Helpers to format the order of tools or the type of variables
format_db_table <- function(db_table){
  require(lubridate)
  db_table %>%
    mutate(
      time=lubridate::as.duration(time),
      memory=rlang::as_bytes(memory),
      cpu=(cpu/100) %>% scales::percent(),
      io=rlang::as_bytes(io)
    )
}
format_tool_table <- function(tool_table){
  format_factor_tool(tool_table) %>% arrange(tool)
}
