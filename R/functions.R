read_stats <- function(path_to_family_dir){
  family <- stringr::str_remove(path_to_family_dir,"benchmark-gtdb-")
  prot <- readr::read_tsv(paste0(path_to_family_dir, "/proteins_statistics.tsv"), show_col_types = FALSE)
  dplyr::mutate(prot, family = family)
}

read_R2 <- function(path_to_family_dir){
  family <- stringr::str_remove(path_to_family_dir,"benchmark-gtdb-")
  pocp <- readr::read_csv(paste0(path_to_family_dir, "/compare_pocp/blast-vs-all-pocp-r2.csv"), show_col_types = FALSE)
  pocpu <- readr::read_csv(paste0(path_to_family_dir, "/compare_pocp/blast-vs-all-pocpu-r2.csv"), show_col_types = FALSE)
  dplyr::bind_rows(pocp, pocpu)%>% select(type, tool, r.squared)%>%dplyr::mutate(family = family)
}
read_pocp <- function(path_to_family_dir){
  family <- stringr::str_remove(path_to_family_dir,"benchmark-gtdb-")
  pocp_values <- readr::read_csv(paste0(path_to_family_dir, "/eval_genus_delineation/comparisons_classification_pocp_rand.csv"),
                                 show_col_types = FALSE)
  pocp_values %>%
    dplyr::filter(tool != "blast_blastp") %>% 
    dplyr::mutate(family = family)
}



plot_pocp_distribution<-function(pocp_values, type = c("POCP", "POCPu")){
  ggplot(data = pocp_values, aes(x = pocp, fill = same_genus_truth))+geom_density(alpha=0.5)+
    facet_wrap(~tool)+theme_cowplot()+scale_fill_okabe_ito()+
    geom_vline(xintercept = 50, linetype="dashed")+
    theme(legend.position = "bottom")+
    labs(x=type)
}


read_compute_stats <- function(path_to_family_dir){
  require(magrittr)
  family <- stringr::str_remove(path_to_family_dir,"benchmark-gtdb-")
  path_to_compute_stats <- list.files(paste0(path_to_family_dir, "/pipeline_info"),
                                      pattern = "execution_trace", full.names = TRUE) %>% 
    sort(decreasing = TRUE)
  compute_stats <- readr::read_tsv(path_to_compute_stats[1], show_col_types = FALSE) %>% 
    filter(str_detect(name, "BLAST|DIAMOND|MMSEQS2") & exit == 0) %>%
    select(name, realtime, `%cpu`, peak_vmem, wchar, rchar) %>%
    separate_wider_delim(name, delim = " ",
                         names = c("name_id", "dataset_id")) %>%
    mutate(dataset_id = str_remove(dataset_id, "^\\(") %>% str_remove("\\)$")) %>%
    separate_wider_delim(name_id, delim = ":",names = c(NA, "category", "tool"))
  compute_stats%>% 
    dplyr::mutate(family = family)
}

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


parse_compute_stats <- function(compute_stats){
  require(magrittr)
  compute_stats %>% mutate(
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

get_db_parsed_stats <- function(compute_stats){
  compute_stats %>%
    parse_compute_stats() %>% 
    filter(tool %in% c("MMSEQS2_CREATEDB", "BLAST_MAKEBLASTDB", "DIAMOND_MAKEDB")) %>% 
    select(-dataset_id)
}
get_tool_parsed_stats <- function(compute_stats){
  compute_stats %>%
    parse_compute_stats() %>% 
    filter(!tool %in% c("MMSEQS2_CREATEDB", "BLAST_MAKEBLASTDB", "DIAMOND_MAKEDB")) %>% 
    mutate(
      query = str_remove(dataset_id, "-.*")
    )
}

generate_table_tool <- function(compute_stats_tool){
  compute_stats_tool %>% group_by(dataset_id) %>%
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

generate_table_db <- function(compute_stats_db){
  compute_stats_db %>%
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
  tool_table %>%
    mutate(
      tool= factor(tool, levels = c(
        "BLAST_BLASTP", "BLAST_BLASTPDB",
        "DIAMOND_FAST", "DIAMOND_SENSITIVE",
        "DIAMOND_VERYSENSITIVE", "DIAMOND_ULTRASENSITIVE",
        "MMSEQS2_S1DOT0","MMSEQS2_S2DOT5",
        "MMSEQS2_S6DOT0", "MMSEQS2_S7DOT5"))
    ) %>% arrange(tool)
}
