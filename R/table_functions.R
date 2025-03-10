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

# Fit linear regressions to blast_vs_all comparisons and extract R2
get_lm_R2 <- function(df, type = c("POCP", "POCPu")){
  df %>% group_by(tool) %>%
    nest() %>%
    mutate(
      fit = map(data, ~ lm(BLAST_BLASTP ~ pocp, data = .)),
      tidied = map(fit, broom::glance)
    ) %>%
    unnest(tidied) %>%
    select(-data, -fit) %>%
    ungroup() %>% select(tool,r.squared,p.value,nobs ) %>% 
    rename("R2"="r.squared", "p"="p.value") %>% 
    mutate(
      type = {{ type }},
      p_label = map_chr(p, scales::label_pvalue()),
      label = glue::glue("atop(italic(R)^2 == {R2}, italic(p){p_label})",
                         R2=round(R2,digits = 3) %>% prettyNum()))
}

# Helper to join and format both R2 values for POCP and POCPu regressions
format_R2_table <- function(pocp_table, pocpu_table){
  rename_R2_table <- function(tbl){
    tbl %>% select(tool, R2, p_label, nobs, type) %>%
      rename_with(~paste(unique(tbl$type), .x,sep = "_", recycle0 = TRUE),
                  c(R2,p_label,nobs)) %>%
      select(-type)
  }
  
  full_join(
    rename_R2_table(pocp_table),
    rename_R2_table(pocpu_table),
    by = "tool")
}

format_optimized_pocp_table <- function(optim_df){
  optim_df %>% 
    arrange(Phylum, desc(mcc)) %>% 
    mutate(
      across(Phylum:Family, ~ str_remove(.x, "[p|f]__")),
      need_increase = if_else(improved_classification,
                              if_else(optimized_threshold == 50,
                                      NA, optimized_threshold > 50),
                              NA),
      is_rescued = if_else(improved_classification,
                           mcc <0.7 & maximum_mcc>0.7,
                           NA)
      ) %>% 
    group_by(Phylum) %>% 
    gt::gt() %>% 
    fmt_tf(
      columns = need_increase,
      tf_style = "arrows", na_val = ""
    ) %>%
    fmt_tf(
      columns = is_rescued,
      tf_style = "squares", na_val = ""
    ) %>%
    fmt_number(
      columns = ends_with("mcc"), decimals = 2,
    ) %>%
    fmt_number(
      columns = optimized_threshold, decimals = 1
    ) %>% 
    tab_style(
      style = cell_text(style = "italic"),
      locations = cells_body(columns = c(Phylum, Family))
    ) %>%
    tab_style(
      style = cell_text(style = "italic"),
      locations = cells_row_groups()
    ) %>% 
    cols_label(
      mcc_change = md("$\\Delta$MCC"), 
      mcc = "MCC",maximum_mcc="Max. MCC",
      need_increase = "",
      is_rescued = "",
      optimized_threshold = md("Threshold")
    ) %>%
    tab_options(column_labels.font.weight = "bold") %>% 
    cols_hide(c(threshold_change,improved_classification)) %>%
    cols_move(need_increase, after = optimized_threshold) %>% 
    cols_move_to_start(is_rescued)
}

format_metrics_table <- function(tbl){
  tbl %>%
    mutate(category = str_remove(tool, "_.*")) %>% 
    group_by(category) %>% 
    gt::gt() %>% 
    cols_hide(n) %>% 
    fmt_number(
      columns = gt::ends_with("fold"), decimals = 3
    ) %>% 
    cols_label(
      tool = "Approach name",
      time_fold =  "Time",
      memory_fold = "Memory",
      cpu_fold = "CPU",
      io_fold = "Disk usage (I/O)") %>% 
    cols_align(align = "left", columns = tool) %>% 
    tab_style_body(
      cell_text(weight = "bold"),
      values = "DIAMOND_VERYSENSITIVE",
      targets = "row"
    ) %>% 
    tab_footnote(
      footnote = md("@camachoBLASTArchitectureApplications2009"),
      locations = cells_row_groups("BLAST")
    ) %>% 
    tab_footnote(
      footnote = md("@buchfinkSensitiveProteinAlignments2021"),
      locations = cells_row_groups("DIAMOND")
    ) %>% 
    tab_footnote(
      footnote = md("@steineggerMMseqs2EnablesSensitive2017"),
      locations = cells_row_groups("MMSEQS2")
    )
}
