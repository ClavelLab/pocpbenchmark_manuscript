# R functions to support the prepare data {targets} workflow

# Read protein stats from seqkit reports
read_protein_stats <- function(path_to_family_dir){
  family <- stringr::str_remove(path_to_family_dir,"benchmark-gtdb-")
  prot <- readr::read_tsv(paste0(path_to_family_dir, "/proteins_statistics.tsv"), show_col_types = FALSE)
  prot %>% dplyr::mutate(
    accession = stringr::str_remove(file,".faa$"),
    family = family
  ) %>%
    dplyr::select(
      accession,family, num_seqs, sum_len, min_len, avg_len, max_len
    )
}

# Convert a tibble into optimized table parquet file
tibble_to_parquet <- function(tibble, parquet_name){
  arrow::write_parquet(tibble, parquet_name)
  return(parquet_name)
}

# Read in the trace files created by nextflow
read_computing_metrics <- function(path_to_family_dir){
  require(magrittr)
  family <- stringr::str_remove(path_to_family_dir,"benchmark-gtdb-")
  path_to_computing_metrics <- list.files(paste0(path_to_family_dir, "/pipeline_info"),
                                      pattern = "execution_trace", full.names = TRUE) %>% 
    sort(decreasing = TRUE)
  computing_metrics <- readr::read_tsv(path_to_computing_metrics[1], show_col_types = FALSE) %>% 
    filter(str_detect(name, "BLAST|DIAMOND|MMSEQS2") & exit == 0) %>%
    select(name, realtime, `%cpu`, peak_vmem, wchar, rchar) %>%
    separate_wider_delim(name, delim = " ",
                         names = c("name_id", "dataset_id")) %>%
    mutate(dataset_id = str_remove(dataset_id, "^\\(") %>% str_remove("\\)$")) %>%
    separate_wider_delim(name_id, delim = ":",names = c(NA, "category", "tool"))
  computing_metrics%>% 
    dplyr::mutate(Family = family) %>% 
    dplyr::relocate(Family)
}


# Read in the POCP/POCPu values produced by the workflow
read_pocp <- function(path_to_family_dir){
  family <- stringr::str_remove(path_to_family_dir,"benchmark-gtdb-")
  pocp_values <- readr::read_csv(paste0(path_to_family_dir, "/eval_genus_delineation/comparisons_classification_pocp_rand.csv"),
                                 show_col_types = FALSE)
  pocp_values %>%
    dplyr::mutate(Family = family,
                  tool = stringr::str_to_upper(tool) %>%
                    factor(levels = c(
                      "BLAST_BLASTP", "BLAST_BLASTPDB",
                      "DIAMOND_FAST", "DIAMOND_SENSITIVE",
                      "DIAMOND_VERYSENSITIVE", "DIAMOND_ULTRASENSITIVE",
                      "MMSEQS2_S1DOT0","MMSEQS2_S2DOT5",
                      "MMSEQS2_S6DOT0", "MMSEQS2_S7DOT5")),
                  is_recommended_tool = tool == "DIAMOND_VERYSENSITIVE"
                  ) %>% 
    dplyr::relocate(type, tool, is_recommended_tool, pocp)
}

