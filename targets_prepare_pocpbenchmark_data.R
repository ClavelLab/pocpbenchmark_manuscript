library(targets)
library(tarchetypes)

tar_option_set(
  packages = c("tidyverse", "lubridate","arrow"),
  
)

tar_source()

# End this file with a list of target objects.
list(
  tar_files_input(archives,list.files("data_benchmark", full.names = T)),
  tar_target(extracted, unzip(archives), pattern = map(archives), format = "file"),
  tar_files(prots, list.files(pattern = "benchmark-gtdb-f*"), format = "file"),
  tar_target(prot_stats, read_protein_stats(prots), pattern = map(prots),iteration = "vector"),
  tar_target(shortlist_path, "shortlisted_genomes.csv", format = "file"),
  tar_target(genome_metadata,
             readr::read_csv(shortlist_path, show_col_types = FALSE) %>%
               full_join(prot_stats, by = "accession")),
  # 1. Parquet file with Genome metadata
  tar_target(genome_metadata_parquet,
             tibble_to_parquet(genome_metadata, "pocpbenchmark_genome_metadata.parquet"),
                                  format = "file"),
  tar_target(prot_n, prot_stats %>% dplyr::select(family, num_seqs)%>%dplyr::group_by(family) %>%
               summarise(n_genomes = n(), median_proteins = median(num_seqs))),
  tar_target(table_overview, readr::read_csv("cpuhours.csv")%>% left_join(prot_n, by="family")),
  tar_target(r2, read_R2(prots), pattern = map(prots),iteration = "vector"),
  
  tar_target(pocp_values, read_pocp(prots), pattern = map(prots), format = "qs"),
  tar_target(all_pocp,   dplyr::filter(pocp_values, type == "POCP"), format = "qs"),
  tar_target(all_pocpu,   dplyr::filter(pocp_values, type == "POCPu"), format = "qs"),
  tar_target(plot_pocp, plot_pocp_distribution(all_pocp, "POCP"), format = "qs"),
  tar_target(plot_pocpu, plot_pocp_distribution(all_pocpu, "POCPu"), format = "qs"),
  tar_target(compute_stats, read_compute_stats(prots),
             pattern = map(prots),iteration = "vector", format = "qs"),
  tar_target(db_parsed, get_db_parsed_stats(compute_stats),
             pattern = map(compute_stats), iteration = "vector", format = "qs"),
  tar_target(tool_parsed, get_tool_parsed_stats(compute_stats),
             pattern = map(compute_stats), iteration = "vector", format = "qs"),
  tar_target(median_db, generate_table_db(db_parsed), format = "qs"),
  tar_target(median_tool, generate_table_tool(tool_parsed), format = "qs"),
  tar_target(db_table, format_db_table(median_db), format = "qs"),
  tar_target(tool_table, format_tool_table(median_tool), format = "qs")
)
