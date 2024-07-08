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
  tar_target(
    cpu_hours_path, "cpu_hours_benchmark_type_per_family.csv", format = "file"
  ),
  tar_target(
    cpu_hours, readr::read_csv(cpu_hours_path, show_col_types = FALSE)
  ),
  tar_target(prot_stats_per_family,
             prot_stats %>% dplyr::select(family, num_seqs) %>%
               dplyr::group_by(family) %>%
               summarise(n_genomes = n(), median_proteins = median(num_seqs),
                         min_proteins = min(num_seqs), max_proteins = max(num_seqs))
             ),
  tar_target(family_metadata, cpu_hours %>% left_join(prot_stats_per_family,
                                                     by=c("Family"="family"))),
  # 2. Parquet file with Family metadata (CPUh, benchmark type, min/median/max proteins)
  tar_target(family_metadata_parquet,
             tibble_to_parquet(family_metadata, "pocpbenchmark_family_metadata.parquet")),
  tar_target(r2, read_R2(prots), pattern = map(prots),iteration = "vector"),
  
  tar_target(pocp_values, read_pocp(prots), pattern = map(prots), format = "qs"),
  tar_target(all_pocp,   dplyr::filter(pocp_values, type == "POCP"), format = "qs"),
  tar_target(all_pocpu,   dplyr::filter(pocp_values, type == "POCPu"), format = "qs"),
  tar_target(plot_pocp, plot_pocp_distribution(all_pocp, "POCP"), format = "qs"),
  tar_target(plot_pocpu, plot_pocp_distribution(all_pocpu, "POCPu"), format = "qs"),
  tar_target(computing_metrics, read_computing_metrics(prots),
             pattern = map(prots),iteration = "vector", format = "qs"),
  # 3. Parquet file for computing metrics
  tar_target(computing_metrics_parquet, 
             tibble_to_parquet(computing_metrics, "pocpbenchmark_computing_metrics.parquet")
             )
)
