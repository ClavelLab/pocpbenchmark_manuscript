library(targets)
library(tarchetypes)

tar_option_set(
  packages = c("tidyverse", "lubridate","arrow", "archive"),
  
)

tar_source("R/prepare_data_functions.R")

# End this file with a list of target objects.
list(
  tar_files_input(archives,list.files("data_benchmark", full.names = T)),
  tar_target(extracted, unzip(archives), pattern = map(archives), format = "file"),
  tar_files_input(prots, list.files(pattern = "benchmark-gtdb-f*"), format = "file"),
  tar_target(prot_stats, read_protein_stats(prots), pattern = map(prots),iteration = "vector"),
  tar_target(shortlist_path, "shortlisted_genomes.csv", format = "file"),
  tar_target(gtdb_metadata_url, "https://data.gtdb.ecogenomic.org/releases/release214/214.0/bac120_metadata_r214.tar.gz",
             format = "url"),
  tar_target(gtdb_metadata, readr::read_tsv(archive_read(gtdb_metadata_url))),
  tar_target(genome_metadata,
             readr::read_csv(shortlist_path, show_col_types = FALSE) %>%
               full_join(prot_stats, by = "accession") %>% 
               left_join(gtdb_metadata %>% select(accession, genome_size), by = "accession")),
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
  tar_target(computing_metrics, read_computing_metrics(prots),
             pattern = map(prots),iteration = "vector", format = "qs"),
  # 3. Parquet file for computing metrics
  tar_target(computing_metrics_parquet, 
             tibble_to_parquet(computing_metrics, "pocpbenchmark_computing_metrics.parquet")
             ),
  tar_target(pocp_values, read_pocp(prots), pattern = map(prots), format = "qs"),
  # 4. Parquet for all computed POCP values
  tar_target(pocp_values_parquet,
             tibble_to_parquet(pocp_values, "pocpbenchmark_pocp_values.parquet")
  )
)
