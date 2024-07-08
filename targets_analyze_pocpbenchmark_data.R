library(targets)
library(tarchetypes)

tar_option_set(
  packages = c("tidyverse", "lubridate", "arrow",
               "ggplot2", "cowplot", "ggokabeito"),
  
)

tar_source()

# End this file with a list of target objects.
list(
  tar_target(pocp_values_parquet,
             "pocpbenchmark_pocp_values.parquet", format = "file"),
  tar_target(genome_metadata_parquet,
             "pocpbenchmark_genome_metadata.parquet", format = "file"),
  tar_target(family_metadata_parquet,
             "pocpbenchmark_family_metadata.parquet", format = "file"),
  tar_target(pocp_values,
             read_parquet(pocp_values_parquet), format = "parquet"),
  tar_target(genome_metadata,
             read_parquet(genome_metadata_parquet), format = "parquet"),
  tar_target(family_metadata,
             read_parquet(family_metadata_parquet), format = "parquet"),
  tar_target(total_proteins,
             setNames(genome_metadata[["num_seqs"]], genome_metadata[["accession"]])),
  tar_target(all_pocp,
             dplyr::filter(pocp_values, type == "POCP" & is_recommended_tool)%>%
               mutate(
                 query_proteins = total_proteins[query],
                 subject_proteins = total_proteins[subject],
                 delta = abs(query_proteins - subject_proteins)
               ), format = "parquet"),
  tar_target(all_pocpu,
             dplyr::filter(pocp_values, type == "POCPu" & is_recommended_tool)%>%
               mutate(
                 query_proteins = total_proteins[query],
                 subject_proteins = total_proteins[subject],
                 delta = abs(query_proteins - subject_proteins))
             , format = "parquet"),
  tar_target(pocp_group_sizes,
             all_pocp %>% count(same_genus_truth) %>% 
               mutate(
                 label = glue::glue(
                   "{type} (n = {n})",
                   type = if_else(same_genus_truth, "Intra-genera","Inter-genera"),
                   n = prettyNum(n, big.mark =" ")
                 )
               ) %>% select(-n) %>% deframe(),
             format = "qs"
             ),
  tar_target(pocpu_group_sizes,
             all_pocpu %>% count(same_genus_truth) %>% 
               mutate(
                 label = glue::glue(
                   "{type} (n = {n})",
                   type = if_else(same_genus_truth, "Intra-genera","Inter-genera"),
                   n = prettyNum(n, big.mark =" ")
                 )
               ) %>% select(-n) %>% deframe(),
             format = "qs"
  ),
  tar_target(family_label, 
             family_metadata %>% mutate(
               label = glue::glue("{fam} (n = {n_genomes})",
                                  fam=stringr::str_remove(Family,"f__"), 
                                  n_genomes = prettyNum(n_genomes, big.mark= " "))
             ) %>% select(Family,label)
  ),
  tar_target(pocpu_ridges, 
             all_pocpu %>% 
               filter(same_genus_truth) %>%
               left_join(family_label, by ="Family") %>%
               select(-Family) %>% rename("Family"="label") %>% 
               mutate(Family = as_factor(Family),
                      Family = fct_reorder(Family, pocp))
  )
)