library(targets)
library(tarchetypes)

tar_option_set(
  packages = c("tidyverse", "lubridate", "arrow","ggdensity",
               "ggplot2", "cowplot", "ggokabeito", "treeio", "ggtree",
               "jsonlite", "broom", "magrittr", "gt", "yardstick"),
  
)

tar_source(c("R/data_manipulation_functions.R",
           "R/plot_functions.R", 
           "R/table_functions.R"))

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
                   type = if_else(same_genus_truth, "Within genus","Between genera"),
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
                   type = if_else(same_genus_truth, "Within genus","Between genera"),
                   n = prettyNum(n, big.mark =" ")
                 )
               ) %>% select(-n) %>% deframe(),
             format = "qs"
  ),
  tar_target(family_label, 
             family_metadata %>% mutate(
               label = glue::glue("italic(\"{fam}\")~(n == {n_genomes})",
                                  fam=stringr::str_remove(Family,"f__"), 
                                  n_genomes = prettyNum(n_genomes, big.mark= " "))
             ) %>%
               left_join(select(genome_metadata, Phylum,Family) %>% unique(),
                         by ="Family") %>% select(Family,label, Phylum)
  ),
  tar_target(mcc_pocpu_family, get_mcc(all_pocpu, per_family = TRUE), format = "qs"),
  tar_target(mcc_pocpu_global, get_mcc(all_pocpu, per_family = FALSE), format = "qs"),
  tar_target(mcc_pocp_global, get_mcc(all_pocp, per_family = FALSE), format = "qs"),
  tar_target(pocpu_ridges, 
             all_pocpu %>% 
               filter(same_genus_truth) %>%
               left_join(genome_metadata %>% select(Phylum,Family) %>% unique(), by ="Family") %>%
               mutate(Family = glue::glue("italic(\"{Family}\")") %>% as_factor(),
                      Family = fct_reorder(Family, pocp))
  ),
  tar_quarto(slides_retreat, "2024-07-10_RetreatSlidesPOCP.qmd"),
  tar_target(blast_vs_all_pocp,
             pivot_pocp(pocp_values,family_metadata, type = "POCP"),
             format = "parquet"
  ),
  tar_target(blast_vs_all_pocpu,
             pivot_pocp(pocp_values,family_metadata, type = "POCPu"),
             format = "parquet"
  ),
  tar_target(blast_vs_all_pocp_R2, get_lm_R2(blast_vs_all_pocp, type = "POCP"),
             format = "qs"),
  tar_target(blast_vs_all_pocpu_R2, get_lm_R2(blast_vs_all_pocpu, type = "POCPu"),
             format = "qs"),
  tar_target(R2_table, format_R2_table(blast_vs_all_pocp_R2,blast_vs_all_pocpu_R2),
             format = "qs"),
  tar_target(fig_blast_vs_all_pocp,
                    blast_vs_all_pocp %>% arrange(desc(pocp)) %>% 
                      plot_pocp_vs_blast("POCP", R2_table),
             format = "qs"
  ),
  tar_target(fig_blast_vs_all_pocpu,
                    blast_vs_all_pocpu %>% arrange(desc(pocp)) %>% 
                      plot_pocp_vs_blast("POCPu", R2_table),
             format = "qs"
  ),
  tar_target(fig_blast_vs_blastdb,
             plot_grid( blast_vs_all_pocp %>% arrange(desc(pocp)) %>% 
                          plot_pocp_blastdb("POCP"),
                        blast_vs_all_pocpu %>% arrange(desc(pocp)) %>% 
                          plot_pocp_blastdb("POCPu"),
                        ncol = 2, labels = "AUTO"), format = "qs"),
  tar_target(tree_file, "shorlisted_genomes.newick", format = "file"),
  tar_target(shortlisted_tree, treeio::read.tree(file = tree_file), format = "qs"),
  tar_target(tree_metadata,
             genome_metadata %>% select(accession,Domain:Species) %>% 
             left_join(
               select(family_metadata, Family, benchmark_type), by = "Family"
               ) %>% 
               rename("label"="accession") %>% 
               mutate(benchmark_type = forcats::as_factor(benchmark_type) %>% 
                        forcats::fct_recode(
                          "All approaches"="full", "Recommended approach"="recommended"
                        ),
                      across(Domain:Species, ~ str_remove(.x, "[dpcofgs]__"))
                      ),
             format = "qs"),
  tar_target(fig_tree, plot_tree(shortlisted_tree,tree_metadata), format = "qs"),
  tar_target(fig_phyla_count, plot_phyla_count(tree_metadata), format = "qs"),
  tar_target(lpsn_stats_file, "lpsn30.json", format = "file"),
  tar_target(lpsn_stats, parse_lpsn_stats(lpsn_stats_file), format = "qs"),
  tar_target(fig_lpsn_stats, plot_lpsn_stats(lpsn_stats), format = "qs"),
  tar_file(computing_metrics_parquet, "pocpbenchmark_computing_metrics.parquet"),
  tar_target(computing_metrics, read_parquet(computing_metrics_parquet), format = "parquet"),
  tar_target(computing_metrics_fullbenchmark,
             computing_metrics %>% left_join(
               select(family_metadata, Family, benchmark_type), by = "Family"
             ) %>% filter(benchmark_type == "full") %>% select(-benchmark_type),
             format = "parquet"),
  tar_target(db_parsed, get_db_parsed_stats(computing_metrics_fullbenchmark), format = "qs"),
  tar_target(tool_parsed, get_tool_parsed_stats(computing_metrics_fullbenchmark), format = "qs"),
  tar_target(median_db, generate_table_db(db_parsed), format = "qs"),
  tar_target(median_tool, generate_table_tool(tool_parsed), format = "qs"),
  tar_target(db_table, format_db_table(median_db), format = "qs"),
  tar_target(tool_table, format_tool_table(median_tool), format = "qs"),
  tar_target(p_pocp, plot_pocp_density(all_pocp), format = "qs"),
  tar_target(p_pocpu, plot_pocpu_density(all_pocpu), format = "qs"),
  tar_target(p_mcc, plot_mcc(mcc_pocpu_family, family_label, mcc_pocpu_global), format = "qs"),
  tar_target(p_mcc_random, plot_mcc_random(mcc_pocpu_family, family_label, mcc_pocpu_global), format = "qs"),
  tar_target(p_genus_delineation, plot_genus_delineation(p_pocp,p_pocpu, p_mcc), format = "qs"),
  tar_target(pocp_confusion, count(all_pocp, class) %>% deframe() %>% prettyNum(big.mark =" "),
             format = "qs"),
  tar_target(pocpu_confusion, count(all_pocpu, class) %>% deframe() %>% prettyNum(big.mark =" "),
             format = "qs"),
  tar_target(pocpu_confusion_by_family, count(all_pocpu, Family, class), format = "qs"),
  tar_target(lactobacillaceae,
             get_family_confusion_matrix(pocpu_confusion_by_family, "f__Lactobacillaceae"),
             format = "qs"),
  tar_target(streptomycetaceae,
             get_family_confusion_matrix(pocpu_confusion_by_family, "f__Streptomycetaceae"),
             format = "qs"),
  tar_quarto(manuscript)
)