library(targets)
library(tarchetypes)

tar_option_set(
  packages = c("tidyverse", "lubridate", "arrow","ggdensity",
               "ggplot2", "cowplot", "ggokabeito", "treeio", "ggtree",
               "jsonlite", "broom", "magrittr", "gt", "yardstick"),
  format = "qs"
)

tar_source(c("R/data_manipulation_functions.R",
           "R/plot_functions.R",
           "R/table_functions.R"))

# End this file with a list of target objects.
list(
  tar_file(pocp_values_parquet, "pocpbenchmark_pocp_values.parquet"),
  tar_file(genome_metadata_parquet, "pocpbenchmark_genome_metadata.parquet"),
  tar_file(family_metadata_parquet, "pocpbenchmark_family_metadata.parquet"),
  tar_parquet(pocp_values, read_parquet(pocp_values_parquet)),
  tar_parquet(genome_metadata, read_parquet(genome_metadata_parquet)),
  tar_parquet(family_metadata, read_parquet(family_metadata_parquet)),
  tar_target(total_proteins,
             setNames(genome_metadata[["num_seqs"]], genome_metadata[["accession"]])),
  tar_target(genome_sizes,
             setNames(genome_metadata[["genome_size"]], genome_metadata[["accession"]])),
  tar_parquet(all_pocp,
             dplyr::filter(pocp_values, type == "POCP" & is_recommended_tool)%>%
               mutate(
                 query_proteins = total_proteins[query],
                 subject_proteins = total_proteins[subject],
                 delta_proteome = abs(query_proteins - subject_proteins),
                 delta_genome = abs(genome_sizes[query] - genome_sizes[subject])
               )),
  tar_parquet(all_pocpu,
             dplyr::filter(pocp_values, type == "POCPu" & is_recommended_tool)%>%
               mutate(
                 query_proteins = total_proteins[query],
                 subject_proteins = total_proteins[subject],
                 delta_proteome = abs(query_proteins - subject_proteins),
                 delta_genome = abs(genome_sizes[query] - genome_sizes[subject])
               )),
  tar_target(pocp_group_sizes,
             all_pocp %>% count(same_genus_truth) %>% 
               mutate(
                 label = glue::glue(
                   "{type} (n = {n})",
                   type = if_else(same_genus_truth, "Within genus","Between genera"),
                   n = prettyNum(n, big.mark =" ")
                 )
               ) %>% select(-n) %>% deframe() 
             ),
  tar_target(pocpu_group_sizes,
             all_pocpu %>% count(same_genus_truth) %>% 
               mutate(
                 label = glue::glue(
                   "{type} (n = {n})",
                   type = if_else(same_genus_truth, "Within genus","Between genera"),
                   n = prettyNum(n, big.mark =" ")
                 )
               ) %>% select(-n) %>% deframe() 
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
  tar_target(mcc_pocpu_family, get_mcc(all_pocpu, per_family = TRUE)),
  tar_target(mcc_pocpu_global, get_mcc(all_pocpu, per_family = FALSE)),
  tar_target(mcc_pocp_global, get_mcc(all_pocp, per_family = FALSE)),
  tar_target(pocpu_ridges,
             all_pocpu %>% 
               filter(same_genus_truth) %>%
               left_join(genome_metadata %>% select(Phylum,Family) %>% unique(), by ="Family") %>%
               mutate(Family = glue::glue("italic(\"{Family}\")") %>% as_factor(),
                      Family = fct_reorder(Family, pocp))
  ),
  tar_parquet(blast_vs_all_pocp, pivot_pocp(pocp_values,family_metadata, type = "POCP")),
  tar_parquet(blast_vs_all_pocpu, pivot_pocp(pocp_values,family_metadata, type = "POCPu")),
  tar_target(blast_vs_all_pocp_R2, get_lm_R2(blast_vs_all_pocp, type = "POCP")),
  tar_target(blast_vs_all_pocpu_R2, get_lm_R2(blast_vs_all_pocpu, type = "POCPu")),
  tar_target(R2_table, format_R2_table(blast_vs_all_pocp_R2,blast_vs_all_pocpu_R2)),
  tar_target(fig_blast_vs_all_pocp,
                    blast_vs_all_pocp %>% arrange(desc(pocp)) %>% 
                      plot_pocp_vs_blast("POCP", R2_table)
  ),
  tar_target(fig_blast_vs_all_pocpu,
                    blast_vs_all_pocpu %>% arrange(desc(pocp)) %>% 
                      plot_pocp_vs_blast("POCPu", R2_table) 
  ),
  tar_target(fig_blast_vs_blastdb,
             plot_grid( blast_vs_all_pocp %>% arrange(desc(pocp)) %>% 
                          plot_pocp_blastdb("POCP"),
                        blast_vs_all_pocpu %>% arrange(desc(pocp)) %>% 
                          plot_pocp_blastdb("POCPu"),
                        ncol = 2, labels = "AUTO")),
  tar_file(fig_blast_vs_blastdb_png, save_png(fig_blast_vs_blastdb,
                                              "figures/fig_blast_vs_blastdb.png", 9, 5)),
  tar_file(tree_file, "shorlisted_genomes.newick"),
  tar_target(shortlisted_tree, treeio::read.tree(file = tree_file)),
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
             ),
  tar_target(fig_tree, plot_tree(shortlisted_tree,tree_metadata)),
  tar_target(fig_phyla_count, plot_phyla_count(tree_metadata)),
  tar_target(fig_tree_phyla_count,
             plot_tree_phyla_count(fig_tree,fig_phyla_count,genome_metadata)),
  tar_file(fig_tree_phyla_count_png, save_png(fig_tree_phyla_count, "figures/fig_tree_phyla_count.png",7,10)),
  tar_file(lpsn_stats_file, "lpsn30.json"),
  tar_target(lpsn_stats, parse_lpsn_stats(lpsn_stats_file)),
  tar_target(fig_lpsn_stats, plot_lpsn_stats(lpsn_stats)),
  tar_file(fig_lpsn_stats_png, save_png(fig_lpsn_stats, "figures/fig_lpsn_stats.png", 6, 9)),
  tar_file(computing_metrics_parquet, "pocpbenchmark_computing_metrics.parquet"),
  tar_parquet(computing_metrics, read_parquet(computing_metrics_parquet)),
  tar_parquet(computing_metrics_fullbenchmark,
             computing_metrics %>% left_join(
               select(family_metadata, Family, benchmark_type), by = "Family"
             ) %>% filter(benchmark_type == "full") %>% select(-benchmark_type)),
  tar_target(db_parsed, get_db_parsed_stats(computing_metrics_fullbenchmark)),
  tar_target(tool_parsed, get_tool_parsed_stats(computing_metrics_fullbenchmark)),
  tar_target(median_db, generate_table_db(db_parsed)),
  tar_target(median_tool, generate_table_tool(tool_parsed)),
  tar_target(db_table, format_db_table(median_db)),
  tar_target(tool_table, format_tool_table(median_tool)),
  tar_target(p_pocp, plot_pocp_density(all_pocp)),
  tar_target(p_pocpu, plot_pocpu_density(all_pocpu)),
  tar_target(p_mcc, plot_mcc(mcc_pocpu_family, family_label, mcc_pocpu_global)),
  tar_target(p_mcc_random, plot_mcc_random(mcc_pocpu_family, family_label, mcc_pocpu_global)),
  tar_target(p_genus_delineation, plot_genus_delineation(p_pocp,p_pocpu, p_mcc)),
  tar_target(pocp_confusion, count(all_pocp, class) %>% deframe() %>% prettyNum(big.mark =" ")),
  tar_target(pocpu_confusion, count(all_pocpu, class) %>% deframe() %>% prettyNum(big.mark =" ")),
  tar_target(pocpu_confusion_by_family, count(all_pocpu, Family, class)),
  tar_target(lactobacillaceae,
             get_family_confusion_matrix(pocpu_confusion_by_family, "f__Lactobacillaceae"),
             ),
  tar_target(streptomycetaceae,
             get_family_confusion_matrix(pocpu_confusion_by_family, "f__Streptomycetaceae"),
             ),
  tar_target(fig_delta_genome_pocpu, plot_pocp_delta(all_pocpu, "POCPu", delta_genome, "Difference in genome size")),
  tar_target(fig_delta_proteome_pocpu, plot_pocp_delta(all_pocpu, "POCPu", delta_proteome, "Difference in proteome size")),
  tar_quarto(manuscript, quiet = FALSE)
)
