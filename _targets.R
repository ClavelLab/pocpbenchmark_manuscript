library(targets)
library(tarchetypes)

tar_option_set(
  packages = c("tidyverse",
                "ggplot2", "cowplot", "ggokabeito"),
  
)

tar_source()
# Set target-specific options such as packages:
# tar_option_set(packages = "utils") # nolint

# End this file with a list of target objects.
list(
  tar_files_input(archives,list.files("data_benchmark", full.names = T)),
  tar_target(extracted, unzip(archives), pattern = map(archives), format = "file"),
  tar_files(prots, list.files(pattern = "benchmark-gtdb-f*"), format = "file"),
  tar_target(prot_stats, read_stats(prots), pattern = map(prots),iteration = "vector"),
  tar_target(prot_n, prot_stats %>% dplyr::select(family, num_seqs)%>%dplyr::group_by(family) %>%
               summarise(n_genomes = n(), median_proteins = median(num_seqs))),
  tar_target(table_overview, readr::read_csv("cpuhours.csv")%>% left_join(prot_n, by="family")),
  tar_target(r2, read_R2(prots), pattern = map(prots),iteration = "vector"),
  
  tar_target(pocp_values, read_pocp(prots), pattern = map(prots)),
  tar_target(plot_pocp, plot_pocp_distribution(pocp_values, "POCP")),
  tar_target(plot_pocpu, plot_pocp_distribution(pocp_values, "POCPu"))
)
