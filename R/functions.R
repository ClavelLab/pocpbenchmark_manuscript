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
