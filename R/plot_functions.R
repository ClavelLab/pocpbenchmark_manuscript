# Plotting functions supporting {targets} analysis workflow

# POCP vs blast plot
plot_pocp_vs_blast <- function(df, pocp_label){
  # Get the min max values to set up matching x and y axes intervals
  extremes<- df %>%
    summarise(
      min = min(BLAST_BLASTP, pocp),
      max = max(BLAST_BLASTP, pocp)) %>%
    as_vector()
  
  p <- df %>%
    # Remove the database implementation for a supplementary figure
    filter(tool != "BLAST_BLASTPDB") %>%
    ggplot(aes(x = BLAST_BLASTP, y = pocp)) +
    ggdensity::geom_hdr_points(size=0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    coord_fixed()+
    scale_y_continuous(limits = extremes)+
    scale_x_continuous(limits = extremes)+
    facet_wrap(~ tool, nrow = 2) +
    labs(x = paste0(pocp_label, " based on BLASTP (in %)"),
         y = paste0(pocp_label, " based on other tools (in %)"),
         color = "Highest density\nregions probability")+
    theme_cowplot(font_size = 12)+
    theme(legend.position = "bottom", strip.text.x = element_text(size = 8))
  return(p)
}

# Subset tree of GTDB
plot_tree <- function(tree, metadata, size_factor=40){
  ggtree::ggtree(tree, layout = "circular") %<+%
    metadata+
    geom_tree(linewidth=0.01)+
    geom_tippoint(size = 0.3, mapping = aes(color = benchmark_type))+
    scale_color_manual(values = c("Recommended approach"="#CC79A7","All approaches"= "#009E73"))+
    theme_tree()+
    theme(legend.position = "none",
          panel.background = element_rect(fill = "transparent", color = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA),
          plot.margin =  grid::unit(c(-size_factor, -size_factor, -size_factor, -size_factor), "mm"))
}

# Counts of genomes per phyla and per benchmark type
plot_phyla_count<-function(tree_metadata){
  tree_metadata %>% count(benchmark_type,Phylum) %>% 
    mutate(Phylum = as_factor(Phylum),
           Phylum = fct_reorder(Phylum, n, .fun = function(x) max(sum(x)))) %>%
    ggplot(aes(x = Phylum, y = n, fill =benchmark_type))+
    scale_fill_manual(values = c("Recommended approach"="#CC79A7","All approaches"= "#009E73"))+
    geom_text(aes(label = n),hjust=-0.1,
              position = position_dodge2(preserve = "single", width = 0.8))+
    geom_col(position = position_dodge2(preserve = "single", width = 0.8))+
    labs(fill = "Benchmark type", x = "Phyla", y = "Number of genomes")+
    scale_y_continuous(expand = expansion(add = c(0,200)))+
    coord_flip()+
    theme_cowplot()
}

plot_lpsn_stats <- function(lpsn_stats){
  p_lpsn <- ggplot(lpsn_stats, aes(x = grouping, y = n))+
    geom_segment(
      aes(xend = grouping,x = grouping, y = n, yend=0),
      color = ifelse(lpsn_stats$grouping %in% c(2014,2024), "#E69F00", "#999999"),
      linewidth = ifelse(lpsn_stats$grouping %in% c(2014,2024), 1.5,1))+
    geom_point(
      color = ifelse(lpsn_stats$grouping %in% c(2014,2024), "#E69F00", "#999999"),
      linewidth = ifelse(lpsn_stats$grouping %in% c(2014,2024), 1.5,1)
    )+
    coord_flip()+
    labs(y="Validly published genus names under ICNP",
         x="Year")+
    theme_minimal_vgrid()+
    theme(axis.line.x =  element_line(color="black",linewidth = 0.5),
          axis.line.y = element_line(color="black",linewidth = 0.5))
  lpsn_annotations <- lpsn_stats %>%
    filter(grouping %in% c(2014,2024)) %>%
    select(grouping,n) %>% 
    mutate(label = glue::glue(
      "n = {n}\n  genera",
      n = prettyNum(n, big.mark = " ")
    )
    )
  
  p_lpsn+
    geom_text(data=lpsn_annotations,aes(x=grouping,y=n, label=label), 
              color="#E69F00", size=4 , angle=0, fontface="bold",
              hjust=-0.1, vjust=0.8
    )+
    scale_y_continuous(expand = expansion(add = c(0,800)))
}