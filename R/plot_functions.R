# Plotting functions supporting {targets} analysis workflow


# Save ggplot2 plots for file, a wrapper around cowplot::ggsave2 w/ sensible defaults
save_png <- function(plot, filename, width, height){
  fs::dir_create("figures")
  cowplot::ggsave2(filename = filename,
                   plot = plot, width = width, height = height,
                   units = "in", dpi = 300, device = grDevices::png)
}
# POCP vs blast plot
plot_pocp_vs_blast <- function(df, pocp_label, R2_table){
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
    labs(x = paste0(pocp_label, " based on BLAST_BLASTP (in %)"),
         y = paste0(pocp_label, " based on other tools (in %)"),
         color = "Highest density\nregions")+
    theme_cowplot(font_size = 12)+
    theme(legend.position = "bottom", strip.text.x = element_text(size = 8))
  
  
  df_R2 <- R2_table %>% filter(tool != "BLAST_BLASTPDB") %>% 
    select(tool, starts_with(paste0(pocp_label,"_"))) %>%
    rename_with(~str_remove(.x,paste0(pocp_label,"_"))) %>%
    mutate(
      label = glue::glue("atop(italic(R)^2 == {R2}, italic(p){p_label})",
                         R2=round(R2,digits = 3) %>% prettyNum())
    )
  
  p <- p + 
    geom_text(data = df_R2,
              aes(x = extremes[1], y=0.75*extremes[2], label = label),
              size=3, hjust=0, vjust=0,parse = TRUE)
  return(p)
}

# POCP vs blast plot
plot_pocp_blastdb <- function(df, pocp_label, with_R2=TRUE){
  # 
  df <- df %>% filter(tool == "BLAST_BLASTPDB")
  # Get the min max values to set up matching x and y axes intervals
  extremes<- df %>%
    summarise(
      min = min(BLAST_BLASTP, pocp),
      max = max(BLAST_BLASTP, pocp)) %>%
    as_vector()
  
  p <- df %>%
    ggplot(aes(x = BLAST_BLASTP, y = pocp)) +
    ggdensity::geom_hdr_points(size=0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    coord_fixed()+
    scale_y_continuous(limits = extremes)+
    scale_x_continuous(limits = extremes)+
    labs(x = paste0(pocp_label, " based on BLAST_BLASTP (in %)"),
         y = paste0(pocp_label, " based on BLAST_BLASTPDB (in %)"),
         color = "Highest density\nregions")+
    theme_cowplot(font_size = 12)+
    theme(legend.position = "bottom", strip.text.x = element_text(size = 8))
  
  if(with_R2){
    df_R2 <- df %>% group_by(tool) %>%
      nest() %>%
      mutate(
        fit = map(data, ~ lm(BLAST_BLASTP ~ pocp, data = .)),
        tidied = map(fit, broom::glance)
      ) %>%
      unnest(tidied) %>%
      select(-data, -fit) %>%
      ungroup() %>% select(tool,r.squared,p.value,nobs ) %>% 
      rename("R2"="r.squared", "p"="p.value") %>% 
      mutate(
        p_label = map_chr(p, scales::label_pvalue()),
        label = glue::glue("atop(italic(R)^2 == {R2}, italic(p){p_label})",
                           R2=round(R2,digits = 3) %>% prettyNum()))
    
    
    p <- p + 
      geom_text(data = df_R2,
                aes(x = extremes[1], y=0.75*extremes[2], label = label),
                size=3, hjust=0, vjust=0,parse = TRUE)
  }
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

plot_tree_phyla_count <- function(tree,phyla_count,genome_metadata){
  library(ggtree)
  leg<-cowplot::get_legend(phyla_count)
  plot_grid(
    ggtree::rotate_tree(tree, 20),
    phyla_count+theme(legend.position = "none"),
    nrow = 2, labels = "AUTO",
    rel_heights = c(0.7,0.3)
  )+draw_plot(leg, x=0.65,y=-0.15)
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

plot_pocp_density <- function(df){
  df %>% 
    ggplot(aes(x = pocp, fill = same_genus_truth))+
    geom_density(alpha = 0.5)+
    theme_cowplot(font_size = 12, rel_small = 10/14)+
    scale_fill_okabe_ito(labels = c("TRUE"="Within genus","FALSE"="Between genera"))+
    geom_vline(xintercept = 50, linetype="dashed")+
    labs(x="POCP", fill = "True category", y = "Density")+
    scale_y_continuous(expand = expansion(mult  = c(0,0.6)))+
    scale_x_continuous(expand = expansion(add = c(25,0)))+
    annotate("segment", x = 55, xend = 70, y=0.029, yend = 0.029,
             color="black", arrow=arrow(type= "closed", length = unit(0.05, "in"))
    )+
    annotate("text",x = 55, y=0.029, hjust=0, vjust=-0.85,
             label="Same\ngenus")+
    annotate("segment", x = 45, xend = 30, y=0.029, yend = 0.029,
             color="black", arrow=arrow(type= "closed", length = unit(0.05, "in"))
    )+
    annotate("text",x = 45, y=0.029, hjust=1, vjust=-0.85,
             label="Different\ngenera")+
    theme(legend.position = "bottom")
}

plot_pocpu_density <- function(df){
  df %>% 
    ggplot(aes(x = pocp, fill = same_genus_truth))+
    geom_density(alpha = 0.5)+
    theme_cowplot(font_size = 12, rel_small = 10/14)+
    scale_fill_okabe_ito(labels = c("TRUE"="Within genus","FALSE"="Between genera"))+
    geom_vline(xintercept = 50, linetype="dashed")+
    labs(x="POCPu", fill = "True category",y = "Density")+
    theme(legend.position = "bottom")+
    scale_y_continuous(expand = expansion(mult = c(0,0.01)))
}

plot_mcc <- function(mcc_df_per_family,family_label, mcc_df_global){
  mcc_df_per_family <- mcc_df_per_family %>%
    left_join(family_label, by = "Family") %>% 
    mutate(label = as_factor(label),
           label = fct_reorder(label, mcc, .fun = max))
  
  p <- ggplot(mcc_df_per_family, aes(x = mcc, y = label)) +
    geom_vline(data = mcc_df_global, aes(xintercept = mcc),
               linetype = "dashed")+
    geom_segment(aes(xend = 0, yend=label),color = "#999999")+
    geom_point(color = "black")+
    scale_y_discrete(labels = function(x) as.character(x) %>% parse(text = .))+
    geom_text(aes(label = prettyNum(mcc, digits = 2)),
              hjust=-0.3, size=3)+
    theme_cowplot(font_size = 12,rel_small = 10/14)+
    theme(legend.position = "none")+
    labs(x = "Matthews Correlation Coefficient (MCC)", y ="Bacterial families" )+
    facet_grid(rows = vars(Phylum), scales = "free_y", space = "free",switch = "y",
               labeller = as_labeller(function(x) str_remove(x,"p__")))+
    scale_x_continuous(expand = expansion(mult =  c(0,0.1)), limits = c(0,1))
  rm(mcc_df_per_family)
  return(p)
}


plot_mcc_random <- function(mcc_df_per_family,family_label, mcc_df_global){
  mcc_df_per_family <- mcc_df_per_family %>%
    left_join(family_label, by = "Family") %>% 
    mutate(label = as_factor(label),
           label = fct_reorder(label, mcc_random, .fun = max))
  
  p <- ggplot(mcc_df_per_family, aes(x = mcc_random, y = label)) +
    geom_vline(data = mcc_df_global, aes(xintercept = mcc_random),
               linetype = "dashed")+
    geom_segment(aes(xend = 0, yend=label),color = "#999999")+
    geom_point(color = "black")+
    scale_y_discrete(labels = function(x) as.character(x) %>% parse(text = .))+
    geom_text(data = function(x) {x[x$mcc_random >= 0, ]},
              aes(label = prettyNum(mcc_random, digits = 2)),
              hjust=-0.3, size=3)+
    geom_text(data = function(x) {x[x$mcc_random < 0, ]},
              aes(label = prettyNum(mcc_random, digits = 2)),
              hjust=1.3, size=3)+
    theme_cowplot(font_size = 12,rel_small = 10/14)+
    theme(legend.position = "none")+
    labs(x = "Matthews Correlation Coefficient", y ="Bacterial families" )+
    facet_grid(rows = vars(Phylum), scales = "free_y", space = "free",switch = "y",
               labeller = as_labeller(function(x) str_remove(x,"p__")))+
    scale_x_continuous(expand = expansion(mult =  c(0,0.1)), limits = c(-0.25,1))
  rm(mcc_df_per_family)
  return(p)
}

plot_genus_delineation <- function(p_pocp_density, p_pocpu_density, p_mcc){
  p_densities <- plot_grid(
    p_pocp_density+theme(legend.position = "none"),
    p_pocpu_density,
    nrow = 2,
    labels = c("A","B")
  )
  plot_grid(p_densities, p_mcc,
            ncol = 2, axis="b",labels = c("","C"),
            rel_widths =  c(0.4, 0.6))
}

plot_pocp_delta <- function(df, pocp_label, delta, delta_label){
  ggplot(df, aes(x = pocp, y = {{ delta }}))+
    ggdensity::geom_hdr_points()+
    scale_y_continuous(labels = scales::label_number())+
    facet_wrap(~same_genus_truth,
               labeller = as_labeller( c("TRUE"="Within genus","FALSE"="Between genera")))+
    theme_cowplot()+
    labs(x = pocp_label ,  y = delta_label,
         color = "Highest density\nregions")
}