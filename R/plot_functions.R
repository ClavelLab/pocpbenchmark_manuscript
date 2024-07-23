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