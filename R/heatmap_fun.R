#' Heatmap
#'
#' Generate an heatmap with top taxa.
#'
#' @param data output from decontam or generate_phyloseq
#' @param output Output directory
#' @param column1 Column name of factor to plot with.
#' @param top Number of top features to plot.
#' @param rank Taxonomic rank name.
#'
#' @return Generate an heatmap with top taxa (html plotly version in output directory)
#' @importFrom htmlwidgets saveWidget
#'
#' @export


# Decontam Function

heatmap_fun <- function(data = data, column1 = "", top = 20, output = "./plot_heatmap/", rank = "Species"){

  output1 <- paste(getwd(),'/',output,'/',sep='')
  if(!dir.exists(output1)){
    dir.create(output1)
  }
  flog.info('Done.')

  psobj.top <- microbiome::aggregate_top_taxa(data, rank, top = top)

  ps.glom.rel <- microbiome::transform(psobj.top, "compositional")

  plot.composition.relAbun <- plot_composition(ps.glom.rel, x.label = column1)
  data.com <- plot.composition.relAbun$data
  colnames(data.com)

  p.heat <- ggplot(data.com, aes(x = Sample, y = Tax)) + geom_tile(aes(fill = Abundance))
  p.heat <- p.heat + scale_fill_distiller("Abundance", palette = "RdYlBu") + theme_bw()

  # Make bacterial names italics
  p.heat <- p.heat + theme(axis.text.y = element_text(colour = 'black',
                                                      size = 10,
                                                      face = 'italic'))
  # Make seperate samples based on main varaible
  p.heat <- p.heat + facet_grid(~xlabel, scales = "free")

  p.heat <- p.heat + ylab(rank)

  #Clean the x-axis
  p.heat <- p.heat + theme(axis.title.x=element_blank(),
                           axis.text.x=element_text(angle = 90),
                           axis.ticks.x=element_blank())

  # Clean the facet label box
  p.heat <- p.heat + theme(legend.key = element_blank(),
                           strip.background = element_rect(colour="black", fill="white"))
  pltly.heat <- ggplotly(p.heat)

  print(paste(output1,"heatmap_",column1,".html",sep=''))
  saveWidget(pltly.heat, file=  paste(output1,"heatmap_",column1,".html",sep=''))

  return(p.heat)
}
