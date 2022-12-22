#' Heatmap
#'
#' Generate an heatmap with top taxa.
#'
#' @param data a phyloseq object (output from decontam or generate_phyloseq)
#' @param output Output directory
#' @param column1 Column name of factor to plot with (among sample_variables(data)).
#' @param top If not NULL, only this number of top features are plotted. Non "top" taxa abundances are aggregated in a new taxa named "Other".
#' @param norm Normalization method ("TSS" or "VST"), needs a phyloseq object with raw abundance.
#' @param rank Taxonomy rank to merge features that have same taxonomy at a certain taxonomic rank (among rank_names(data), or 'ASV' for no glom)
#' @param clust If TRUE, taxa are reordered with standard clustering. 
#' @param legend Legend title ("Abundance")
#'
#' @return Generate an heatmap with top taxa (html plotly version in output directory)
#' @importFrom htmlwidgets saveWidget
#'
#' @export


heatmap_fun <- function(data = data, column1 = "", top = NULL, output = "./plot_heatmap/", 
  rank = "Species", norm = "TSS", legend = "Abundance", clust = TRUE){
  LL = list()

  output1 <- paste(getwd(),'/',output,'/',sep='')
  if(!dir.exists(output1)){
    dir.create(output1, recursive = TRUE)
  }
  flog.info('Done.')

  if( is.na(match(column1, sample_variables(data))) ){stop("This factor is not in sample_variables(data).")}

  ps.glom.rel <- psobj.top <- dataglom <- data

  if(!is.null(rank)){
    dataglom <- tax_glom(data, rank)
    tt <- tax_table(dataglom)
    taxa <- tt[,rank]
    taxa_names(dataglom) <- taxa

    if(!is.null(top)){dataglom <- aggregate_top_taxa(dataglom, rank, top = top)}
  }


  if(norm == "TSS"){
    normf = function(x){ x/sum(x) }
    dataglom <- transform_sample_counts(dataglom, normf)
  }

  if(norm == "CLR"){ 
    clr = function(x){log(x+1) - rowMeans(log(x+1))}
    otable <- otu_table(dataglom)
    otableCLR <- clr(otable)
    otu_table(dataglom) <- otableCLR
  }

  #VST deseq2
  if(norm == "VST"){
      otable <- dataglom@otu_table@.Data+1
      otableVST <- DESeq2::varianceStabilizingTransformation(otable, fitType='local')
      dataglom@otu_table@.Data <- otableVST
  }


  otable <- otu_table(dataglom)
  sdata <- as.data.frame(as.matrix(sample_data(dataglom)))


  if(clust){
    h1 <- hclust(dist(otable))
    otable <- otable[h1$order,]
  }


  data.com <- reshape2::melt(otable)
  data.com$xlabel <- as.factor(sdata[as.character(data.com$Var2),match(column1,names(sdata))])
  names(data.com) <- c("Tax", "Sample", "Abundance", "xlabel")
  if(!clust){
    data.com$Tax = factor(data.com$Tax, levels = sort(unique(as.character(data.com$Tax))))
  }

<<<<<<< HEAD
=======



# PLOT

>>>>>>> 163e324275801571fc931402b3b6a463f9473cdb
  p.heat <- ggplot(data.com, aes(x = Sample, y = Tax)) + geom_tile(aes(fill = Abundance))
  p.heat <- p.heat + scale_fill_distiller(legend, palette = "RdYlBu") + theme_bw()

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

  flog.info(paste(output1,"heatmap_",column1,".html",sep=''))
  saveWidget(pltly.heat, file=  paste(output1,"heatmap_",column1,".html",sep=''))


  LL$plot <- p.heat
  LL$table <- data.com
  return(LL)
}
