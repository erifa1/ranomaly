#' Diversity Beta
#'
#' Provides ordination (PCoA, NMDS) of beta diversity indices (BrayCurtis, UniFrac and WeightedUnifrac) and statistical tests like PERMANOVA and pairwisePERMANOVA.
#'
#' @param psobj a phyloseq object (output from decontam or generate_phyloseq)
#' @param rank Taxonomy rank to merge features that have same taxonomy at a certain taxonomic rank (among rank_names(data), or 'ASV' for no glom)
#' @param col A metadata column (among sample_variables(data)).
#' @param cov Covariable names comma separated vector (last covariable is used to plot samples with different shape).
#' @param dist0 Dissimilarity index, partial match to "unifrac", "wunifrac", "dpcoa", "jsd", "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao" or "mahalanobis".
#' @param ord0 Currently supported method options are: c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA")
#' @param output The output file directory.
#' @param tests Whether to compute tests or not (TRUE/FALSE)
#' @param axes Axes to plot (c(1,2))
#' @param ellipse Plot ellipse (TRUE)
#' @param verbose Verbose level. (1: quiet, 2: print infos, 3: print infos + debug)
#'
#' @return Return specific plots and tests in list and output them in the output directory.
#'
#' @import phyloseq
#' @import ggplot2
#' @import pairwiseAdonis
#' @import vegan
#' @importFrom plotly ggplotly
#'
#' @export


# Decontam Function

diversity_beta_light <- function(psobj, rank = "ASV", col = NULL, cov = NULL, dist0 = "bray", ord0 = "MDS", output="./plot_div_beta/", axes = c(1,2), tests = TRUE, verbose = 2, ellipse = TRUE) {

  if(verbose == 3){
  invisible(flog.threshold(DEBUG))
  }
  if(verbose == 2){
    invisible(flog.threshold(INFO))
  }
  if(verbose == 1){
    invisible(flog.threshold(ERROR))
  }


  if(!dir.exists(output)){
    dir.create(output, recursive=TRUE)
  }

  if( (dist0 == "unifrac" | dist0 == "wunifrac") & is.null(phy_tree(psobj, errorIfNULL=FALSE)) ){return(print("Error: unifrac distances needs phylogenetic tree."))}

  flog.debug(cov)
  if(!is.null(cov)){
    cov1 <- unlist(strsplit(cov, ","))
  }

  col1 <- col
  if(length(col1)==1){
    fun <- paste("psobj = subset_samples(psobj, !is.na(",col,"))",sep="")
    eval(parse(text=fun))
  }else{
    stop("Only one factor in 'col' argument, add more covariable in 'cov' argument.")
    }

  if(rank=="ASV"){
    flog.info('No glom ...')
    data_rank <- psobj
  }else{
    data_rank <- tax_glom(psobj, rank)
  }


  # Figure
  flog.info('Plot ...')
  resBeta <- list()
  if(!is.null(cov)){

  sdata = sample_data(data_rank)
  fun = glue::glue("sdata${cov1[length(cov1)]}_{col} = factor(paste(sdata${cov1[length(cov1)]}, sdata${col}, sep='_'))")
  eval(parse(text=fun))
  sample_data(data_rank) = sdata

  p1 <- plot_samples(data_rank, ordinate(data_rank, ord0, dist0), color = col, shape = cov1[length(cov1)], axes = axes ) +
  theme_bw() + ggtitle(glue::glue("{ord0} + {dist0}")) + scale_shape_manual(values = 0:10)
  if(ellipse){p1 <- p1 + stat_ellipse()}

  p2 <- plot_samples(data_rank, ordinate(data_rank, ord0, dist0), color = glue::glue("{cov1[length(cov1)]}_{col}"), shape = NULL, axes = axes ) +
  theme_bw() + ggtitle(glue::glue("{ord0} + {dist0}")) + scale_shape_manual(values = 0:10)
  if(ellipse){p2 <- p2 + stat_ellipse()}

  resBeta$plot2 <- p2 + theme(axis.text.x = element_text(angle = 45, hjust=1),
  ,axis.text=element_text(size=18),
  axis.title=element_text(size=16,face="bold"),
  strip.text.x = element_text(size = 18,face="bold"),
  title=element_text(size=16,face="bold"))
  ggsave(glue::glue("{output}/beta_diversity2.eps"), plot=resBeta$plot2, height = 20, width = 30, units="cm", dpi = 500, device="eps")


}else{
  p1 <- plot_samples(data_rank, ordinate(data_rank, ord0, dist0), color = col, axes = axes ) +
  theme_bw() + ggtitle(glue::glue("{ord0} + {dist0}"))
  if(ellipse){p1 <- p1 + stat_ellipse()}
}
  # plot(p1)
  flog.info('Plot ok...')
  resBeta$plot <- p1 + theme(axis.text.x = element_text(angle = 45, hjust=1),
  ,axis.text=element_text(size=18),
  axis.title=element_text(size=16,face="bold"),
  strip.text.x = element_text(size = 18,face="bold"),
  title=element_text(size=16,face="bold"))

  if(tests){
    otable <- otu_table(data_rank)

    flog.info(glue::glue('Tests on {dist0} ...'))
    mdata <- data.frame(sample_data(data_rank))
    mdata$Depth <- sample_sums(data_rank)
    if( any(grepl(dist0, c("unifrac", "wunifrac", "dpcoa", "jsd") )) ){
      dist1 <-phyloseq::distance(data_rank, dist0)
    }else{
      dist1 <<- vegdist(t(otable), method = dist0)
    }
    if(!is.null(cov)){
      form1 <- as.formula(paste('dist1 ~ Depth +', paste(cov1, collapse="+"), "+", col))
      resBC <- adonis(form1, data = mdata, permutations = 1000)
        
    }else{
      form1 <- as.formula(paste('dist1 ~ Depth +', col))
      resBC <- adonis(form1, data = mdata, permutations = 1000)
    }

    #PairwiseAdonis
    if(!is.null(cov)){
      fun = glue::glue("resBC2 <- pairwise.adonis(dist1, mdata${cov1[length(cov1)]}_{col}, p.adjust.m='fdr')" )
      eval(parse(text = fun))
    } else {
      resBC2 <- pairwise.adonis(dist1, mdata[,c(col1)], p.adjust.m='fdr')
    }

    write.table(resBC$aov.tab, file=paste0(output,'/',col,'_permANOVA.txt'), sep="\t")
    write.table(resBC2, file=paste0(output,'/',col,'pairwisepermANOVA.txt'), sep="\t")

    ggsave(glue::glue("{output}/beta_diversity.eps"), plot=resBeta$plot, height = 20, width = 30, units="cm", dpi = 500, device="eps")

    resBeta$permanova <- resBC$aov.tab
    resBeta$permanova_formula <- format(form1)
    resBeta$pairwisepermanova <- resBC2
    resBeta$test_table <- mdata
    resBeta$dist <- dist1

  }
  return(resBeta)
}