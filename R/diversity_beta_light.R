#' Diversity Beta
#'
#' Provides ordination (PCoA, NMDS) of beta diversity indices (BrayCurtis, UniFrac and WeightedUnifrac) and statistical tests like PERMANOVA and pairwisePERMANOVA.
#'
#' @param psobj a phyloseq object (output from decontam or generate_phyloseq)
#' @param rank Taxonomy rank to merge features that have same taxonomy at a certain taxonomic rank (among rank_names(data), or 'ASV' for no glom)
#' @param col A metadata column (among sample_variables(data)).
#' @param cov Covariable names comma separated vector.
#' @param dist0 Dissimilarity index, partial match to "unifrac", "wunifrac", "dpcoa", "jsd", "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao" or "mahalanobis".
#' @param ord0 Currently supported method options are: c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA")
#' @param output The output file directory.
#' @param tests Whether to compute tests or not (TRUE/FALSE)
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

diversity_beta_light <- function(psobj, rank = "ASV", col = NULL, cov = NULL, dist0 = "bray", ord0 = "MDS", output="./plot_div_beta/", tests = TRUE) {

  if(!dir.exists(output)){
    dir.create(output, recursive=TRUE)
  }

  if( (dist0 == "unifrac" | dist0 == "wunifrac") & is.null(phy_tree(psobj, errorIfNULL=FALSE)) ){return(print("Error: unifrac distances needs phylogenetic tree."))}

  flog.debug(cov)
  if(!is.null(cov)){
    cov1 = unlist(strsplit(cov, ","))
  }

  # metrics <- sapply(strsplit(measures,","), '[')
  col1 = unlist(strsplit(col, "[+]"))
  if(length(col1)==1){
    fun <- paste("psobj = subset_samples(psobj, !is.na(",col,"))",sep="")
    eval(parse(text=fun))
  }else{
    fun <- paste("psobj = subset_samples(psobj, !is.na(",col1[1],"))",sep="")
    eval(parse(text=fun))
    fun <- paste("psobj = subset_samples(psobj, !is.na(",col1[2],"))",sep="")
    eval(parse(text=fun))
  }
  if(rank=="ASV"){
    flog.info('No glom ...')
    data_rank = psobj
  }else{
    data_rank = tax_glom(psobj, rank)
  }


  # Figure
  resBeta = list()
  if(!is.null(cov)){
  p1 <- plot_samples(data_rank, ordinate(data_rank, ord0, dist0), color = col, shape = cov1[1] ) +
  theme_bw() + ggtitle(glue::glue("{ord0} + {dist0}")) + stat_ellipse()
}else{
  p1 <- plot_samples(data_rank, ordinate(data_rank, ord0, dist0), color = col) +
  theme_bw() + ggtitle(glue::glue("{ord0} + {dist0}")) + stat_ellipse()
}
  # plot(p1)

  resBeta$plot = p1 + theme(axis.text.x = element_text(angle = 45, hjust=1),
  legend.position = "none",axis.text=element_text(size=18),
  axis.title=element_text(size=16,face="bold"),
  strip.text.x = element_text(size = 18,face="bold"),
  title=element_text(size=16,face="bold"))

  if(tests){
    otable = otu_table(data_rank)

    flog.info(glue::glue('Tests on {dist0} ...'))
    mdata = data.frame(sample_data(data_rank))
    mdata$Depth <- sample_sums(data_rank)
    if( any(grepl(dist0, c("unifrac", "wunifrac", "dpcoa", "jsd") )) ){
      dist1 <-phyloseq::distance(data_rank, dist0)
    }else{
      dist1 <<- vegdist(t(otable), method = dist0)
    }
    if(!is.null(cov)){
      resBC = adonis(as.formula(paste('dist1 ~ Depth +', paste(cov1, collapse="+"), "+", col)), data = mdata, permutations = 1000)
    }else{
      resBC = adonis(as.formula(paste('dist1 ~ Depth +', col)), data = mdata, permutations = 1000)
    }

    #PairwiseAdonis
    # resBC2 = pairwise.adonis2(as.formula( paste('BC.dist ~ ', col,sep="") ), data = mdata)
    if(length(col1)>1){
      fact1 <- apply( mdata[,c(col1)] , 1 , paste , collapse = "-" )
      resBC2 = pairwise.adonis(dist1, fact1, p.adjust.m='fdr')
    } else {
      resBC2 = pairwise.adonis(dist1, mdata[,c(col1)], p.adjust.m='fdr')
    }

    write.table(resBC$aov.tab, file=paste0(output,'/',col,'_permANOVA.txt'), sep="\t")
    write.table(resBC2, file=paste0(output,'/',col,'pairwisepermANOVA.txt'), sep="\t")

    ggsave(glue::glue("{output}/beta_diversity.eps"), plot=resBeta$plot, height = 15, width = 30, units="cm", dpi = 500, device="eps")

    resBeta$permanova = resBC$aov.tab
    resBeta$pairwisepermanova = resBC2

  }
  return(resBeta)
}
