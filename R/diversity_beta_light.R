#' Diversity Beta
#'
#' Provides ordination (PCoA, NMDS) of beta diversity indices (BrayCurtis, UniFrac and WeightedUnifrac) and statistical tests like PERMANOVA and pairwisePERMANOVA.
#'
#' @param data output from decontam or generate_phyloseq
#' @param output Output directory
#' @param glom Taxonomic rank to agglomerate data (one of rank_names(data) )
#' @param column1 Column name of main factor to test
#' @param covar One or more factor name to integrate as covariable in permanova. (for multiple covariable, provide as vector)
#' @param column2 Column name to split dataset with.
#' @param supp Supplementary plot (jaccard indexes + CCA/RDA ordination)
#'
#' @return Export plots and tests in the output directory.
#'
#' @import phyloseq
#' @import ggplot2
#' @import pairwiseAdonis
#' @import vegan
#' @importFrom plotly ggplotly
#'
#' @export


# Decontam Function

diversity_beta_light <- function(psobj, rank = "ASV", col = NULL, cov=NULL, dist0 = "bray", ord0 = "MDS", output="./plot_div_beta/", tests = TRUE) {

  if(!dir.exists(output)){
    dir.create(output, recursive=TRUE)
  }

  if( (dist0 == "unifrac" | dist0 == "wunifrac") & is.null(phy_tree(psobj, errorIfNULL=FALSE)) ){return(print("Error: unifrac distances needs phylogenetic tree."))}

  print(cov)
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
  otable = otu_table(data_rank)

  flog.info(glue::glue('Tests on {dist0} ...'))
  mdata = data.frame(sample_data(data_rank))
  mdata$Depth <- sample_sums(data_rank)
  print(table(mdata[,col1]))
  dist1 <<- vegdist(t(otable), distance=dist0)
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


  # Figure
  p1 <- plot_samples(data_rank, ordinate(data_rank, ord0, dist0), color = col ) + theme_bw() + ggtitle(glue::glue("{ord0} + {dist0}")) + stat_ellipse()
  # plot(p1)


  if(tests){
    ### Print tests
    sink(paste(output,'/',col,'_permANOVA.txt',sep=''), split = TRUE)
    cat("\n#####################\n##PERMANOVA on BrayCurtis distances\n#####################\n")
    print(resBC)
    cat("\n#####################\n##pairwisePERMANOVA on BrayCurtis distances\n#####################\n")
    print(resBC2)
    sink()

    facts=unlist(strsplit(col,"[+]"))
    if(length(facts)==2){
      fact1 = paste(mdata[,facts[1]],mdata[,facts[2]],sep="_")
      mdata$fact1 <- fact1
      col <- "fact1"
    }
  }



  return(p1)
}
