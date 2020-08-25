#' Plot beta diversity
#'
#' @param psobj a phyloseq object (output from decontam or generate_phyloseq)
#' @param rank Taxonomy rank to merge features that have same taxonomy at a certain taxonomic rank (among rank_names(data), or 'ASV' for no glom)
#' @param col A column name from metadata (among sample_variables(data))
#' @param cov One or more column names from metadata to treat as a covariable (among sample_variables(data)), provided as comma separated vector.
#' @param path Output directory
#' @param var Name of analysis for output customization


plot_beta <- function(psobj, rank, col, path, var = 'total', cov=covar) {
  cov1 = unlist(strsplit(cov, ","))
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

  flog.info('Bray ...')
  mdata = data.frame(sample_data(data_rank))
  mdata$Depth <- sample_sums(data_rank)
  print(table(mdata[,col1]))
  BC.dist <<- vegdist(t(otable), distance="bray")
  resBC = adonis(as.formula(paste('BC.dist ~ Depth +', paste(cov1, collapse="+"), "+", col)), data = mdata, permutations = 1000)
  #PairwiseAdonis
  # resBC2 = pairwise.adonis2(as.formula( paste('BC.dist ~ ', col,sep="") ), data = mdata)
  if(length(col1)>1){
    fact1 <- apply( mdata[,c(col1)] , 1 , paste , collapse = "-" )
    resBC2 = pairwise.adonis(BC.dist, fact1, p.adjust.m='fdr')
    resBC2 = pairwise.adonis(BC.dist, mdata[,c(col1)], p.adjust.m='fdr')
  }
  flog.info('Done')

  if(!is.null(phy_tree(psobj, errorIfNULL=FALSE))){
    flog.info('Unifrac ...')
    #UF wUF sur otus
    UF.dist <<- phyloseq::distance(psobj, "unifrac")
    resUF = adonis( as.formula(paste('UF.dist ~ Depth +', paste(cov1, collapse="+"), "+" ,col)), data = mdata, permutations = 1000)
    if(length(col1)>1){
      fact1 <- apply( mdata[,c(col1)] , 1 , paste , collapse = "-" )
      resUF2 = pairwise.adonis(UF.dist, fact1, p.adjust.m='fdr')
    } else {
      resUF2 = pairwise.adonis(UF.dist, mdata[,c(col1)], p.adjust.m='fdr')
    }
    flog.info('Done')

    flog.info('wunifrac ...')
    wUF.dist <<- phyloseq::distance(psobj, "wunifrac")
    reswUF = adonis( as.formula(paste("wUF.dist ~ Depth +", paste(cov1, collapse="+"), "+" , col)), data = mdata, permutations = 1000)

    if(length(col1)>1){
      fact1 <- apply( mdata[,c(col1)] , 1 , paste , collapse = "-" )
      reswUF2 = pairwise.adonis(wUF.dist, fact1, p.adjust.m='fdr')
    } else {
      reswUF2 = pairwise.adonis(wUF.dist, mdata[,c(col1)], p.adjust.m='fdr')
    }
    flog.info('Done')
  }
  sink(paste(path,'/',var,'_permANOVA.txt',sep=''), split = TRUE)
  cat("\n#####################\n##PERMANOVA on BrayCurtis distances\n#####################\n")
  print(resBC)
  cat("\n#####################\n##pairwisePERMANOVA on BrayCurtis distances\n#####################\n")
  print(resBC2)
  if(!is.null(phy_tree(data, errorIfNULL=FALSE))){
    cat("\n#####################\n##PERMANOVA on UniFrac distances\n#####################\n")
    print(resUF)
    cat("\n#####################\n##pairwisePERMANOVA on UniFrac distances\n#####################\n")
    print(resUF2)
    cat("\n#####################\n##PERMANOVA on Weighted UniFrac distances\n#####################\n")
    print(reswUF)
    cat("\n#####################\n##pairwisePERMANOVA on Weighted UniFrac distances\n#####################\n")
    print(reswUF2)
  }
  sink()


  facts=unlist(strsplit(col,"[+]"))
  if(length(facts)==2){
    fact1 = paste(mdata[,facts[1]],mdata[,facts[2]],sep="_")
    mdata$fact1 <- fact1
    col <- "fact1"
  }

  flog.info('Plotting ...')
  sample_data(data_rank) <- mdata
  p1 <- plot_samples(data_rank, ordinate(data_rank, "MDS", "bray"), color = col ) + theme_bw() + ggtitle(paste("MDS + BC")) + stat_ellipse()
  p2 <- plot_samples(data_rank, ordinate(data_rank, "NMDS", "bray"), color = col ) + theme_bw() + ggtitle(paste("NMDS + BC")) + stat_ellipse()
  p1bis <- plot_samples(data_rank, ordinate(data_rank, "MDS", "jaccard"), color = col ) + theme_bw() + ggtitle(paste("MDS + jaccard")) + stat_ellipse()
  p2bis <- plot_samples(data_rank, ordinate(data_rank, "NMDS", "jaccard"), color = col ) + theme_bw() + ggtitle(paste("NMDS + jaccard")) + stat_ellipse()
  # mdata1 = sample_data(psobj)
  sample_data(psobj) <- mdata
  if(!is.null(phy_tree(data, errorIfNULL=FALSE))){
    p3 <- plot_samples(psobj, ordinate(psobj, "MDS", "unifrac"), color = col) + theme_bw() + ggtitle(paste("MDS + UF")) + stat_ellipse()
    p4 <- plot_samples(psobj, ordinate(psobj, "NMDS", "unifrac"), color = col) + theme_bw() + ggtitle(paste("NMDS + UF")) + stat_ellipse()
    p5 <- plot_samples(psobj, ordinate(psobj, "MDS", "wunifrac"),color = col) + theme_bw() + ggtitle(paste("MDS + wUF")) + stat_ellipse()
    p6 <- plot_samples(psobj, ordinate(psobj, "NMDS","wunifrac"), color = col) + theme_bw() + ggtitle(paste("NMDS + wUF"))+ stat_ellipse()
    # save(p1,p2,p3,p4,p5,p6, file = "diversity_plots.rdata")
    flog.info('Done.')

    flog.info('Saving ...')
    if(var!=''){
      png(paste(path,'/',var,"_beta.png",sep=''), width=40,height=50, units="cm", res=200)
    }else{
      png(paste(path,"/beta.png",sep=''), width=40,height=50, units="cm", res=200)
    }
    ppp = grid.arrange(p1 +  theme(legend.position = "none"), # + stat_ellipse()
                       p2 + theme(legend.position = "none" ),
                       p1bis + theme(legend.position = "none" ),
                       p2bis + theme(legend.position = "none" ),
                       p3 + theme(legend.position = "none"),
                       p4 + theme(legend.position = "none"),
                       p5 + theme(legend.position = "none"),
                       p6 + theme(legend.position = "none"), ncol = 2)
    dev.off()
    betaRes$plot <- ppp
  }else{
    flog.info('No phy_tree ...')
    ppp = grid.arrange(p1 +  theme(legend.position = "none"), #
                       p2 +  theme(legend.position = "none"),
                       p1bis +  theme(legend.position = "none"),
                       p2bis +  theme(legend.position = "none"), ncol = 2)
    betaRes$plot <- ppp
  }

  flog.info('Supplement Beta plots ...')
  if(supp == TRUE){
    p1 <- plot_samples(data_rank, ordinate(data_rank, "CCA", "bray"), color = col ) + theme_bw() + ggtitle(paste("CCA + BC")) + stat_ellipse()
    p2 <- plot_samples(data_rank, ordinate(data_rank, "RDA", "bray"), color = col ) + theme_bw() + ggtitle(paste("RDA + BC")) + stat_ellipse()
    p1bis <- plot_samples(data_rank, ordinate(data_rank, "CCA", "jaccard"), color = col ) + theme_bw() + ggtitle(paste("CCA + jaccard")) + stat_ellipse()
    p2bis <- plot_samples(data_rank, ordinate(data_rank, "RDA", "jaccard"), color = col ) + theme_bw() + ggtitle(paste("RDA + jaccard")) + stat_ellipse()
    # mdata1 = sample_data(psobj)
    sample_data(psobj) <- mdata
    if(!is.null(phy_tree(data, errorIfNULL=FALSE))){
      p3 <- plot_samples(psobj, ordinate(psobj, "CCA", "unifrac"), color = col) + theme_bw() + ggtitle(paste("CCA + UF")) + stat_ellipse()
      p4 <- plot_samples(psobj, ordinate(psobj, "RDA", "unifrac"), color = col) + theme_bw() + ggtitle(paste("RDA + UF")) + stat_ellipse()
      p5 <- plot_samples(psobj, ordinate(psobj, "CCA", "wunifrac"),color = col) + theme_bw() + ggtitle(paste("CCA + wUF")) + stat_ellipse()
      p6 <- plot_samples(psobj, ordinate(psobj, "RDA","wunifrac"), color = col) + theme_bw() + ggtitle(paste("RDA + wUF"))+ stat_ellipse()
      # save(p1,p2,p3,p4,p5,p6, file = "diversity_plots.rdata")
      flog.info('Done.')

      flog.info('Saving ...')
      if(var!=''){
        png(paste(path,'/',var,"_betasupp.png",sep=''), width=40,height=50, units="cm", res=200)
      }else{
        png(paste(path,"/betasupp.png",sep=''), width=40,height=50, units="cm", res=200)
      }
      ppp = grid.arrange(p1 +  theme(legend.position = "none"), # + stat_ellipse()
                         p2 + theme(legend.position = "none" ),
                         p1bis + theme(legend.position = "none" ),
                         p2bis + theme(legend.position = "none" ),
                         p3 + theme(legend.position = "none"),
                         p4 + theme(legend.position = "none"),
                         p5 + theme(legend.position = "none"),
                         p6 + theme(legend.position = "none"), ncol = 2)

      dev.off()
      betaRes$plot <- ppp

    }else{
      flog.info('No phy_tree ...')
      ppp = grid.arrange(p1 +  theme(legend.position = "none"), #
                         p2 +  theme(legend.position = "none"),
                         p1bis +  theme(legend.position = "none"),
                         p2bis +  theme(legend.position = "none"), ncol = 2)
      betaRes$plot <- ppp
    }
  }
  flog.info('Done.')
}



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
#'
#' @export


# Decontam Function

diversity_beta_fun <- function(data = data, output = "./plot_div_beta/", glom = "ASV", column1 = "", column2 = "", covar ="", supp = FALSE){

  # suppressMessages(source(system.file("supdata", "phyloseq_extended_graphical_methods.R", package="ranomaly")))
  if(!dir.exists(output)){
    dir.create(output, recursive=TRUE)
  }

  mdata <- as.matrix(sample_data(data))
  fact1 <- paste(na.omit(mdata[,column1]))
  # print(levels(as.factor(na.omit(mdata[,column1]))))
  #too few levels

  if(length(levels(as.factor(fact1))) > 1 ){
    # If too few observations.
    if(min(table(fact1)) < 3 & column2 != ""){
      p <- plot_beta(data, glom, paste(column1, column2, sep="+") , output, var="global")
      p <- plot_beta(data, glom, column1, output, var=column1)
      p <- plot_beta(data, glom, column2, output, var=column2)
    }else{
      if(column2 != ""){
        flog.info('Option1...')
        fact1 <- paste(mdata[,column1], mdata[,column2],sep="_")
        vector <- levels(data.frame(sample_data(data)[,column1])[,1])
        print(vector)
        for (var in vector){
          flog.info('Split table %s...', var)
          fun <- paste('psobj <- subset_samples(data, ',column1,'=="',var,'")',sep='')
          eval(parse(text=fun))
          if(all(is.na(sample_data(psobj)[,column2]))){
            flog.info(paste("For",column1,"=", var, column2,'is all NA...')); next
          }
          psobj <- prune_taxa(taxa_sums(psobj) >= 1, psobj)
          flog.info('Done.')
          p <- plot_beta(psobj, glom, column2, output, var)

          ## Global permanova + ordination
        }
        flog.info('Global1...')
        p <- plot_beta(data, glom, paste(column1, column2, sep="+") , output, var="global")
        flog.info('Global2...')
        p <- plot_beta(data, glom, column1, output, var=column1)
        flog.info('Global3...')
        p <- plot_beta(data, glom, column2, output, var=column2)
      }else{
        flog.info('Option2...')
        p <- plot_beta(data, glom, column1, output)
      }

    }

  }else{flog.info('Too few levels < 2.'); print(levels(as.factor(fact1)))
    if(column2 != "" & !all(is.na(mdata[,column2]))){
      print(sample_data(data))
      fact1=paste(mdata[,column1],mdata[,column2], sep="_")
      mdata=cbind(mdata,fact1)
      sample_data(data)=mdata
      p <- plot_beta(data, glom, fact1, output, var=column2)
    } else{flog.info('No test.')}
  }


  save.image(paste(output,"/beta_env_bak.rdata", sep=""))
  flog.info('Finish')
  return(betaRes)


}
