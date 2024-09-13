#' MEtaGenomeSeq Differential Analysis function.
#'
#'
#' @param data a phyloseq object (output from decontam or generate_phyloseq)
#' @param output Output directory
#' @param rank Taxonomy rank to merge features that have same taxonomy at a certain taxonomic rank (among rank_names(data), or 'ASV' for no glom)
#' @param column1 Column name of factor to test (among sample_variables(data))
#' @param verbose Verbose level. (1: quiet, 3: verbal)
#' @param comp Comma separated list of comparison to test. Comparisons are informed with a tilde (A~C,A~B,B~C). If empty, test all combination
#'
#'
#' @return Returns list with table of features and plots for each comparison. Export CSV files with significant differentialy abondant ASV.
#'
#' @import phyloseq
#' @import ggplot2
#' @import metagenomeSeq
#'
#' @export


# Decontam Function

metagenomeseq_fun <- function(data = data, output = "./metagenomeseq/", column1 = "",
                              verbose = 1, rank = "Species", comp = ""){

  if(verbose == 3){
    invisible(flog.threshold(DEBUG))
  } else {
    invisible(flog.threshold(INFO))
  }

  if(!dir.exists(output)){
    dir.create(output, recursive = TRUE)
  }
  ranks <- rank_names(data)
  if(comp == ''){
    fun <- paste('combinaisons <- combn(na.omit(unique(sample_data(data)$',column1,')),2) ',sep='')
    eval(parse(text=fun))
  }else{
    comp_list <- unlist(strsplit(comp,","))
    combinaisons <- matrix(, nrow = 2, ncol = length(comp_list))
    for (i in 1:length(comp_list)){
      tmp <- unlist(strsplit(comp_list[i],"\\~"))
      combinaisons[1,i] <- tmp[1]
      combinaisons[2,i] <- tmp[2]
    }
  }

if(rank != 'ASV'){
  data.glom <- tax_glom(data, taxrank=rank)
} else {
  data.glom <- data
}

  # save(list = ls(all.names = TRUE), file = "debug.rdata", envir = environment())
  MGdata <- phy2MGseq(data.glom) #phyloseq_to_metagenomeSeq

# featureData(MGdata)
# fData(MGdata)
  '%!in%' <- function(x,y)!('%in%'(x,y))

  outF = list()
  for (col in (1:ncol(combinaisons))){
    fun <- paste('tmp <- sample_data(data)$',column1,sep='')
    eval(parse(text=fun))
    if((combinaisons[1,col] %!in% tmp) || combinaisons[2,col] %!in% tmp){
      flog.warn(paste(combinaisons[1,col],' not in sample_data. Next;'),sep='')
      next
    }
    flog.info(paste('Combinaison ',combinaisons[1,col], ' ' , combinaisons[2,col],sep=''))
    fun <- paste('nb_cond1 <- nsamples(subset_samples(data.glom, ',column1, ' == "',combinaisons[1,col],'"))',sep='')
    eval(parse(text=fun))

    fun <- paste('nb_cond2 <- nsamples(subset_samples(data.glom, ',column1, ' == "',combinaisons[2,col],'"))',sep='')
    eval(parse(text=fun))
    flog.info(paste(combinaisons[1,col], 'has', nb_cond1, 'samples', '&' , combinaisons[2,col], 'has', nb_cond2, 'samples',sep=' '))

    flog.info(paste('Combinaison ',combinaisons[1,col], ' ' , combinaisons[2,col],sep=''))
    fun <- paste('tmp <- subset_samples(data.glom, ',column1, ' %in% c("',combinaisons[1,col],'","',combinaisons[2,col],'"))',sep='')
    eval(parse(text=fun))

    tmp <- prune_taxa(taxa_sums(tmp) >= 1, tmp)

    mdata <- sample_data(tmp)
    if(min(table(mdata[,column1])) == 1) {
      flog.info('Not enough samples')
      TAB <- data.frame(matrix(ncol=5))
      colnames(TAB) = c(paste(combinaisons[,col],collapse="_vs_"), "logFC", "se", "pvalues", "adjPvalues")
      write.table(na.omit(TAB), paste(output,'/signtab_',column1,'_',paste(combinaisons[,col],collapse="_vs_"),'.csv',sep=''), row.names=FALSE, quote=FALSE, sep="\t")
      next
    }

    featuresToKeep = which(rowSums(MGdata@assayData$counts) > 0)
    samplesToKeep = which(Biobase::pData(MGdata)[,column1] == combinaisons[1,col] | Biobase::pData(MGdata)[,column1] == combinaisons[2,col])
    obj_f = MGdata[featuresToKeep, samplesToKeep]

    #FitFeature model : zero-inflated log-normal model
    flog.info('Fitzig Model')
    pd <- Biobase::pData(obj_f)
    mod <- model.matrix(as.formula(paste("~", column1)), data = pd)

    res1 = NULL
    tryCatch( {res1 = fitFeatureModel(obj_f, mod)} ,
              error=function(e){e;cat("ERROR :",conditionMessage(e), "\n")})

    if(is.null(res1)) {flog.info(paste("Error: no table for condition test", paste(combinaisons[,col],collapse="_vs_"))); next}

    ( top = MRcoefs(res1, number=30) )
    diffTaxa <- row.names(top)
    flog.info('Save significant features (fdr padj < 0.05)')
    TAB = MRcoefs(res1) #fdr adjustment
    TAB = TAB[TAB$adjPvalues<=0.05,]

    nrow(TAB)
    TAB=cbind(row.names(TAB) ,TAB)

    TABF = na.omit(TAB)
    TABF <- cbind(TABF, Biobase::fData(MGdata)[row.names(TABF),1:match(rank,ranks)+1])
    colnames(TABF)[1]=paste(combinaisons[,col],collapse="_vs_")

    TABF <- TABF[order(TABF$logFC, decreasing = FALSE), ]
    flog.debug(head(TABF))
    write.table(na.omit(TABF), paste(output,'/signtab_',column1,'_',paste(combinaisons[,col],collapse="_vs_"),'.csv',sep=''), row.names=FALSE, quote=FALSE, sep="\t")

    #Plots
    fun <- paste('p <- ggplot(TABF,
    aes(x=reorder(',rank,', -logFC),
    y=logFC,
    color=',rank,')) + geom_point(size=6) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))',sep='')
    eval(parse(text=fun))

    outF[[paste(combinaisons[,col],collapse="_vs_")]] = list(plot = p, table = TABF)
  }
  return(outF)
}
