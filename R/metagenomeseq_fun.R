#' MEtagenomeSeq Analyse Diff
#'
#'
#' @param dada_res output from dada2_fun
#'
#' @return Return raw otu table in phyloseq object.
#' @import phyloseq
#' @import ggplot2
#' @import metagenomeSeq
#'
#' @export


# Decontam Function

metagenomeseq_fun <- function(data = data, output = "./metagenomeseq/", column1 = "", verbose = 1, rank = "Species", comp = ""){


  if(!dir.exists(output)){
    dir.create(output)
  }
  flog.info('Done.')

  if(rank != 'ASV'){
    data.glom <- tax_glom(data, taxrank=rank)
  } else {
    data.glom <- data
  }
  if(comp == ''){
    fun <- paste('combinaisons <- combn(na.omit(unique(sample_data(data.glom)$',column1,')),2) ',sep='')
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
  # save.image("debug.rdata")
  # quit()

  MGdata <- phyloseq_to_metagenomeSeq(data)
  '%!in%' <- function(x,y)!('%in%'(x,y))

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
    samplesToKeep = which(pData(MGdata)[,column1] == combinaisons[1,col] | pData(MGdata)[,column1] == combinaisons[2,col])
    obj_f = MGdata[featuresToKeep, samplesToKeep]

    #FitFeature model : zero-inflated log-normal model
    flog.info('Fitzig Model')
    pd <- pData(obj_f)
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
    TAB=cbind(row.names(TAB),TAB)
    colnames(TAB)[1]=paste(combinaisons[,col],collapse="_vs_")
    write.table(na.omit(TAB), paste(output,'/signtab_',column1,'_',paste(combinaisons[,col],collapse="_vs_"),'.csv',sep=''), row.names=FALSE, quote=FALSE, sep="\t")

  }


}
