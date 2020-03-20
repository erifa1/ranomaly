#' DESEQ2 Analyse Diff
#'
#'
#' @param dada_res output from dada2_fun
#'
#' @return Return raw otu table in phyloseq object.
#' @import phyloseq
#' @import ggplot2
#' @import DESeq2
#' @import ggpubr
#' @import gridExtra
#' @import taxa
#'
#' @export


# Decontam Function

deseq2_fun <- function(data = data, output = "./deseq/", column1 = "", verbose = 1, rank = "Species", comp = ""){


  if(verbose == 3){
    invisible(flog.threshold(DEBUG))
  } else {
    invisible(flog.threshold(INFO))
  }


  if(!dir.exists(output)){
    dir.create(output, recursive = TRUE)
  }
  flog.info('Done.')
  if(rank != 'ASV'){
    data.glom <- tax_glom(data, taxrank=rank)
  } else {
    data.glom <- data
  }

  # save.image("debug.rdata")
  # quit()
  if(is.null(comp)){
    fun <- paste('combinaisons <- combn(na.omit(unique(sample_data(data.glom)$',column1,')),2) ',sep='')
    eval(parse(text=fun))
  }else{
    comp_list <- unlist(strsplit(comp,","))
    # combinaisons <- combn(comp_list,2)
    combinaisons <- matrix(, nrow = 2, ncol = length(comp_list))
    for (i in 1:length(comp_list)){
      tmp <- unlist(strsplit(comp_list[i],"\\~"))
      # cbind(combinaisons,tmp)
      combinaisons[1,i] <- tmp[1]
      combinaisons[2,i] <- tmp[2]
    }
  }

  pdf(file=paste(output,'/deseq2_',column1,'.pdf',sep=''),width=15,height=16)


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
    print(nb_cond1)
    fun <- paste('nb_cond2 <- nsamples(subset_samples(data.glom, ',column1, ' == "',combinaisons[2,col],'"))',sep='')
    eval(parse(text=fun))
    print(nb_cond2)

    if(nb_cond1 < 2){
      flog.warn(paste('Condition: ',combinaisons[1,col],' only have ',nb_cond1,' sample(s). Skip...',sep=''))
      tab0 = data.frame(matrix(ncol = 14, nrow = 0))
      write.table(tab0, file = paste(output,'/signtab_',column1,'_',paste(combinaisons[,col],collapse="_vs_"),'.csv',sep=''),quote=FALSE,sep="\t", row.names=FALSE)
      next
    } else if (nb_cond2 < 2){
      flog.warn(paste('Condition: ',combinaisons[2,col],' only have ',nb_cond2,' sample(s). Skip...',sep=''))
      tab0 = data.frame(matrix(ncol = 14, nrow = 0))
      write.table(tab0, file = paste(output,'/signtab_',column1,'_',paste(combinaisons[,col],collapse="_vs_"),'.csv',sep=''),quote=FALSE,sep="\t", row.names=FALSE)
      next
    } else if (nb_cond1 < 2 & nb_cond2 < 2){
      tab0 = data.frame(matrix(ncol = 14, nrow = 0))
      write.table(tab0, file = paste(output,'/signtab_',column1,'_',paste(combinaisons[,col],collapse="_vs_"),'.csv',sep=''),quote=FALSE,sep="\t", row.names=FALSE)
      flog.warn(paste('Condition: ',combinaisons[1,col],' only have ',nb_cond1,' sample(s). Skip...',sep=''))
      flog.warn(paste('Condition: ',combinaisons[2,col],' only have ',nb_cond2,' sample(s). Skip...',sep=''))
      next
    } else {
      flog.warn(paste('Condition: ',combinaisons[1,col],' contains ',nb_cond1,' samples.',sep=''))
      flog.warn(paste('Condition: ',combinaisons[2,col],' contains ',nb_cond2,' samples.',sep=''))
    }



    fun <- paste('tmp <- subset_samples(data.glom, ',column1, ' %in% c("',combinaisons[1,col],'","',combinaisons[2,col],'"))',sep='')
    eval(parse(text=fun))

    # if(nsamples(tmp)<4){
    # 	flog.warn('Less than 4 samples')
    # 	next;
    # }

    tmp <- prune_taxa(taxa_sums(tmp) >= 1, tmp)

    flog.info('DESeq2...')
    fun <- paste('deseq <- phyloseq_to_deseq2(tmp, ~ ',column1,')',sep='')
    eval(parse(text = fun))

    gm_mean = function(x, na.rm=TRUE){
      exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    }
    geoMeans = apply(DESeq2::counts(deseq), 1, gm_mean)
    deseq = estimateSizeFactors(deseq, geoMeans = geoMeans)

    flog.info('DESeq2...')
    deseq = DESeq2::DESeq(deseq, test="Wald", fitType="parametric")
    flog.debug(show(deseq))

    res = results(deseq, cooksCutoff = FALSE)
    flog.debug(show(res))
    alpha = 0.05

    if(length(which(res$padj < alpha)) > 0){
      # sigtab = res[which(res$padj < alpha), ]
      sigtab = res
      sigtab = cbind(row.names(sigtab),as(sigtab, "data.frame"), as(tax_table(data)[rownames(sigtab), ], "matrix"))
      colnames(sigtab)[1]=resultsNames(deseq)[2]
      save.image("debug.rdata")
      write.table(sigtab, file = paste(output,'/signtab_',column1,'_',paste(combinaisons[,col],collapse="_vs_"),'.csv',sep=''),quote=FALSE,sep="\t", row.names=FALSE)

      if(rank != 'ASV'){
        fun <- paste('x = tapply(sigtab$log2FoldChange, sigtab$',rank,', function(x) max(x))',sep='')
        fun2 <- paste('sigtab$',rank,' = factor(as.character(sigtab$',rank,'), levels=names(x))')
      } else {
        fun <- paste('x = tapply(sigtab$log2FoldChange, row.names(sigtab), function(x) max(x))',sep='')
        fun2 <- paste('sigtab$',rank,' = factor(as.character(row.names(sigtab)), levels=names(x))')
      }
      eval(parse(text=fun))
      x = sort(x, TRUE)
      eval(parse(text=fun2))

      fun <- paste('p <- ggplot(sigtab,
			aes(x=',rank,',
			y=log2FoldChange,
			color=',rank,')) + geom_point(size=6) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))',sep='')
      eval(parse(text=fun))

      ggtable <- ggtexttable(sigtab[,c("baseMean","log2FoldChange","stat","pvalue","padj",rank)],theme = ttheme("mOrange"),rows=NULL)
      grid.arrange(p,ggtable,top=text_grob(paste('Combination ',combinaisons[1,col], ' VS ' , combinaisons[2,col],sep=''), size=20))
    } else{
      flog.info(paste('No significant results for comparison ',combinaisons[1,col], ' ' , combinaisons[2,col], sep=''))
      tab0 = data.frame(matrix(ncol = 14, nrow = 0))
      # print(length(c(resultsNames(deseq)[2], colnames(res), colnames(tax_table(data)))))
      # print(c(resultsNames(deseq)[2], colnames(res), colnames(tax_table(data))))
      colnames(tab0) <- c(resultsNames(deseq)[2], colnames(res), colnames(tax_table(data)))
      write.table(tab0, file = paste(output,'/signtab_',column1,'_',paste(combinaisons[,col],collapse="_vs_"),'.csv',sep=''),quote=FALSE,sep="\t", row.names=FALSE)

    }
  }

  invisible(dev.off())




}
