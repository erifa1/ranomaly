#' DESEQ2 Differential Analysis function.
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
#' @return Returns list with table of features and plots for each comparison. Exports CSV files with significant differentialy abondant ASV.
#'
#' @import futile.logger
#' @import phyloseq
#' @importFrom BiocGenerics estimateSizeFactors
#' @importFrom DESeq2 DESeq
#' @importFrom DESeq2 results
#' @importFrom DESeq2 resultsNames
#' @importFrom gridExtra grid.arrange
#' @importFrom gridExtra tableGrob
#' @importFrom grid textGrob
#' @importFrom plotly ggplotly
#' @import ggplot2

#'
#' @export


# Decontam Function

deseq2_fun <- function(data = data, output = "./deseq/", column1 = "", verbose = 1,
                       rank = "Species", comp = ""){


  if(verbose == 3){
    invisible(flog.threshold(DEBUG))
  } else {
    invisible(flog.threshold(INFO))
  }


  if(!dir.exists(output)){
    dir.create(output, recursive = TRUE)
  }

  if(rank != 'ASV'){
    flog.info('Glom rank...')
    data.glom <- tax_glom(data, taxrank=rank)
  } else {
    data.glom <- data
  }

  # save.image("debug.rdata")
  # quit()
  flog.info('Defining comparison...')
  if(comp==""){
    fun <- paste('combinaisons <- utils::combn( as.character( na.omit(unique(sample_data(data.glom)$',column1,')) ),2)',sep='')
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

  print(combinaisons)
  flog.info('Done...')

  outF = list()
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
    flog.debug( print(nb_cond1) )
    fun <- paste('nb_cond2 <- nsamples(subset_samples(data.glom, ',column1, ' == "',combinaisons[2,col],'"))',sep='')
    eval(parse(text=fun))
    flog.debug( print(nb_cond2) )

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
    # return(tmp)

    flog.info('DESeq2...')
    fun <- paste('dseq <- phyloseq_to_deseq2(tmp, ~ ',column1,')',sep='')
    eval(parse(text = fun))

    gm_mean = function(x, na.rm=TRUE){
      exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    }
    geoMeans = apply(DESeq2::counts(dseq), 1, gm_mean)
    dseq2 = estimateSizeFactors(dseq, geoMeans = geoMeans)

    # flog.info('DESeq2...')
    # print(dseq2)
    dseq3 = try(DESeq2::DESeq(dseq2, test="Wald", fitType="parametric"))
    if(class(dseq3) == "try-error"){next}
    flog.debug(show(dseq3))

    res = results(dseq3, cooksCutoff = FALSE, contrast = c(column1,combinaisons[1,col] , combinaisons[2,col]))
    flog.debug(show(res))

    # outlist1 <- list()

    flog.info('Plot...')
    alpha = 0.05
    if(length(which(res$padj < alpha)) > 0){
      sigtab = as.data.frame(res[which(res$padj < alpha), ])
      sigtab = cbind(taxon_id = row.names(sigtab),as(sigtab, "data.frame"), as(tax_table(data)[rownames(sigtab), ], "matrix"))
      colnames(sigtab)[1]=resultsNames(dseq3)[2]
      
      # Output Global Table
      tabOUT = cbind(taxon_id = row.names(res),as(res, "data.frame"), as(tax_table(data)[rownames(res), ], "matrix"))
      colnames(res)[1]=resultsNames(dseq3)[2]
      write.table(tabOUT, file = paste(output,'/signtab_',column1,'_',paste(combinaisons[,col],collapse="_vs_"),'.csv',sep=''),quote=FALSE,sep="\t", row.names=FALSE)

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
			color=',rank,')) + geom_point(size=6) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5), legend.position = "none")',sep='')
      eval(parse(text=fun))

      Ftable = sigtab[,c("baseMean","log2FoldChange","stat","pvalue","padj",rank)]
      ggtable <- gridExtra::tableGrob(Ftable)

      pdf(file=paste(output,'/deseq2_',column1,"_", paste(combinaisons[,col],collapse="_vs_"), '.pdf',sep=''),width=15,height=16)
      grid.arrange(p,ggtable,top=grid::textGrob(paste('Combination ',combinaisons[1,col], ' VS ' , combinaisons[2,col],sep=''), size=20))
      invisible(dev.off())


      outF[[paste(combinaisons[,col],collapse="_vs_")]] = list(plot = p, table = Ftable)


    } else{
      flog.info(paste('No significant results for comparison ',combinaisons[1,col], ' ' , combinaisons[2,col], sep=''))
      headers <- c(resultsNames(dseq3)[2], colnames(res), colnames(tax_table(data)))
      tab0 = data.frame(matrix(ncol = length(headers), nrow = 0))
      # print(length(c(resultsNames(deseq)[2], colnames(res), colnames(tax_table(data)))))
      # print(c(resultsNames(deseq)[2], colnames(res), colnames(tax_table(data))))
      colnames(tab0) <- headers
      write.table(tab0, file = paste(output,'/signtab_',column1,'_',paste(combinaisons[,col],collapse="_vs_"),'.csv',sep=''),quote=FALSE,sep="\t", row.names=FALSE)

    }
  }
  return(outF)
  flog.info('Done...')
}
