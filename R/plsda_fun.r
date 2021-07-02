#' PLSDA from MixOmics package
#'
#' Multivariate methods for discriminant analysis.
#'
#' @param data A phyloseq object (output from decontam or generate_phyloseq)
#' @param output Output directory
#' @param column1 Factor to test.
#' @param rank Taxonomy rank to merge features that have same taxonomy at a certain taxonomic rank (among rank_names(data), or 'ASV' for no glom)
#'
#'
#' @return Return a list object with plots, and loadings table.
#'
#' @import futile.logger
#' @import dada2
#' @import phyloseq
#' @import DECIPHER
#' @import ShortRead
#' @import Biostrings
#' @import mixOmics
#' @export


plsda_fun <- function(data = data, output = "./plsda/", column1 = "",
            rank = "ASV"){

  outF <- list()

  if(!dir.exists(output)){
  	dir.create(output,recursive=T)
  }


  flog.info('Preparing tables...')
  if(rank != 'ASV'){
    flog.info('Glom rank...')
    physeqDA <- tax_glom(data, taxrank=rank)
    otable <- otu_table(physeqDA)
    row.names(otable) <- tax_table(physeqDA)[,rank]
  } else {
    physeqDA <- data
    otable <- otu_table(physeqDA)
  }
  mdata <- sample_data(physeqDA)
  flog.info('Done.')

  flog.info('Removing levels with only one sample...')
  fun <- paste('lvl_to_remove <- names(table(mdata$',column1,')[table(mdata$',column1,') <= 1])',sep='')
  eval(parse(text=fun))
  flog.debug(paste0('Level "',lvl_to_remove, '" will be removed.'))
  if(length(lvl_to_remove) > 0){
  	fun <- paste('sample_to_remove <- rownames(mdata[mdata$',column1,' == lvl_to_remove])')
    eval(parse(text=fun))
  	flog.info(paste('Removing ',sample_to_remove,sep=''))
  	otable <- otable[,colnames(otable)!=sample_to_remove]
  	mdata <- mdata[rownames(mdata)!=sample_to_remove]

  }
  flog.info('Done.')

  flog.info('PLSDA...')
  fun <- paste('plsda.res = plsda(t(otable+1), mdata$',column1,', ncomp = 5, logratio="CLR")',sep='')
  eval(parse(text=fun))
  flog.info('Done.')

  flog.info('Plotting PLSDA individuals...')
  png(paste(output,'/plsda_indiv_',column1,'_',rank,'.png',sep=''))
  background = background.predict(plsda.res, comp.predicted=2, dist = "max.dist")
  fun <- paste('outF$plsda.plotIndiv <- plotIndiv(plsda.res,
  	comp= 1:2,
  	group = mdata$',column1,',
  	ind.names=FALSE,
  	ellipse=TRUE,
  	legend=TRUE,
  	title= "PLSDA plot of individuals",
  	background = background)',sep='')
  eval(parse(text=fun))
  dev.off()
  flog.info('Done.')

  plot_plsda_perf <- function(){
  	flog.info('Plotting PLSDA performance...')
  	perf.plsda <- perf(plsda.res, validation = "Mfold", folds = 5, progressBar = FALSE, auc = TRUE, nrepeat = 10)

  	plotperf <- plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, legend.position = "horizontal", title = "PLSDA performance plot")
  	# plotperf <- recordPlot()
  	# invisible(dev.off())
  	#
  	# png(paste(output,'/plsda_perf_',column1,'_',rank,'.png',sep=''))
  	# replayPlot(plotperf)
  	# dev.off()
  	flog.info('Done.')

    return(plotperf)
  }

  # tryCatch(plot_perf(), error = function(e) { flog.warn("PLSDA perf function not working.")})
  # outF$splsda.plotPerf <- plot_plsda_perf()



  tune_splsda <- function(){
  	flog.info('Tune SPLDA...')
  	fun <- paste('tune.splsda <- tune.splsda(t(otable+1),
  	mdata$',column1,
  	', ncomp = 4,
  	validation = "Mfold",
  	folds = 4,
  	progressBar = FALSE,
  	dist = "max.dist",
  	nrepeat = 10)',sep='')
  	eval(parse(text=fun))

  # 	plot(tune.splsda, col = color.jet(4), title = "Error rates SPLSDA")
  # 	outF$splsda.plotError <- recordPlot()
  # 	invisible(dev.off())
  #
  # 	png(paste(output,'/splsda_error_',column1,'_',rank,'.png',sep=''))
  #   replayPlot(outF$splsda.plotError)
  # 	dev.off()

  	ncomp <- tune.splsda$choice.ncomp$ncomp + 1
  	select.keepX <- tune.splsda$choice.keepX[1:ncomp-1]
  	r_lst <- list("ncomp" = ncomp, "selectkeepX" = select.keepX)
  	flog.info('Done.')
  	return(r_lst)
  }

  # tryCatch(tune_splsda(), error = function(e) { stop("At least one class is not represented in one fold.") })
  r_list <- tune_splsda()
  ncomp <- r_list$ncomp
  select.keepX <- r_list$selectkeepX

  flog.info(paste('keepX: ',select.keepX,sep=''))

  flog.info('SPLSDA...')
  fun <- paste('splsda.res <- splsda(t(otable+1), mdata$',column1,', ncomp = ncomp, keepX = select.keepX, logratio= "CLR")',sep='')
  eval(parse(text=fun))
  flog.info('Done.')
  flog.info('Plot Individuals...')
  png(filename=paste(output,'/splsda_indiv_',column1,'_',rank,'.png',sep=''),width=480,height=480)
  fun <- paste('outF$splsda.plotIndiv <- plotIndiv(splsda.res, comp= c(1:2), group = mdata$',column1,', ind.names = FALSE, ellipse = TRUE, legend = TRUE, title = "sPLS-DA on ',column1,'")',sep='')
  eval(parse(text=fun))
  dev.off()
  flog.info('Done.')

  flog.info('SPLSDA performance...')
  perf.splsda <- perf(splsda.res, validation = "Mfold", folds = 5,
                   dist = 'max.dist', nrepeat = 10,
                   progressBar = FALSE)

  png(paste(output,'/splsda_perf_',column1,'_',rank,'.png',sep=''))
  plot(perf.splsda, col = color.mixo(5))
  dev.off()


  outF[["loadings"]] = list()
  for (comp in 1:ncomp){
  	plotLoadings(splsda.res, comp = comp, title = paste('Loadings on comp',comp,sep=''), contrib = 'max', method = 'mean')
    outF$loadings[[glue::glue("comp{comp}")]] <- recordPlot()
    invisible(dev.off())

    png(paste(output,'/splsda_loadings_',column1,'_',rank,'_comp',comp,'.png',sep=''),width=800,height=800)
    replayPlot(outF$loadings[[glue::glue("comp{comp}")]])
  	dev.off()
  }

  # plotArrow(splsda.res, legend=T)
  # outF$splsda.plotArrow <- recordPlot()
  # invisible(dev.off())
  #
  # png(paste(output,'/splsda_arrow_',column1,'_',rank,'.png',sep=''))
  # replayPlot(outF$splsda.plotArrow)
  # dev.off()

  outF$splsda.loadings_table = splsda.res$loadings$X
  write.table(splsda.res$loadings$X,paste(output,'/splsda_table_',column1,'_',rank,'.csv',sep=''),quote=FALSE,sep="\t",col.names=NA)
  flog.info('Done.')

  return(outF)

}
