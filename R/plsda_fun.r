#' PLSDA from MixOmics package
#'
#' Multivariate methods for discriminant analysis.
#'
#' @param data A phyloseq object (output from decontam or generate_phyloseq)
#' @param output Output directory
#' @param column1 Factor to test.
#' @param rank Taxonomy rank to merge features that have same taxonomy at a certain taxonomic rank (among rank_names(data), or 'ASV' for no glom)
#' @param axis Select the axis to plot
#' @param multilevel Factor (metadata column name) used within matrix decomposition for repeated measurements
#' @param ind.names either a character vector of names for the individuals to be plotted, or FALSE for no names. If TRUE, the row names of the first (or second) data matrix is used as names (see mixOmics Details).
#' @param ellipse Logical indicating if ellipse plots should be plotted (see mixOmics Details).
#' @param progressBar Silence the progress bar.
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


plsda_fun <- function (data = data, output = "./plsda/", column1 = "",
                       rank = "ASV", axis = c(1, 2), multilevel = NULL,
                       ind.names = FALSE,
                       ellipse = FALSE, progressBar = FALSE)
{
  outF <- list()
  if (!dir.exists(output)) {
    dir.create(output, recursive = T)
  }
  flog.info("Preparing tables...")
  if (rank != "ASV") {
    flog.info("Glom rank...")
    physeqDA <- tax_glom(data, taxrank = rank)
    otable <- otu_table(physeqDA)
    row.names(otable) <- tax_table(physeqDA)[, rank]
  } else {
    physeqDA <- data
    otable <- otu_table(physeqDA)
  }
  mdata <- sample_data(physeqDA)
  flog.info("Done.")
  flog.info("Removing levels with only one sample...")
  fun <- paste("lvl_to_remove <- names(table(mdata$", column1,
               ")[table(mdata$", column1, ") <= 1])", sep = "")
  eval(parse(text = fun))
  flog.debug(paste0("Level \"", lvl_to_remove, "\" will be removed."))
  if (length(lvl_to_remove) > 0) {
    fun <- paste("sample_to_remove <- rownames(mdata[mdata$",
                 column1, " == lvl_to_remove])")
    eval(parse(text = fun))
    flog.info(paste("Removing ", sample_to_remove, sep = ""))
    otable <- otable[, colnames(otable) != sample_to_remove]
    mdata <- mdata[rownames(mdata) != sample_to_remove]
  }
  flog.info("Done.")
  flog.info("PLSDA...")
  if (taxa_are_rows(physeqDA)) {
    otable <- t(otable + 1)
  } else {
    otable <- otable + 1
  }
  if (!is.null(multilevel)) {
    multilevel <- mdata[, multilevel] |> data.frame()
    Y <- mdata[, column1]
    multilevel <- data.frame(multilevel, Y)
    multilevel[, 1] <- multilevel[, 1] |> factor() |> as.numeric()
    otable <- withinVariation(otable, multilevel)
    fun <- paste("plsda.res <- plsda(otable, as.numeric(factor(multilevel[, 2])),
                       ncomp = 5, logratio=\"none\")", sep = "")
  } else {
    fun <- paste("plsda.res <- plsda(otable, mdata$", column1,
                 ", ncomp = 5, logratio=\"CLR\")", sep = "")
  }
  eval(parse(text = fun))
  flog.info("Done.")
  flog.info("Plotting PLSDA individuals...")
  png(paste(output, "/plsda_indiv_", column1, "_", rank, ".png",
            sep = ""))
  background = background.predict(plsda.res, comp.predicted = 2,
                                  dist = "max.dist")
  if (all(axis == c(1, 2))) {
    fun <- paste("outF$plsda.plotIndiv <- plotIndiv(plsda.res,\n  \tcomp= 1:2,\n  \tgroup = mdata$",
                 column1,
                 ",\n   \tind.names=",
                 ind.names,
                 ",\n  \tellipse=",
                 ellipse,
                 ",\n  \tlegend=TRUE,\n  \ttitle= \"PLSDA plot of individuals\",\n  \tbackground = background)",
                 sep = "")
  } else if (all(axis == c(2, 3))) {
    fun <- paste("outF$plsda.plotIndiv <- plotIndiv(plsda.res,\n  \tcomp= 2:3,\n  \tgroup = mdata$",
                 column1,
                 ",\n   \tind.names=",
                 ind.names,
                 ",\n  \tellipse=",
                 ellipse,
                 ",\n  \tlegend=TRUE,\n  \ttitle= \"PLSDA plot of individuals\",\n  \tbackground = background)",
                 sep = "")
  }
  eval(parse(text = fun))
  dev.off()
  flog.info("Done.")
  plot_plsda_perf <- function() {
    flog.info("Plotting PLSDA performance...")
    perf.plsda <- perf(plsda.res, validation = "Mfold", folds = 5,
                       progressBar = progressBar, auc = TRUE, nrepeat = 10)
    plotperf <- plot(perf.plsda, col = color.mixo(1:3), sd = TRUE,
                     legend.position = "horizontal", title = "PLSDA performance plot")
    flog.info("Done.")
    return(plotperf)
  }
  tune_splsda <- function() {
    flog.info("Tune SPLSDA...")
    if (!is.null(multilevel)) {
      fun <- paste("tune.splsda <- tune.splsda(otable, as.numeric(factor(multilevel[, 2])),
                   ncomp = 4,\n  \tvalidation = \"Mfold\",\n  \tfolds = 4,\n  \tprogressBar = ",
                   progressBar,
                   ",\n  \tdist = \"max.dist\",\n  \tnrepeat = 10)",
                   sep = "")
    } else {
      fun <- paste("tune.splsda <- tune.splsda(otable,\n  \tmdata$",
                   column1,
                   ", ncomp = 4,\n  \tvalidation = \"Mfold\",\n  \tfolds = 4,\n  \tprogressBar = ",
                   progressBar,
                   ",\n  \tdist = \"max.dist\",\n  \tnrepeat = 10)",
                   sep = "")
    }
    eval(parse(text = fun))
    ncomp <- tune.splsda$choice.ncomp$ncomp + 1
    select.keepX <- tune.splsda$choice.keepX[1:ncomp - 1]
    r_lst <- list(ncomp = ncomp, selectkeepX = select.keepX)
    flog.info("Done.")
    return(r_lst)
  }
  r_list <- tune_splsda()
  ncomp <- r_list$ncomp
  select.keepX <- r_list$selectkeepX
  flog.info(paste("keepX: ", select.keepX, sep = ""))
  flog.info("SPLSDA...")
  if (!is.null(multilevel)) {
    fun <- paste("splsda.res <- splsda(otable, as.numeric(factor(multilevel[, 2])),
                 ncomp = ncomp, keepX = select.keepX, logratio= \"none\")",
                 sep = "")
  } else {
    fun <- paste("splsda.res <- splsda(otable, mdata$",
                 column1, ", ncomp = ncomp, keepX = select.keepX, logratio= \"CLR\")",
                 sep = "")
  }
  eval(parse(text = fun))
  flog.info("Done.")
  flog.info("Plot Individuals...")
  png(filename = paste(output, "/splsda_indiv_", column1, "_",
                       rank, ".png", sep = ""), width = 480, height = 480)
  if (all(axis == c(1, 2))) {
    fun <- paste("outF$splsda.plotIndiv <- plotIndiv(splsda.res, comp= c(1:2), group = mdata$",
                 column1,
                 ",\n   \tind.names=",
                 ind.names,
                 ",\n  \tellipse=",
                 ellipse,
                 ",\n   \tlegend = TRUE, title = \"sPLS-DA on ",
                 column1, "\")", sep = "")
  } else if (all(axis == c(2, 3))) {
    fun <- paste("outF$splsda.plotIndiv <- plotIndiv(splsda.res, comp= c(2:3), group = mdata$",
                 column1,
                 ",\n   \tind.names=",
                 ind.names,
                 ",\n  \tellipse=",
                 ellipse,
                 ",\n   \tlegend = TRUE, title = \"sPLS-DA on ",
                 column1, "\")", sep = "")
  }
  eval(parse(text = fun))
  dev.off()
  flog.info("Done.")
  flog.info("SPLSDA performance...")
  perf.splsda <- perf(splsda.res, validation = "Mfold", folds = 5,
                      dist = "max.dist", nrepeat = 10, progressBar = progressBar)
  png(paste(output, "/splsda_perf_", column1, "_", rank, ".png",
            sep = ""))
  plot(perf.splsda, col = color.mixo(5))
  dev.off()
  outF[["loadings"]] = list()
  for (comp in 1:ncomp) {
    plotLoadings(splsda.res, comp = comp, title = paste("Loadings on comp",
                                                        comp, sep = ""), contrib = "max", method = "mean")
    outF$loadings[[glue::glue("comp{comp}")]] <- recordPlot()
    invisible(dev.off())
    png(paste(output, "/splsda_loadings_", column1, "_",
              rank, "_comp", comp, ".png", sep = ""), width = 800,
        height = 800)
    replayPlot(outF$loadings[[glue::glue("comp{comp}")]])
    dev.off()
  }

  outF[["var"]] = list()
  for (comp in 2:3) {
    plotVar(splsda.res, comp = c(1,comp), title = paste("Loadings on comp 1 and ",
                                                        comp, sep = ""))
    outF$var[[glue::glue("comp{comp}")]] <- recordPlot()
    invisible(dev.off())
    png(paste(output, "/splsda_var_", column1, "_",
              rank, "_comp1-", comp, ".png", sep = ""), width = 800,
        height = 800)
    replayPlot(outF$var[[glue::glue("comp{comp}")]])
    dev.off()
  }

  outF$splsda.loadings_table = splsda.res$loadings$X
  write.table(splsda.res$loadings$X, paste(output, "/splsda_table_",
                                           column1, "_", rank, ".csv", sep = ""), quote = FALSE,
              sep = "\t", col.names = NA)
  flog.info("Done.")
  outF$splsda.res = splsda.res
  outF$plsda.res = plsda.res
  return(outF)
}
