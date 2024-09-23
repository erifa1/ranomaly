#' @title Remove non-ASCII characters from a string
#' @keywords internal
#' @param text A string
remove_non_ascii <- function(text) {
  cleaned_text <- iconv(text, from = "UTF-8", to = "ASCII//TRANSLIT", sub = "")
  return(cleaned_text)
}

# Raised by devtools::check()
globalVariables(c("<<-", ".", "Abundance", "Condition", "DESeqLFC", "Sample", "Tax", "aggregate", "aov", "as.formula", "as.phylo", "cLS", "contaminant", "data", "data.pa.neg", "data.pa.pos", "dev.off", "dist", "dseq", "freqASV", "ggs", "graphics.off", "hclust", "i", "lme1", "lvl_to_remove", "mclapply", "median", "model.matrix", "na.omit", "nb_cond1", "nb_cond2", "opt_parser", "p", "pa.neg", "pa.pos", "padj", "pander", "pdf", "plsda.res", "png", "print_help", "psobj", "rainbow", "read.table", "recordPlot", "reorder", "replayPlot", "sample_to_remove", "setNames", "setTxtProgressBar", "splsda.res", "sumASV", "tax", "trainingSet", "txtProgressBar", "update", "where", "wilcox.test", "wilcox_res", "wilcox_res1", "write.table"))

# Dataset Documentation
#' Lookup-table for IDs of taxonomic ranks (metacoder)
#'
#' Composed of two columns:
#' \itemize{
#'  \item rankid - the ordered identifier value. lower values mean higher rank
#'  \item ranks - all the rank names that belong to the same level, with
#'  different variants that mean essentially the same thing
#' }
#'
#' @name ranks_ref
#' @docType data
#' @keywords data
NULL


.onLoad <- function(libname, pkgname) {
  packageStartupMessage(paste("Package", pkgname, "version", utils::packageVersion(pkgname), "loaded."))
  packageStartupMessage(paste0("To cite this package in publications, please use ", "citation('",pkgname, "').") )
}

