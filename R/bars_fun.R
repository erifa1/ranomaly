#' Barplots community
#'
#'
#' @param dada_res output from dada2_fun
#'
#' @return Return raw otu table in phyloseq object.
#' @import phyloseq
#' @import ggplot2
#' @import gridExtra
#' @import grid
#' @import microbiome
#' @import viridis
#' @import plotly
#'
#'
#'
#' @export


# Decontam Function

bars_fun <- function(data = data, bar = TRUE, compo = TRUE, output = "./plot_bar/", column1 = "", column2 = "",
                     sname = FALSE, num = 10, rare = NULL, rank = "Genus"){

  suppressMessages(source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/load-extra-functions.R"))


}
