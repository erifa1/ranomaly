#' Rarefaction plotly
#'
#' @param data a phyloseq object (output from decontam or generate_phyloseq)
#' @param col The metadata column name you want to split the graph.
#' @param step The number of steps you want your rarefaction plots to be done.
#'
#' @return A plotly graph.
#'
#' @import plotly
#' @importFrom ranacapa ggrare
#'
#' @export


rarefaction <- function(data = data, col = NULL, step = 100){
  plot_rare <- ggrare(data, step = step, color = col, plot = FALSE)
  plot_rare <- plot_rare + facet_wrap(col, ncol = 4) + theme_bw()
  return(ggplotly(plot_rare))
}


#' Barplots plotly
#'
#'
#' @param data a phyloseq object (output from decontam or generate_phyloseq)
#' @param rank Taxonomy rank to merge features that have same taxonomy at a certain taxonomic rank (among rank_names(data), or 'ASV' for no glom)
#' @param top Number of top taxa to plot
#' @param Ord1 Variable used to order sample (X axis)
#' @param Fact1 Variable used to change X axis tick labels and color
#' @param relative Plot relative (TRUE, default) or raw abundance plot (FALSE)
#' @param output Output directory
#'
#' @return Exports barplots in an interactive plotly community plot
#'
#'
#' @import plotly
#' @importFrom microbiome aggregate_top_taxa
#' @importFrom reshape2 melt
#'
#'
#'
#' @export


bars_fun <- function(data = data, rank = "Genus", top = 10, Ord1 = NULL, Fact1 = NULL, relative = TRUE){

    # Fdata <- prune_samples(sample_names(r$data16S())[r$rowselect()], r$data16S())
    # Fdata <- prune_taxa(taxa_sums(Fdata) > 0, Fdata)
    # if( r$RankGlom() == "ASV"){
    #   Fdata <- prune_taxa(r$asvselect(), Fdata)
    # }

    Fdata = data
    # print("top")
    psobj.top <- aggregate_top_taxa(Fdata, rank, top = top)

    # print("get data")
    sdata = as.data.frame(sample_data(psobj.top))
    sdata$sample.id = sample_names(psobj.top)
    otable = as.data.frame(otu_table(psobj.top))
    row.names(otable) = tax_table(psobj.top)[,rank]

    # print("melt data")
    dat <- as.data.frame(t(otable))
    dat <- cbind.data.frame(sdata, dat)
    meltdat <- reshape2::melt(dat, id.vars=1:ncol(sdata))
    tt <- levels(meltdat$variable)
    meltdat$variable <- factor(meltdat$variable, levels= c("Other", tt[tt!="Other"]))

    LL=list()
    # print(head(meltdat))
    # print(levels(meltdat$sample.id))

    fun = glue( "xform <- list(categoryorder = 'array',
                    categoryarray = unique(meltdat$sample.id[order(meltdat${Ord1})]),
                    title = 'Samples',
                    tickmode = 'array',
                    tickvals = 0:nrow(sdata),
                    ticktext = sdata[unique(meltdat$sample.id[order(meltdat${Ord1})]), '{Fact1}']@.Data[[1]],
                    tickangle = -45)")
    eval(parse(text=fun))

    # subplot to vizualize groups
    # print(head(sdata))
    df1 <- cbind.data.frame(x=sdata[unique(meltdat$sample.id[order(meltdat[,Ord1])]), "sample.id"]@.Data[[1]],
                            g=sdata[unique(meltdat$sample.id[order(meltdat[,Ord1])]), Fact1]@.Data[[1]],
                            y=1)

    subp1 <- df1 %>% plot_ly(
      type = 'bar',
      x = ~x,
      y = ~y,
      color = ~g,
      legendgroup = ~g,
      showlegend = FALSE
    ) %>% layout(xaxis = list(zeroline = FALSE,showline = FALSE, showgrid = FALSE),
                 yaxis=list(showticklabels = FALSE,title = "",showgrid = FALSE))


if(relative){
  #relative abondance
  otable=apply(otable,2, function(x){Tot=sum(x); x/Tot})
  dat= as.data.frame(t(otable))
  dat <- cbind.data.frame(sdata, dat)
  meltdat = data.table::melt(dat, id.vars=1:ncol(sdata))
  tt=levels(meltdat$variable)
  meltdat$variable = factor(meltdat$variable, levels= c("Other", tt[tt!="Other"]))

  p1=plot_ly(meltdat, x = ~sample.id, y = ~value, type = 'bar', name = ~variable, color = ~variable) %>% #, color = ~variable
    layout(title="Relative abundance", yaxis = list(title = 'Relative abundance'), xaxis = xform, barmode = 'stack')

  p1 <- subplot(p1, subp1, nrows = 2, shareX = T, heights=c(0.95,0.05)) %>%
    layout(xaxis = xform)
}else{
  #raw abundance
  p1=plot_ly(meltdat, x = ~sample.id, y = ~value, type = 'bar', name = ~variable, color = ~variable) %>% #, color = ~variable
    layout(title="Raw abundance", yaxis = list(title = 'Raw abundance'), xaxis = xform, barmode = 'stack')

  p1 <- subplot(p1, subp1, nrows = 2, shareX = T, heights=c(0.95,0.05)) %>%
    layout(xaxis = xform)
}

return(p1)

}
