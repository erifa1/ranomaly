#' Rarefaction plotly
#'
#' @param data a phyloseq object (output from decontam or generate_phyloseq)
#' @param col The metadata column name you want to split the graph.
#' @param step The number of steps you want your rarefaction plots to be done.
#' @param ggplotly Output ggplotly version if TRUE.
#'
#' @return A plotly graph.
#'
#' @import plotly
#' @importFrom ranacapa ggrare
#'
#' @export


rarefaction <- function(data = data, col = NULL, step = 100, ggplotly = TRUE){
  plot_rare <- ggrare(data, step = step, color = col, plot = FALSE)
  plot_rare <- plot_rare + facet_wrap(col, ncol = 4) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust=1),
          legend.position = "none",axis.text=element_text(size=18, angle = 45, hjust=1),
          axis.title=element_text(size=16,face="bold"),
          strip.text.x = element_text(size = 18,face="bold"),
          title=element_text(size=16,face="bold"))


  if(ggplotly){
    return(ggplotly(plot_rare))
  }else{
    return(plot_rare)
  }
}


#' aggregate_top_taxa from microbiome package
#'
#'
#' @param x phyloseq object
#' @param top Keep the top-n taxa, and merge the rest under the category 'Other'. Instead of top-n numeric this can also be a character vector listing the groups to combine.
#' @param level Summarization level (from ‘rank_names(pseq)’)
#'
#' @importFrom microbiome aggregate_taxa
#' @importFrom microbiome top_taxa
#'
#' @export

aggregate_top_taxa <- function (x, top, level){
    x <- aggregate_taxa(x, level)
    tops <- top_taxa(x, top)
    tax <- tax_table(x)
    inds <- which(!rownames(tax) %in% tops)
    tax[inds, level] <- "Other"
    tax_table(x) <- tax
    tt <- tax_table(x)[, level]
    tax_table(x) <- tax_table(tt)
    aggregate_taxa(x, level)
}



#' Barplots plotly
#'
#'
#' @param data a phyloseq object (output from decontam or generate_phyloseq)
#' @param rank Taxonomy rank to merge features that have same taxonomy at a certain taxonomic rank (among rank_names(data), or 'ASV' for no glom)
#' @param top Number of top taxa to plot
#' @param Ord1 Variable used to order sample (X axis) or split the barplot if split = TRUE
#' @param sample_labels If true, x axis labels are sample IDS, if false labels displayed are levels from Ord1 argument. Ignored if split = TRUE (FALSE)
#' @param split if TRUE make a facet_wrap like grouped by Ord1 (default FALSE)
#' @param relative Plot relative (TRUE, default) or raw abundance plot (FALSE)
#' @param autoorder Automatic ordering xaxis labels based on Ord1 factor levels with gtools::mixedorder function (TRUE).
#' @param ylab Y axis title ("Abundance")
#' @param outfile Output html file.
#'
#' @return Returns barplots in an interactive plotly community plot
#'
#' @import plotly
#' @importFrom reshape2 melt
#' @importFrom gtools mixedsort
#' @importFrom dplyr group_map group_by across
#'
#'
#' @export


bars_fun <- function(data = data, rank = "Genus", top = 10, Ord1 = NULL, sample_labels = FALSE, split = FALSE,
                     relative = TRUE, autoorder = TRUE, ylab = "Abundance", outfile="plot_compo.html", verbose = TRUE){

  if(verbose){
    invisible(flog.threshold(INFO))
  } else {
    invisible(flog.threshold(ERROR))
  }


if( all(Ord1 != sample_variables(data))){
  stop(paste("Wrong value in Ord1, please use variables existing in the phyloseq object:", toString(sample_variables(data))))
}

  flog.info('Preprocess...')
  Fdata = data
  psobj.top <- aggregate_top_taxa(Fdata, rank, top = top)

  sdata <- as.data.frame(sample_data(psobj.top), stringsAsFactors = TRUE)
  sdata$sample.id = sample_names(psobj.top)
  otable = as.data.frame(otu_table(psobj.top))
  row.names(otable) = tax_table(psobj.top)[,rank]

  # print("melt data")
  dat <- as.data.frame(t(otable))
  dat <- cbind.data.frame(sdata, dat)

  flog.info('  Melting table...')
  meltdat <- reshape2::melt(dat, id.vars=1:ncol(sdata))
  tt <- levels(meltdat$variable)
  meltdat$variable <- factor(meltdat$variable, levels= c("Other", tt[tt!="Other"]))

  LL=list()
  # print(head(meltdat))
  # print(levels(meltdat$sample.id))
  # save(list = ls(all.names = TRUE), file = "debug.rdata", envir = environment())


  if(autoorder){
  flog.info('  Ordering samples...')
      fun = glue( "labs = gtools::mixedorder(as.character(meltdat${Ord1}))" )
      eval(parse(text=fun))

      orderedIDS <- unique(meltdat$sample.id[gtools::mixedorder(as.character(meltdat[,Ord1]))])
      orderedOrd1 <- meltdat[,Ord1][gtools::mixedorder(as.character(meltdat[,Ord1]))]
      orderedOrd1 <- factor(orderedOrd1, levels = gtools::mixedsort(levels(orderedOrd1)))
    }else{
      labs = 1:nrow(meltdat)

      orderedIDS <- unique(meltdat$sample.id)
      orderedOrd1 <- meltdat[,Ord1]
    }

  if(sample_labels){
      lab1 = "sample.id"
    }else{
      lab1 = Ord1
    }

  flog.info('  Set labels...')
  xform <- list(categoryorder = 'array',
                categoryarray = unique(meltdat$sample.id[labs]),
                title = 'Samples',
                tickmode = 'array',
                tickvals = 0:nrow(sdata),
                ticktext = sdata[as.character(unique(meltdat$sample.id[labs])), lab1]@.Data[[1]],
                tickangle = -90)

  # subplot to vizualize groups

  flog.info('  Subplot...')
  df1 <- cbind.data.frame(x=sdata[orderedIDS, "sample.id"]@.Data[[1]],
                          g=sdata[orderedIDS, Ord1]@.Data[[1]],
                          y=1)

  fun = glue( "df1$g <- factor(df1$g, levels = as.character(unique(orderedOrd1)))")
  eval(parse(text=fun))
  fun = glue( "meltdat${Ord1} <- factor(meltdat${Ord1}, levels = as.character(levels(orderedOrd1)))")
  eval(parse(text=fun))

  subp1 <- df1 %>% plot_ly(
    type = 'bar',
    x = ~x,
    y = ~y,
    color = ~g,
    legendgroup = ~g,
    showlegend = FALSE
  ) %>% layout(xaxis = list(zeroline = FALSE,showline = FALSE, showgrid = FALSE),
               yaxis=list(showticklabels = FALSE,title = ylab, showgrid = FALSE))


  if(relative){
  flog.info('Plotting relative...')
    #relative abondance
    otable=apply(otable,2, function(x){Tot=sum(x); x/Tot})
    dat= as.data.frame(t(otable))
    dat <- cbind.data.frame(sdata, dat)

    meltdat <- reshape2::melt(dat, id.vars=1:ncol(sdata))
    tt <- levels(meltdat$variable)
    meltdat$variable <- factor(meltdat$variable, levels= c("Other", tt[tt!="Other"]))

    fun = glue( "meltdat${Ord1} <- factor(meltdat${Ord1}, levels = as.character(levels(orderedOrd1)))")
    eval(parse(text=fun))

    p1=plot_ly(meltdat, x = ~sample.id, y = ~value, type = 'bar', name = ~variable, color = ~variable) %>% #, color = ~variable
      layout(title="", yaxis = list(title = ylab), xaxis = xform, barmode = 'stack')

    if(length(df1$x) != length(unique(df1$g))){
      p1 <- subplot(p1, subp1, nrows = 2, shareX = T, heights=c(0.95,0.05)) %>%
      layout(xaxis = xform)
    }
  }else{
  flog.info('Plotting raw...')
    #raw abundance
    p1=plot_ly(meltdat, x = ~sample.id, y = ~value, type = 'bar', name = ~variable, color = ~variable) %>% #, color = ~variable
      layout(title="", yaxis = list(title = ylab), xaxis = xform, barmode = 'stack')

    if(length(df1$x) != length(unique(df1$g))){
      p1 <- subplot(p1, subp1, nrows = 2, shareX = T, heights=c(0.95,0.05)) %>%
      layout(xaxis = xform)
    }
  }

  # Splitted plot output
  if(!split) {
    if(!is.null(outfile)){
      htmlwidgets::saveWidget(p1, outfile)
    }
  flog.info('Finish...')
    return(p1)
  } else {
  flog.info('Splitted plot...')
    p1 = meltdat %>% group_by(across({Ord1})) %>%
      dplyr::group_map(~ plot_ly(data=., x = ~sample.id, y = ~value, type = 'bar',
                                 name = ~variable,
                                 color = ~variable, legendgroup = ~variable,
                                 showlegend = (.y == levels(meltdat[, Ord1])[1])),
                       keep = TRUE)  %>%
      subplot(nrows = 1, shareX = TRUE, shareY=TRUE, titleX = FALSE) %>%
      layout(title="",
             xaxis = list(title = glue("{Ord1} = {levels(meltdat[, Ord1])[1]}")),
             yaxis = list(title = ylab),
             barmode = 'stack')

    for (i in 2:length(unique(meltdat[, Ord1]))) {
      p1$x$layoutAttrs[[1]][[paste0("xaxis", i)]] = NULL
      p1$x$layoutAttrs[[1]][[paste0("xaxis", i)]]$title <- glue("{Ord1} = {levels(meltdat[, Ord1])[i]}")
    }
    if(!is.null(outfile)){
      htmlwidgets::saveWidget(p1, outfile)
    }
  flog.info('Finish...')
    return(p1)
  }
}
