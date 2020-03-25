#' Barplots community
#'
#'
#' @param data output from decontam or generate_phyloseq
#' @param output Output directory
#' @param bar Plotting bar plot with raw reads number.
#' @param compo Plotting with relative composition.
#' @param column1 Column name of factor used to sort sample
#' @param column2 Column name of factor used to split barplot
#' @param sname Change sample.id by the corresponding factor levels in graph.
#' @param num Number of top taxon to display.
#' @param rare Column name for splitting rare curves.
#' @param rank Taxonomic rank name. You can provide multiple ranks seperated by comma.
#'
#' @return Export barplots in an html file.
#'
#' @import phyloseq
#' @import ggplot2
#' @import gridExtra
#' @import grid
#' @importFrom microbiome aggregate_top_taxa
#' @importFrom microbiome plot_composition
#' @importFrom microbiome transform
#' @importFrom ggpubr as_ggplot
#' @importFrom ggpubr get_legend
#' @importFrom plotly ggplotly
#' @import viridis
#' @importFrom phyloseq.extended ggrare
#'
#'
#'
#' @export


# Decontam Function

bars_fun <- function(data = data, bar = TRUE, compo1 = TRUE, output = "./plot_bar/", column1 = "", column2 = "",
                     sname = FALSE, num = 10, rare = NULL, rank = "Genus"){

  out1 <- paste(getwd(),'/',output,'/',sep='')
  if(!dir.exists(output)){

    dir.create(out1, recursive=TRUE)
  }else{
    unlink(out1, recursive=TRUE)
    dir.create(out1, recursive=TRUE)
  }

    taxonomy <- sapply(strsplit(rank,","), '[')


    #Gestion NA
    fun <- paste("data <- subset_samples(data, !is.na(",column1,"))",sep="")
    eval(parse(text=fun))

    # Bars
    rmd_data=list()
    if(bar==TRUE){

      for(i in 1:length(taxonomy)){
        j <- taxonomy[i]; print(j)

        psobj.top <- aggregate_top_taxa(data, j, top = num)
        tn = taxa_names(psobj.top)
        tn[tn=="Other"] = tn[length(tn)]
        tn[length(tn)] = "Other"
        if(column1 != '' & column2 == ''){
          flog.info('Plotting bar (%s)...',j)
          if(sname == TRUE){
            plot.composition.COuntAbun <- plot_composition(psobj.top, x.label = column1, sample.sort=column1, otu.sort=tn)
          } else{
            plot.composition.COuntAbun <- plot_composition(psobj.top, x.label = "sample.id", sample.sort=column1, otu.sort=tn)
          }

          plot.composition.COuntAbun <- plot.composition.COuntAbun + theme(legend.position = "bottom") +
            theme_bw() + scale_fill_viridis(discrete = TRUE, direction=-1) +
            theme(axis.text.x = element_text(angle = 90)) +
            ggtitle("Raw abundance") + theme(legend.title = element_text(size = 18))


          flog.info('Done.')
        }
        else if(column1 != '' & column2 != ''){
          flog.info('Plotting bar (%s)...',j)

          sdata = psobj.top@sam_data@.Data
          names(sdata) = psobj.top@sam_data@names
          FACT1=sdata[[column2]]
          (lvls = levels(as.factor(sdata[[column2]])))
          LL <- list()
          for(i in 1:length(lvls)){
            fun  <- paste("ppp <- subset_samples(psobj.top, ",column2," %in% '",lvls[i],"')",sep="")
            eval(parse(text=fun))
            tn = taxa_names(ppp)
            tn[tn=="Other"] = tn[length(tn)]
            tn[length(tn)] = "Other"
            if(sname){
              p1 <- plot_composition(ppp, x.label = column1, verbose=TRUE,sample.sort=column1, otu.sort=tn)
            } else{
              p1 <- plot_composition(ppp, x.label = "sample.id", verbose=TRUE,sample.sort=column1, otu.sort=tn)
            }


            p1 <- p1 + theme(legend.position = "bottom") +
              theme_bw() + scale_fill_viridis(discrete = TRUE, direction=-1) +
              theme(axis.text.x = element_text(angle = 90), legend.title = element_text(size = 18)) +
              ggtitle(paste("Raw abundance",column2,"=", lvls[i]))

            ###ici legende separee
            if(i < length(lvls)){
              p2 <- p1 + theme(legend.position = "none")
              LL[[i]] <- p2
            }else{
              p2 <- p1 + theme(legend.position = "none")
              LL[[i]] <- p2
              leg = ggpubr::get_legend(p1)
              LL[[i+1]] <- ggpubr::as_ggplot(leg)
            }
          }

          plot.composition.COuntAbun <- p3 <- gridExtra::grid.arrange(grobs=LL, ncol=length(lvls)+1)
          flog.info('Done.')
        }

        fun <- paste('rmd_data$bars$',j,' <- plot.composition.COuntAbun',sep='')
        eval(parse(text=fun))
      }
    }




    # Composition
    compo <- function (rankList, psobj, var = 'sample.id', col2='') {
      if(col2 != ''){
        sdata = psobj@sam_data@.Data
        names(sdata) = psobj@sam_data@names
        FACT1=sdata[[col2]]
        lvls = levels(as.factor(sdata[[col2]]))

        for(i in 1:length(rankList)){
          j <- rankList[i]
          LL <- list()
          for(k in 1:length(lvls)){
            LVL=lvls[k]
            fun  <- paste("ppp <- subset_samples(psobj, ",column2," %in% '",LVL,"')",sep="")
            eval(parse(text=fun))
            flog.info('Plotting %s composition (%s)...',var, j)
            psobj.fam <- aggregate_top_taxa(ppp, j, top = num)
            psobj.rel <-  transform(psobj.fam, "compositional")
            tn = taxa_names(psobj.rel)
            tn[tn=="Other"] = tn[length(tn)]
            tn[length(tn)] = "Other"
            p1 <- plot_composition(psobj.rel, x.label = column1, sample.sort=column1, otu.sort=tn) +
              theme() + theme_bw() + scale_fill_viridis(discrete = TRUE, direction=-1) +
              theme(axis.text.x = element_text(angle = 90), legend.title = element_text(size = 18)) +
              ggtitle(paste("Relative abundance",column2,"=", LVL))

            if(k < length(lvls)){
              p2 <- p1 + theme(legend.position = "none")
              LL[[k]] <- p2
            }else{
              p2 <- p1 + theme(legend.position = "none")
              LL[[k]] <- p2
              leg = ggpubr::get_legend(p1)
              LL[[k+1]] <- ggpubr::as_ggplot(leg) + theme( plot.background = element_blank())
            }

          }
          p3 <- gridExtra::grid.arrange(grobs=LL, ncol=length(lvls)+1)
          flog.info('Done.')

          # fun <- paste('rmd_data$compo$',j,'[["',column2,'_',LVL,'"]]',' <- p3',sep='')
          fun <- paste('rmd_data$compo$',j,' <- p3',sep='')
          eval(parse(text= fun))
          flog.info('Done.')
        }


      }else{

        for(i in 1:length(rankList)){
          j <- rankList[i]
          flog.info('Plotting %s composition (%s)...',var, j)
          psobj.fam <- aggregate_top_taxa(psobj, j, top = num)
          psobj.rel <-  transform(psobj.fam, "compositional")
          tn = taxa_names(psobj.rel)
          tn[tn=="Other"] = tn[length(tn)]
          tn[length(tn)] = "Other"

          p1 <- plot_composition(psobj.rel, x.label = column1, sample.sort=column1, otu.sort=tn) +
            theme() + theme_bw() + scale_fill_viridis(discrete = TRUE, direction=-1) +
            theme(axis.text.x = element_text(angle = 90), legend.title = element_text(size = 18)) +
            ggtitle("Relative abundance")
          fun <- paste('rmd_data$compo$',j,' <- p1',sep='')
          eval(parse(text= fun))
          flog.info('Done.')
        }
      }

      return(rmd_data)
    }
    #

    if(!is.null(rare)){
      flog.info('Plotting rarefaction ...')
      GROUPE=rare
      plot_rare <- ggrare(data, step = 100, color = rare, plot = FALSE)
      plot_rare <- plot_rare + facet_wrap(GROUPE, ncol = 4) + theme_bw()
      rmd_data$rare <- ggplotly(plot_rare)
      flog.info('Done.')
    }

    #Coupage selon le facteur 2
    if(compo1==TRUE){
      if(column1 != '' & column2 != ''){
        vector <- levels(data.frame(sample_data(data)[,column2])[,1])
        rmd_data <- compo(taxonomy,data,col2=column2)

      }else{
        rmd_data <- compo(taxonomy,data,column1)
      }
    }


    # Generating rmd template for report
    # PATHanomaly="/home/erifa/Repository/LRF/anomaly/"
    sink(paste(out1,'/bars2.Rmd', sep=""))
    cat("---
title: Bar plot
fig_width: 24
params:
  rmd_data: p
  col1: col1
---

```{r message=FALSE, warning=FALSE, include=FALSE, results='hide'}
rmd_data <- params$rmd_data

```

```{r hold=TRUE, echo=FALSE, comment = FALSE, message= FALSE, warning = FALSE, results='asis',fig.keep='all', fig.align='left', fig.width = 10, fig.height = 10}
if('rare' %in% names(rmd_data)){
  cat('# Rarefaction plot\\n')
  rmd_data$rare
}
```

```{r hold=TRUE, echo=FALSE, comment = FALSE, message= FALSE, warning = FALSE, results='asis',fig.keep='all', fig.align='left', fig.width = 20, fig.height = 10}
if('bars' %in% names(rmd_data)){
  cat('# Plot raw value composition')
  cat('\\n')
}
```

")

    for(Nplot in names(rmd_data$bars)){
      if(Nplot %in% names(rmd_data$bars)){
        cat(paste("
```{r hold=TRUE, echo=FALSE, comment = FALSE, message= FALSE, warning = FALSE, results='asis',  fig.keep='all', fig.align='left', fig.width = 20, fig.height = 10}
    cat('\\n')
    cat('## ",Nplot,"')
    cat('\\n')
    grid.draw(rmd_data$bars[['",Nplot,"']])
```
  ", sep=""))
      }

    }


    cat("
```{r hold=TRUE, echo=FALSE, comment = FALSE, message= FALSE, warning = FALSE, results='asis',fig.keep='all', fig.align='left', fig.width = 20, fig.height = 10}
if('compo' %in% names(rmd_data)){
  cat('# Composition plot')
  cat('\\n')
}
```
")

    for(Nplot in names(rmd_data$compo)){
      if(Nplot %in% names(rmd_data$compo)){
        cat(paste("
```{r hold=TRUE, echo=FALSE, comment = FALSE, message= FALSE, warning = FALSE, results='asis',  fig.keep='all', fig.align='left', fig.width = 20, fig.height = 10}
    cat('\\n')
    cat('## ",Nplot,"')
    cat('\\n')
    grid.draw(rmd_data$compo[['",Nplot,"']])
```
  ", sep=""))
      }

    }

    sink()
    # cat(paste('## ',",Nplot,",sep=''))

    rmarkdown::render(paste(out1,'/bars2.Rmd', sep=""),params= list('rmd_data' = rmd_data, 'col1' = column1),output_file=paste(out1,'/','bars.html',sep=''))  ## determiner automatiquement le path du md
    flog.info('Finish.')




  }




