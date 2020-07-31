#' Diversity Alpha
#'
#' Provides boxplots with multiples diversity indices and statistical tests like ANOVA with post hoc test and non parametric Wilcoxon tests.
#'
#' @param data output from decontam or generate_phyloseq
#' @param output Output directory
#' @param column1 Column name of first factor to test (covariable in ANOVA).
#' @param column2 Column name of second factor to test (last factor in ANOVA).
#' @param column3 Column name of subjects in case of repeated mesures (mixed model)
#' @param supcovs One or more supplementary covariables to test in anova (provided as vector)
#' @param measures Diversity indices (provided as vector)
#'
#' @return Export plots and tests in the output directory.
#'
#' @import phyloseq
#' @import ggplot2
#' @importFrom nlme lme
#' @importFrom glue glue
#'
#' @export


# Decontam Function

diversity_alpha_fun <- function(data = data, output = "./plot_div_alpha/", column1 = "", column2 = "",
                                column3 = "", supcovs = "", measures = c("Observed","Shannon","Simpson","InvSimpson")){


  if(!dir.exists(output)){
    dir.create(output, recursive=TRUE)
  }


  ## Gestion des NA
  if(!all(is.na(sample_data(data)[,column1]))){
    fun <- paste("data <- subset_samples(data, !is.na(",column1,"))",sep="")
    eval(parse(text=fun))

    if(column2 != ''){
      fun <- paste("data <- subset_samples(data, !is.na(",column2,"))",sep="")
      eval(parse(text=fun))
    }

    if(column3 != ''){
      fun <- paste("data <- subset_samples(data, !is.na(",column3,"))",sep="")
      eval(parse(text=fun))
    }

    #alpha diversité tableau

    resAlpha = list()
    flog.info('Alpha diversity tab ...')
    alpha.diversity <- estimate_richness(data, measures = c("Observed","Shannon","Simpson","InvSimpson") )
    # row.names(alpha.diversity) <- gsub("X","",row.names(alpha.diversity))
    resAlpha$alphatable = alpha.diversity
    write.table(alpha.diversity,paste(output,'/alphaDiversity_table.csv',sep=''), sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
    flog.info('Done.')


    #alpha diversité graphique
    alphaPlot <- function() {
      flog.info('Plotting ...')
      if(column2 == ''){
        p <- plot_richness(data,x=column1, color=column1, measures=measures)
      } else{
        p <- plot_richness(data,x=column1, color=column2, measures=measures)
      }
      p$layers <- p$layers[-1]
      p <- p + ggtitle('Alpha diversity indexes') +  geom_boxplot(position = position_dodge(width = 0.5),alpha = 0.7, outlier.shape = NA) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=11), plot.title = element_text(hjust = 0.5)) + theme_bw()
      flog.info('Done.')
      return(p)
    }
    p <- alphaPlot()
    if(length(measures) == 1){
      ggplotly(p)
    }else{
      plot(p)
    }
    resAlpha$plot = p
    ggsave(paste(output,'/alpha_diversity.png',sep=''), plot=p, height = 15, width = 30, units="cm")


    anova_data <- cbind(sample_data(data), alpha.diversity)
    anova_data$Depth <- sample_sums(data)
    if(length(levels(as.factor(anova_data[,column1])))>1){
      flog.info('ANOVA ...')
      # variables <- paste(sep=" + ", "Depth", var1)
      sink(paste(output,'/all_ANOVA.txt', sep=''), split = FALSE)
      for (m in measures){

        cat(paste("\n\n############\n",m,"\n############\n"))
        if(supcovs != ""){
          COVS <- sapply(strsplit(supcovs,","), '[')
          print(COVS)
          if(column2 != ''){
            f <- paste(m," ~ ", "Depth + ", paste(COVS, collapse="+"), "+", column1," + ",column2)
            anova_data$fact1 <- paste( anova_data[,column1],  anova_data[,column2], sep="_")
          } else {
            f <- paste(m," ~ ", "Depth + ", paste(COVS, collapse="+"), "+", column1)
          }
        } else {
          if(column2 != ''){
            f <- paste(m," ~ ", "Depth + ", column1," + ",column2)
            anova_data$fact1 <- paste( anova_data[,column1],  anova_data[,column2], sep="_")
          } else {
            f <- paste(m," ~ ", "Depth + ", column1)
          }
        }

        cat("############\nANOVA + pairwise wilcox test\n")
        print(f)
        anova_res1 <- aov( as.formula(paste(f)), anova_data)
        # print(anova_data)
        # write.table(anova_data, paste(output,"/anovatable.csv", sep=""), sep="\t", row.names=FALSE)
        # print(anova_res1)
        res1 <- summary(anova_res1)
        print(res1)

        resAlpha$anova = res1

        # # post hoc test  commented du to conflict between LSD.test() and DESeq() function. #' @importFrom agricolae LSD.test
        # cat("############\npost hoc LSD.test\n")
        # if(column2 == ''){
        #   fun = glue("lsd1 <- LSD.test(anova_res1, '{column1}', p.adj='fdr')")
        #   eval(parse(text = fun))
        #   print(lsd1)
        # } else {
        #   fun = glue("lsd1 <- LSD.test(anova_res1, '{column2}', p.adj='fdr')")
        #   eval(parse(text = fun))
        #   print(lsd1)
        # }

        cat(paste("\n##pvalues of pairwise wilcox test on ", m, "with FDR correction \n"), sep="")
        fun <- paste("wilcox_res1 <- pairwise.wilcox.test(anova_data$",m,", anova_data[,column1], p.adjust.method='fdr')", sep="")
        eval(parse(text = fun))
        print(round(wilcox_res1$p.value,3))

        if(column2 != ''){
          cat(paste("\n##pvalues of pairwise wilcox test on ", m, "with FDR correction \n"), sep="")
          fun <- paste("wilcox_res1 <- pairwise.wilcox.test(anova_data$",m,", anova_data[,column2], p.adjust.method='fdr')", sep="")
          eval(parse(text = fun))
          print(round(wilcox_res1$p.value,3))

          cat("\n\n#######################\n")
          cat(paste("##pvalues of pairwise wilcox test on ", m, " with collapsed factors (no correction)\n"), sep="")
          fun <- paste("wilcox_res <- pairwise.wilcox.test(anova_data$",m,", anova_data$fact1, p.adjust.method='none')", sep="")
          eval(parse(text = fun))
          print(round(wilcox_res$p.value,3))
        }

        resAlpha$wilcox = round(wilcox_res$p.value,3)

        if(column3 != ''){
          cat("\n############\nANOVA repeated measures\n")
          print(paste(f,"+ Error(",column3,")"))
          anova_res <- aov( as.formula(paste(f,"+ Error(",column3,")")), anova_data)
          res <- summary(anova_res)
          print(res)

          resAlpha$anovarepeat = res

          cat("\n############\nMixed effects models (nlme::lme)\n")
          print(f)
          print(paste("lme1 = anova(lme(as.formula(",f,"), random = ~1 | ",column3,", data = anova_data, method = 'ML') )", sep=""))
          fun <- paste( "lme1 = anova(lme(as.formula(",f,"), random = ~1 | ",column3,", data = anova_data, method = 'ML') ) ", sep="")
          eval(parse(text=fun))
          print(lme1)

          resAlpha$mixedeffect = lme1
        }

      }
      sink()

    }else{flog.info('Factor with less than 2 levels, no test ...'); print(levels(as.factor(anova_data[,column1])))}
    flog.info('Done.')


  }else{flog.info(paste(column1, 'is all NA.'))}

  flog.info('Finish.')

  return(resAlpha)

}
