#' ASVenn
#'
#' Function to create Venn Diagram of shared features
#'
#' @param data a phyloseq object (output from decontam or generate_phyloseq)
#' @param output Output directory
#' @param rank Taxonomy rank to merge features that have same taxonomy at a certain taxonomic rank (among rank_names(data), or 'ASV' for no glom)
#' @param column1 Column name of factor to test (among sample_variables(data))
#' @param subset Subset sample, please provide as c(FACTOR,LEVEL).
#' @param lvls Vector comma separated list levels of factor to print in venn diagram (max. 5).
#' @param krona Krona of exclusive ASV or shared with informed level and others. Must be among levels of column1 argument.
#' @param shared shared [TRUE] or exclusive [FALSE] mode.
#'
#'
#' @return Returns list with venn diagram and table with shared features. Exports a venn diagram with corresponding tabulated file.
#'
#' @importFrom glue glue
#' @importFrom venn venn
#' @export



ASVenn_fun <- function(data = data, output = "./ASVenn/", rank = "ASV",
                            column1 = NULL, subset = "", lvls = NULL, krona = "",
                            shared = TRUE){

  invisible(flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger"))

  if(!is.null(output)){
    if(!dir.exists(output)){
      dir.create(output, recursive=TRUE)
    }
  }


  #Check phyloseq object named data exists
  if(!any(ls()=="data")){
    for(i in ls()){
      fun <- paste("cLS <- class(",i,")")
      eval(parse(text=fun))
      # print(c(i,cLS))
      if(cLS == "phyloseq"){
        fun <- paste("data = ", i)
        eval(parse(text=fun))
      }
    }
  }


  if (is.null(column1)){
    print_help(opt_parser)
    flog.info("You must provide a factor:")
    print(names(sample_data(data)))
    stop()
  }

  #Subset data
  if(subset!=""){
    flog.info('Subset phyloseq object ...')
    args1 <- unlist(strsplit(subset,","))
    fun <- paste("data <- subset_samples(data, ",args1[1]," %in% '",args1[2],"')",sep="")
    eval(parse(text=fun))
    TITRE=paste(column1, args1[1], args1[2], sep="-")
  }else{
    TITRE=paste(column1)
  }

  if(length(lvls)==0){
    flog.error('You must provide levels...')
    stop()
  } else if(length(lvls)>5){
    flog.error('Venn diagram is limited to 5 levels')
    stop()
  } else{
    if(!all(lvls %in% na.omit(levels(as.factor(sample_data(data)[,column1]@.Data[[1]])) ))){
      flog.error('Your levels are not present in metadata...')
      stop()
    }
  }


  if(length(lvls)==0){
    flog.error('You must provide levels...')
    stop()
  } else if(length(lvls)>5){
    flog.error('Venn diagram is limited to 5 levels')
    stop()
  } else{
    if(!all(lvls %in% na.omit(levels(as.factor(sample_data(data)[,column1]@.Data[[1]])) ))){
      flog.error('Your levels are not present in metadata...')
      stop()
    }
  }


  #Nombre d'espèce par matrice
  flog.info('Parsing factor ...')
  level1 <- na.omit(levels(as.factor(sample_data(data)[,column1]@.Data[[1]])) )
  TFdata <- list()
  TFtax <- list()
  if(!is.null(refseq(data, errorIfNULL=FALSE))){
    refseq1 <- as.data.frame(refseq(data)); names(refseq1)="seq"
  }else{refseq1 = NULL; flog.info('No Tree ...')}
  if(rank!="ASV"){
    flog.info(glue::glue("Glom phyloseq object to {rank} rank..."))
    ;data <- tax_glom(data, rank)
  }
  # print(data)

  databak <- data
  for(i in 1:length(lvls)){
    databak -> data
    LOC=as.character(lvls[i])
    print(LOC)
    fun <- paste("data <- subset_samples(data, ",column1," %in% '",LOC,"')",sep="")
    eval(parse(text=fun))

    # print(rank)
    sp_data <- prune_taxa(taxa_sums(data) > 0, data)
    # cat(LOC,ntaxa(sp_data)," ", rank, " \n")
    print(data)
    ttable <- sp_data@tax_table@.Data
    otable <- as.data.frame(otu_table(sp_data))
    # print(nrow(ttable))

    if(!any(rownames(ttable) == rownames(otable))){flog.info("Different order in otu table and tax table");quit()}

    TT = cbind(otable,ttable)

    TFdata[[lvls[i]]] <- TT
    TFtax[[lvls[i]]] <- cbind(row.names(TT), as.character(apply(TT[,colnames(ttable)], 1, paste, collapse=";") ) ) #c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
    row.names(TFtax[[lvls[i]]]) = TFtax[[lvls[i]]][,1]

    # write.table(TT, paste(output,"/otu_table_sp_",LOC,".csv",sep=""), sep="\t", quote=FALSE, col.names=NA)

  }
  # print(names(TFtax))
  # print(TFdata)
  print(str(TFtax))
  # stop()
  ## Venn diag
  flog.info('Defining unique taxa ...')
  alltax <- do.call(rbind, TFtax)
  alltax <- alltax[!duplicated(alltax[,1]),]
  row.names(alltax)=alltax[,1]
  # print(alltax)
  # stop()
  # print(alltax)
  flog.info('Plotting ...')

  # Specific use to screen taxonomic composition of shared taxa...
  if(krona != ""){
    TF <- TFbak
    flog.info('Krona ...')
    env1 <- TF[[krona]]
    others1 <- unique( unlist( TF[level1[level1!=krona]] ) )

    TF2 <- list(env1, others1)
    names(TF2) <- c(krona, "others")
    #Venn 2
    venn.plot <- venn.diagram(TF2, filename = NULL, col = "black",
                              fill = rainbow(length(TF2)), alpha = 0.50,
                              cex = 1.5, cat.col = 1, lty = "blank",
                              cat.cex = 1.8, cat.fontface = "bold",
                              margin = 0.1, main=TITRE, main.cex=2.5,
                              fontfamily ="Arial",main.fontfamily="Arial",cat.fontfamily="Arial") #cat.dist = 0.09,
    venn_tab=paste(output,"/",TITRE,"_",krona, "_kronaVenn.png", sep="")
    png(venn_tab, width=20, height=20, units="cm", res=200)
    grid.draw(venn.plot)
    invisible(dev.off())

    #Krona
    if(shared==TRUE){
      L1 = intersect(env1, others1)
    }else{
      L1 = setdiff(env1, others1)
    }

    subtaxdata=prune_taxa(L1,data)
    ttable <- as.data.frame(subtaxdata@tax_table@.Data)
    fttable = cbind.data.frame(rep(1,nrow(ttable)),ttable)
    fttable$ASV = row.names(fttable)
    flog.info('Write Krona table ...')
    krona_tab=paste(output,"/",TITRE,"_",krona, "_krona.txt", sep="")
    write.table(fttable, krona_tab, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

    flog.info('Generate Krona html ...')
    output <- paste(output,"/",TITRE,"_",krona, "_krona.html", sep = "")
    system(paste("ktImportText", krona_tab, "-o", output, sep = " "))
    # browseURL(output)
  }



  TFbak <- TF <- sapply(TFtax, row.names, simplify = FALSE)
  # print( length(unique(unlist(TF))) )

  if(length(lvls)>5){
    flog.info('Too much levels (max. 5) ...')
    # print(lvls)
    if(is.null(lvls)){
      flog.info('Selecting 5 first levels ...')
      TF <- TF[1:5]
      res1 = VENNFUN(TF = TF, mode=1, TITRE = TITRE, output = output, refseq1 = refseq1, alltax = alltax)
    }else{
      # print(names(TF))
      flog.info(glue('Selecting {lvls} ...'))
      LVLs <- unlist(strsplit(lvls,","))
      TF <- TF[match(LVLs, names(TF))]
      if(length(TF) <= 5){
        flog.info(glue('mode 1 ...'))
        res1 = VENNFUN(TF = TF, mode = 1, TITRE = TITRE, output = output, refseq1 = refseq1, alltax = alltax)
      }else{
        flog.info(glue('mode 2 ...'))
        res1 = VENNFUN(TF = TF, mode = 2, TITRE = TITRE, output = output, refseq1 = refseq1, alltax = alltax)
      }
    }
  } else {
      print("< 5 levels")
    if(!is.null(lvls)){
      flog.info(glue('Selecting {lvls} ...'))
      LVLs <- unlist(strsplit(lvls,","))
      TF <- TF[match(LVLs, names(TF))]
    }
    res1 = VENNFUN(TF = TF, mode = 1, TITRE = TITRE, output = output, refseq1 = refseq1, alltax = alltax)
  }

  flog.info('End ...')
#save(list = ls(all.names = TRUE), file = "debug_asvenn.rdata", envir = environment())
  return(res1)
}


#' VENNFUN
#'
#' Plotting Venn Function and table of shared feature generation
#'
#' @param TF list containing vectors to plot.
#' @param mode = 1: length(TF)<=5, mode = 2 5<length(TF)<7
#' @param TITRE
#' @param output
#' @param refseq1
#' @param alltax
#'
#'
#' @return Exports a venn diagram with corresponding tabulated file.
#'
#' @importFrom glue glue
#' @importFrom venn venn
#' @export


VENNFUN <- function(TF = TF, mode = 1, TITRE = TITRE, output = "./", refseq1 = NULL, alltax=NULL){
  if(mode==1){
    venn::venn(TF, zcol = rainbow(7), ilcs = 2, sncs = 2) #, col=rainbow(7)
    venn.plot <- recordPlot()
    invisible(dev.off())

    if(!is.null(output)){
      png(paste(output,'/',TITRE,'_venndiag.png',sep=''), width=20, height=20, units="cm", res=200)
      replayPlot(venn.plot)
      dev.off()
    }

    print("plotOK")


    ov <- calculate.overlap(TF)
    # print(sapply(ov, length))

    flog.info('Calculating lists ...')
    uniqTax = TABf = unique(do.call(c,TF))
    for (j in 1:length(TF)){
      TABtest = TF[[j]]
      TABtest_filt=rep(0, length(uniqTax))
      for (i in 1:length(uniqTax)) {
        featureI = uniqTax[i]
        res=grep( paste('^',featureI,'$', sep="" ) , TABtest)
        if(length(res)>0){TABtest_filt[i]=length(res)
        }
      }
      TABf=cbind.data.frame( TABtest_filt, TABf )
      names(TABf)[1] = names(TF)[j]
    }

    if(!is.null(alltax)){
      if(!is.null(refseq1)){
        TABf <- cbind(TABf,alltax[as.character(TABf$TABf),2], refseq1[as.character(TABf$TABf),])
        names(TABf) <- c(rev(names(TF)), "ASV", "taxonomy", "seq")
      }else{
        TABf <- cbind(TABf,alltax[as.character(TABf$TABf),2])
        names(TABf) <- c(rev(names(TF)), "ASV", "taxonomy")
      }
    }

    if(!is.null(output)){
      write.table(TABf, paste(output,"/",TITRE,"_venn_table.csv",sep=""), sep="\t", quote=FALSE, row.names=FALSE)
    }
  } else if(mode == 2){ # more than 5 environments

    venn.plot <- venn::venn(TF, zcol = rainbow(7), ilcs = 2, sncs = 2) #, col=rainbow(7)
    venn.plot <- recordPlot()
    invisible(dev.off())

    if(!is.null(output)){
      png(paste(output,'/',TITRE,'_venndiag.png',sep=''), width=20, height=20, units="cm", res=200)
      replayPlot(venn.plot)
      dev.off()
    }

    ENVS = na.omit(names(TF)[1:7])  #maximum 7
    Tabf <- NULL; Tab1 <- NULL
    # Exclusive
    for(i in ENVS){
      tt = c(i, ENVS[ENVS != i])
      # print(tt)
      yy = Reduce(setdiff, TF[tt])  # setdiff(setdiff(tt[1], tt[2]), tt[3] )
      # print(length(yy))
      if(!is.null(alltax)){
        Tab1 <- cbind(yy, rep(i, length(yy)), alltax[yy,1])
      }else{
        Tab1 <- cbind(yy, rep(i, length(yy)))
      }
      Tabf <- rbind(Tabf, Tab1)
    }
    #Core
    yy <- Reduce(intersect, TF[ENVS]) #maximum 7
    if(!is.null(alltax)){
      Core <- cbind(yy, rep("core", length(yy)), alltax[yy,1])
    }else{
      Core <- cbind(yy, rep("core", length(yy)))
    }
    TABf <- as.data.frame(rbind(Core, Tabf))
    names(TABf) = paste("V", 1:ncol(TABf), sep="")
    if(!is.null(output)){
      write.table(TABf, paste(output,"/",TITRE,"_venn_table.csv",sep=""), sep="\t", quote=FALSE, row.names=FALSE)
    }
  }

  LL=list()
  LL$venn_plot = venn.plot
  LL$TABf = TABf
  return(LL)
}
