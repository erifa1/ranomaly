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
#' @param shared shared [TRUE] or exclusive [FALSE] mode for krona plot.
#' @param verbose Verbose level. (1: quiet, 2: print infos, 3: print infos + debug)
#' @param ggplotmode if TRUE plot the Venn diagram using ggplot
#'
#'
#' @return Returns list with venn diagram and table with shared features. Exports a venn diagram with corresponding tabulated file.
#'
#' @importFrom glue glue
#' @importFrom venn venn
#' @export



ASVenn_fun <- function(data = data, output = "./ASVenn/", rank = "ASV",
                            column1 = NULL, subset = NULL, lvls = NULL, krona = "",
                            shared = TRUE, verbose = 2, ggplotmode = FALSE){

  invisible(flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger"))

  if(verbose == 3){
    invisible(flog.threshold(DEBUG))
  }
  if(verbose == 2){
    invisible(flog.threshold(INFO))
  }
  if(verbose == 1){
    invisible(flog.threshold(ERROR))
  }



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
  if(!is.null(subset)){
    flog.info('Subset phyloseq object ...')
    args1 <- unlist(strsplit(as.character(subset),","))
    fun <- paste("data <- subset_samples(data, ",args1[1]," %in% '",args1[2],"')",sep="")
    eval(parse(text=fun))
    TITRE=paste(column1, args1[1], args1[2], sep="-")
  }else{
    TITRE=paste(column1)
  }

  if(length(lvls)==0 & length(levels(as.factor(sample_data(data)[,column1]@.Data[[1]]))) > 7 ){
    flog.error('More than 7 levels in the provided factor, you must specify levels...')
    stop()
  } else if(length(lvls)==0 & length(levels(as.factor(sample_data(data)[,column1]@.Data[[1]]))) <= 7 ){
    lvls <- unique(na.omit(sample_data(data)[,column1]@.Data[[1]]))
  } else if(length(lvls)>7){
    flog.error('Venn diagram is limited to 7 levels')
    stop()
  } else {
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
    flog.info(LOC)
    fun <- paste("data <- subset_samples(data, ",column1," %in% '",LOC,"')",sep="")
    eval(parse(text=fun))

    # print(rank)
    sp_data <- prune_taxa(taxa_sums(data) > 0, data)
    # cat(LOC,ntaxa(sp_data)," ", rank, " \n")
    # print(data)
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

  ## Venn diag
  flog.info('Defining unique taxa ...')
  alltax <- do.call(rbind, TFtax)
  alltax <- alltax[!duplicated(alltax[,1]),]
  row.names(alltax)=alltax[,1]

  flog.info('Plotting ...')
  # Specific use to screen taxonomic composition of shared taxa...
  TFbak <- TF <- sapply(TFtax, row.names, simplify = FALSE)
  names(TFbak) = names(TF) = lvls
  # print( length(unique(unlist(TF))) )


  if(krona != ""){
    flog.info('Krona ...')
    env1 <- TF[[krona]]
    others1 <- unique( unlist( TF[level1[level1!=krona]] ) )

    TF2 <- list(env1, others1)
    names(TF2) <- c(krona, "others")
    #Venn 2
    # pdf(NULL)
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
    outputkrona <- paste(output,"/",TITRE,"_",krona, "_krona.html", sep = "")
    system(paste("ktImportText", krona_tab, "-o", outputkrona, sep = " "))
    # browseURL(output)
  }

  res1 = VENNFUN(TF = TF, TITRE = TITRE, output = output, refseq1 = refseq1, alltax = alltax)

  flog.info('End ...')
  return(res1)
}


#' VENNFUN
#'
#' Plotting Venn Function and table of shared feature generation
#'
#' @param TF list containing vectors to plot.
#' @param TITRE Plot title.
#' @param output Output path.
#' @param refseq1 Reference sequences.
#' @param alltax Taxonomy table.
#' @param ggplotmode if TRUE plot the Venn diagram using ggplot
#'
#'
#' @return Exports a venn diagram with corresponding tabulated file.
#'
#' @importFrom glue glue
#' @importFrom venn venn
#' @importFrom qdapTools mtabulate
#' @export


VENNFUN <- function(TF = TF, TITRE = TITRE, output = "./", refseq1 = NULL, alltax=NULL, ggplotmode = FALSE){

    venn.plot <- venn::venn(TF, zcol = rainbow(7), ilcs = 2, sncs = 2, ggplot = FALSE, ilabels = "counts") #, col=rainbow(7)
    venn.plot <- recordPlot()

    invisible(dev.off())

    if(!is.null(output)){
      png(paste(output,'/',TITRE,'_venndiag.png',sep=''), width=20, height=20, units="cm", res=200)
      replayPlot(venn.plot)
      dev.off()
    }

    # Generating table of shared features
    v.table <- tibble::as_tibble(t(qdapTools::mtabulate(TF)), rownames = "taxa")
    colnames(alltax) <- c("taxa", "taxonomy")
    alltaxDF <- as.data.frame(alltax)
    refseq1DF <- tibble::as_tibble(refseq1, rownames = "taxa")
    v.table <- dplyr::left_join(v.table, alltaxDF, by = 'taxa')
    TABf<- dplyr::left_join(v.table, refseq1DF, by = 'taxa')

    if(!is.null(output)){
      write.table(TABf, paste(output,"/",TITRE,"_venn_table.csv",sep=""), sep="\t", quote=FALSE, row.names=FALSE)
    }


  LL=list()
  LL$venn_plot <- venn.plot
  LL$TABf <- TABf

  flog.info('Finish ...')
  return(LL)
}
