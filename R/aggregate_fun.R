#' Aggregate
#'
#' Aggregate results of the tree differential analysis methods in one single file.
#'
#' @param data a phyloseq object (output from decontam or generate_phyloseq)
#' @param output Output directory
#' @param metacoder Path to the metacoder CSV file
#' @param deseq Path to deseq results folder
#' @param mgseq Path to metagenomeseq results folder
#' @param column1 Column name of factor to test (among sample_variables(data))
#' @param column2 Column name on which table were splitted
#' @param verbose Verbose level. (1: quiet, 2: print infos, 3: print infos + debug)
#' @param rank Taxonomy rank to merge features that have same taxonomy at a certain taxonomic rank (among rank_names(data), or 'ASV' for no glom)
#' @param comp Comparison to test. Comma separated and comparisons are informed with a tilde (A~C,A~B,B~C). If empty, test all combination.
#' @param returnval Boolean for function to return values.
#'
#' @return Export final CSV files, barplot with top significant ASV and Venn Digramm.
#'
#' @import phyloseq
#' @import VennDiagram
#' @import ggplot2
#'
#' @export


# Decontam Function

aggregate_fun <- function(data = data, metacoder = NULL, deseq = NULL, mgseq = NULL, output = "./aggregate_diff/",
                          column1 = NULL, column2 = NULL, verbose = 1, rank = "Species", comp = "", returnval = TRUE){

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

  if(!dir.exists(output)){
    dir.create(output,recursive=T )
  }

  flog.info('Glom Tax.')
  data.glom <- tax_glom(data, taxrank=rank)

  if(comp == ''){
    fun <- paste('combinaisons <- combn(na.omit(unique(sample_data(data)$',column1,')),2)',sep='')
    eval(parse(text=fun))
  }else{
    comp_list <- unlist(strsplit(comp,","))
    combinaisons <- matrix(, nrow = 2, ncol = length(comp_list))
    for (i in 1:length(comp_list)){
      tmp <- unlist(strsplit(comp_list[i],"~"))
      combinaisons[1,i] <- tmp[1]
      combinaisons[2,i] <- tmp[2]
    }
  }

  #Metacoder table
  otable <- otu_table(data)
  ttax <- tax_table(data)
  ssample <- as.matrix(sample_data(data))
  sseq <- refseq(data)
  if(file.exists(paste(metacoder))){
    mcoderTab <- read.table(paste(metacoder), h=TRUE)
  }

  outF = list()
  TABfinal <- data.frame()
  for (col in (1:ncol(combinaisons))){
    flog.info(paste('Combinaison ',combinaisons[1,col], ' ' , combinaisons[2,col],sep=''))
    comp1 = paste(combinaisons[1,col], '_vs_' , combinaisons[2,col],sep='')
    flog.info('Deseq2.')
    if(file.exists(paste(deseq,"/signtab_",column1,"_",combinaisons[1,col],"_vs_",combinaisons[2,col],".csv",sep=""))){
      deseqT <- read.table(paste(deseq,"/signtab_",column1,"_",combinaisons[1,col],"_vs_",combinaisons[2,col],".csv",sep=""), h=TRUE,sep="\t")
    } else{
      flog.warn('File does not exists.')
      deseqT <- data.frame()
    }
    flog.info('MetagenomeSeq.')
    if(file.exists(paste(mgseq,"/signtab_",column1,"_",combinaisons[1,col],"_vs_",combinaisons[2,col],".csv",sep=""))){
      mgseqT <- read.table(paste(mgseq,"/signtab_",column1,"_",combinaisons[1,col],"_vs_",combinaisons[2,col],".csv",sep=""), h=TRUE,sep="\t")
    } else{
      if(file.exists(paste(mgseq,"/signtab_",column1,"_",combinaisons[2,col],"_vs_",combinaisons[1,col],".csv",sep=""))){
        mgseqT <- read.table(paste(mgseq,"/signtab_",column1,"_",combinaisons[2,col],"_vs_",combinaisons[1,col],".csv",sep=""), h=TRUE,sep="\t")
      }
      else{
        flog.warn('File does not exists.')
        mgseqT <- data.frame()
      }
    }
    flog.debug(pander::pander(mgseqT, split.tables=2000))
    # print(head(mgseqT))
    flog.info('Metacoder.')
    if(file.exists(paste(metacoder))){
      mcoderT <- mcoderTab[mcoderTab$treatment_1 == as.character(combinaisons[1,col]) & mcoderTab$treatment_2 == as.character(combinaisons[2,col]),]
      if(nrow(mcoderT) == 0){
        flog.warn("Metacoder table empty...")
        mcoderT <- mcoderTab[mcoderTab$treatment_2 == as.character(combinaisons[1,col]) & mcoderTab$treatment_1 == as.character(combinaisons[2,col]),]
      }
      mcoderT$wilcox_p_value[is.na(mcoderT$wilcox_p_value)] = 1
      mcoderTsignif <- mcoderT[mcoderT$wilcox_p_value <= 0.05,]
    }else{
      flog.warn('File does not exists.')
      mcoderT <- data.frame()
      mcoderTsignif <- data.frame()
    }
    # print(mcoderTsignif)

    flog.info("Retrieving IDS of differential features...")
    if(length(deseqT) > 0){
      deseq_ids <- as.character(deseqT[deseqT$padj <= 0.05 & !is.na(deseqT$padj),1])
    }else{
      deseq_ids <- NULL
    }
    if(length(mgseqT) > 0){
      mgseq_ids <- as.character(mgseqT[,1])
    }else{
      mgseq_ids <- NULL
    }
    if(length(mcoderTsignif) > 0){
      mcoder_ids <- as.character(mcoderTsignif$otu_id)
    }else{
      mcoder_ids <- NULL
    }

    # Unique IDS
    flog.info('Filtering unique IDS...')
    ListAllOtu = unique(c(deseq_ids,mgseq_ids,mcoder_ids))

    # flog.info('Deseq2.')
    # print(deseq_ids)
    # flog.info('MetagenomeSeq.')
    # print(mgseq_ids)
    # flog.info('Metacoder.')
    # print(mcoder_ids)
    # flog.info('All.')
    # print(ListAllOtu)

    #Construction de la table
    flog.info('Building table...')
    col_comp = rep(comp1, length(ListAllOtu))
    if(!is.null(column2)){
      col_env = rep(unique(ssample[,column2]), length(ListAllOtu))
      TABf = cbind.data.frame(ListAllOtu,col_env,col_comp)
    } else {TABf = cbind.data.frame(ListAllOtu, col_comp)}

    TF = list(x1=deseq_ids, x2=mgseq_ids, x3=mcoder_ids)
    names(TF) <- c("DESeq", "metagenomeSeq", "metacoder")
    # print(TF)
    if(length(unlist(TF)) == 0){next}

    # Test si chaque ASV est diff dans les méthdes.
    for (j in 1:length(TF)){
      TABtest = TF[[j]]
      # TABtest=gsub("\\[|\\]", "", TF[[j]]) # cherche les ASVids

      TABtest_signif=rep(0, length(ListAllOtu))
      for (i in 1:length(ListAllOtu)) {
        featureI = ListAllOtu[i]
        res=grep( paste('^',featureI,'$', sep="" ) , TABtest)
        #print(res)
        if(length(res)>0){TABtest_signif[i]=length(res)
        #print(c(featureI, TABtest[res]))
        }
      }

      TABf=cbind.data.frame(TABf, TABtest_signif)
      names(TABf)[ncol(TABf)] = names(TF)[j]
    }

    TABfbak0 <- TABf
    # add new columns, sumMethods, DeseqLFC, Mean Relative Abundance (TSS) condition 1 & 2
    if(nrow(deseqT)==0){
      TABf <- cbind(TABf, sumMethods = apply(TABf[3:5], 1, sum, na.rm=TRUE),
                    DESeqLFC = rep(NA, nrow(TABf)),
                    absDESeqLFC = rep(NA, nrow(TABf)))
    }else{
      row.names(deseqT) = deseqT[,1]
      TABf <- cbind( TABf, sumMethods = apply(TABf[3:5], 1, sum, na.rm=TRUE),
                     DESeqLFC = deseqT[as.character(TABf[,1]),"log2FoldChange"],
                     absDESeqLFC = abs(deseqT[as.character(TABf[,1]),"log2FoldChange"]) )
    }
    # clr = function(x){log(x+1) - rowMeans(log(x+1))}
    # otableNORM <- clr(otable)
    normf = function(x){ x/sum(x) }
    data.norm <- transform_sample_counts(data.glom, normf)
    otableNORM <- otu_table(data.norm)

    Gtab <- cbind(as.data.frame(ssample), t(otableNORM))
    MeanRelAbcond1=NULL
    for(i in TABf$ListAllOtu){
      tt=mean(Gtab[Gtab[,column1]==combinaisons[1,col],i], na.rm=TRUE)
      MeanRelAbcond1=c(MeanRelAbcond1,tt)
    }
    MeanRelAbcond2=NULL
    for(i in TABf$ListAllOtu){
      tt=mean(Gtab[Gtab[,column1]==combinaisons[2,col],i], na.rm=TRUE)
      MeanRelAbcond2=c(MeanRelAbcond2,tt)
    }
    TABfbak <- TABf <- cbind(TABf, MeanRelAbcond1, MeanRelAbcond2)

    #Adjust table
    TABf <- TABf[!is.na(TABf$DESeqLFC),]
    TABf$Condition = rep(NA, nrow(TABf))
    TABf[TABf$DESeqLFC>0, "Condition"] = as.character(combinaisons[1,col])
    TABf[TABf$DESeqLFC<0, "Condition"] = as.character(combinaisons[2,col])
    TABf$Condition = factor(TABf$Condition,
                            levels=c(as.character(combinaisons[1,col]),as.character(combinaisons[2,col])) )

    TABfinal <- rbind(TABfinal,TABf)

    # Barplot
    ## differentialy abundant features on 2 methods
    ## Top abs(LogFolchange)
    ## feature with relative abondance > 0.1% (0.001)
    TABbar = TABf[TABf$DESeq ==1 | TABf$metagenomeSeq ==1 |  TABf$sumMethods >=2, ]
    TABbar = TABbar[TABbar$MeanRelAbcond1>=0.001 | TABbar$MeanRelAbcond2>=0.001, ]
    TABbar = tail(TABbar[order(abs(TABbar$DESeqLFC)),],50)

    if(nrow(TABbar)){
      TABbar$tax = ttax[as.character(TABbar$ListAllOtu),rank]

      # png(paste(output,'/topDiffbarplot_',column1,'_',paste(combinaisons[,col],collapse="_vs_"),'.png',sep=''), width=20, height=20, units="cm", res=200)
      pbarplot <- p <-ggplot(data=TABbar, aes(x=reorder(tax, -abs(DESeqLFC)), y=DESeqLFC, fill=Condition ) ) +
        geom_bar(stat="identity", alpha = 1) + ggtitle(paste(combinaisons[,col],collapse="_vs_")) + labs(x='Features') +
        coord_flip() + theme_bw() +
        scale_y_continuous(minor_breaks = seq(-1E4 , 1E4, 1), breaks = seq(-1E4, 1E4, 5)) +
        theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold"),
        strip.text.x = element_text(size = 18,face="bold"),
        title=element_text(size=16,face="bold"))
      # print(p)
      # dev.off()

      ggsave(paste(output,'/topDiffbarplot_',column1,'_',paste(combinaisons[,col],collapse="_vs_"),'.eps',sep='')
        , plot=pbarplot, height = 20, width = 20, units="cm", dpi = 500, device="eps")



      outF[[paste(combinaisons[,col],collapse="_vs_")]] = list(plot = p)

    }else{flog.info('No ASV to plot...')}

    # Krona ? Diversity of differentialy abundant features.

    ## Venn diag pour chaque comparaison.
    # flog.info('Plotting Venn diagrams...')
    # TF = Filter(length, TF)
    # venn.plot <- venn.diagram(TF, filename = NULL, col = "black",
    #                           fill = rainbow(length(TF)), alpha = 0.50,
    #                           cex = 2, cat.col = 1, , lty = "blank",
    #                           cat.cex = 2.5, cat.fontface = "bold",
    #                           margin = 0.07, main=paste(combinaisons[,col],collapse="_vs_"), main.cex=2.5,
    #                           fontfamily ="Arial",main.fontfamily="Arial",cat.fontfamily="Arial");
    # png(paste(output,'/venndiag_',column1,'_',paste(combinaisons[,col],collapse="_vs_"),'.png',sep=''), width=20, height=20, units="cm", res=200)
    # grid.draw(venn.plot)
    # dev.off()

    # grid.draw(venn.plot)
  }


  # print(TABfinal)
  if(length(TABfinal$ListAllOtu) > 0){
    flog.info('Building csv file...')

    TAX <- cbind(seqid = rownames(ttax[as.character(TABfinal$ListAllOtu),]),as.data.frame(ttax[as.character(TABfinal$ListAllOtu),]))

    sseq <- cbind(seqid = row.names(as.data.frame(sseq)), as.data.frame(sseq))
    SEQ <- sseq[as.character(TABfinal$ListAllOtu),]


    colnames(TABfinal)[1:2] = c("seqid", "Comparaison")
    TABfinal = cbind(TABfinal, TAX[,-1], sequence = SEQ[,-1])

    write.table(TABfinal, paste(output,"/aggregate_diff_",column1,'.csv', sep=""), row.names=FALSE, sep="\t", quote = FALSE)



    outF[["table"]] = TABfinal
    if(returnval){return(outF)}

    flog.info('Done.')

    # flog.info('Heatmap of significant %s.', rank)
    # signif_data <- phyloseq::prune_taxa(as.character(TABfinal[,1]), data)
    # if(rank != 'ASV'){
    #   signif_data <- tax_glom(signif_data, rank)
    #   taxa_names(signif_data) <- as.character(tax_table(signif_data)[,rank])
    # }
    # clr = function(x){log(x+1) - rowMeans(log(x+1))}
    # otable <- otu_table(signif_data)
    # otableCLR <- clr(otable)
    # data.norm <- signif_data; otu_table(data.norm) <- otableCLR
    #
    # hgroup=hclust(dist(otable) , method="ward.D2")
    # # plot(hgroup)
    #
    # #Heatmap
    # dataF <- data.norm
    # adf = psmelt(dataF)
    # adf[,"OTU"] = factor(adf[,"OTU"], levels=hgroup$labels)
    #
    # p = ggplot(adf, aes(x = Sample, y = OTU, fill = Abundance)) +
    #     geom_raster() + facet_grid(as.formula(paste("~",column1)), scales = "free", space = "free") +
    #     scale_fill_distiller(direction=-1, palette='Spectral') +
    #     theme(axis.text.x = element_text(angle=90, hjust=1)) + ylab(rank) +
    #     labs(fill='CLR Normalized\nabundance')
    # png(paste(output,'/heatmap_signif_',rank,'.png',sep=''), width=30, height=20, units="cm", res=200)
    # plot(p)
    # invisible(dev.off())
    # #Corrplot?
    # data <- signif_data
    # save(data, file=paste(output,"/signif_phyloseq.rdata",sep=""))
    #

  }else{flog.info("No significant results in the three methods...")}


}
