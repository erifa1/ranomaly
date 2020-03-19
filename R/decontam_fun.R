#' Decontam Function
#'
#' Decontam Function
#'
#' @param dada_res output from dada2_fun
#'
#' @return Return raw otu table in phyloseq object.
#' @import decontam
#' @import phyloseq
#' @import ggplot2
#' @import futile.logger
#' @import psadd
#' @import VennDiagram
#' @importFrom grid grid.draw
#'
#' @export


# Decontam Function

decontam_fun <- function(data = data, output = "./decontam_out/", number = 4000, prev = 2, freq = 0.00005,
                         column = "type", ctrl_identifier = "control", spl_identifier = "sample", batch = NULL, plot = FALSE,
                         method = "prevalence", threshold = 0.1, concentration = NULL, verbose = 1, domain = TRUE, unassigned = FALSE,
                         skip = FALSE, manual_cont_rank = "Genus", manual_cont = NULL, returnval=TRUE){

  invisible(flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger"))


  if(!dir.exists(output)){
    flog.info('Creating output directory...')
    dir.create(output)
    flog.info('Done.')
  }

  listInit <- rownames(otu_table(data))
  data_no_filtering <- data

  if(plot){
    flog.info('Plotting...')
    df <- as.data.frame(sample_data(data))
    df$LibrarySize <- sample_sums(data)
    df <- df[order(df$LibrarySize),]
    df$Index <- seq(nrow(df))
    p <- ggplot(data=df, aes_string(x="Index", y="LibrarySize", color=column)) + geom_point()
    ggsave(paste(output,'/lib_size.png',sep=''), plot=p)
    flog.info('Done.')

    flog.info('Generating Krona...')
    plot_krona(data, paste(output,'/krona_no_filtering',sep=''),'sample.id')
    flog.info('Done.')
  }

  # CHECKING CONTROL SAMPLES
  if(skip == FALSE){
    df <- as.data.frame(sample_data(data))
    print(column)
    samplesCol = df[, column];
    if(!is.null(batch)){
      df2 = df[,c(batch,column)]
      df3 = df2[df2[,column]==ctrl_identifier,]
      flog.info('Per batch control samples counts:')
      print(table(df3[,batch]))
      batchCtrl <- median(table(df3[,batch]))
    } else {batchCtrl <- 1000}
    #DECONTAM STEP
    if(skip == FALSE & nrow(samplesCol[samplesCol==ctrl_identifier]) >= 3 & batchCtrl >= 3){
      flog.info('Decontam step...')
      if(method == 'frequency'){
        flog.info('Method frequency...')
        flog.debug(batch)
        if(is.null(batch)){
          fun <- paste('contamdf <- isContaminant(data, method="frequency", conc=as.numeric(sample_data(data)$',concentration,'))')
          eval(parse(text=fun))
        } else{
          fun <- paste('contamdf <- isContaminant(data, method="frequency", conc=as.numeric(sample_data(data)$',concentration,'), batch=',batch,')')
          eval(parse(text=fun))
        }
        head(contamdf)
        table(contamdf$contaminant)

        head(which(contamdf$contaminant))
        flog.info('Done.')
        if(plot){
          flog.info('Plotting...')
          set.seed(100)
          p <- plot_frequency(data, taxa_names(data)[sample(which(contamdf$contaminant),3)], conc="quant_reading") + xlab("DNA Concentration (PicoGreen fluorescent intensity)")
          ggs
          save(paste(output,'/freq_conta_exemple.png',sep=''), plot=p)
          flog.info('Done.')
        }
      } else if(method == 'prevalence'){
        # sample_data(data)$is.neg <- sample_data(data)$status == 'control'
        flog.info('Method prevalence...')
        tmp <- paste('sample_data(data)$is.neg <- sample_data(data)$',column,' == \'',ctrl_identifier,'\'',sep='')
        eval(parse(text=tmp))
        flog.debug(batch)
        if(is.null(batch)){
          contamdf <- isContaminant(data, method='prevalence', neg='is.neg')
        }else{
          contamdf <- isContaminant(data, method='prevalence', neg='is.neg', batch=batch)
        }
        flog.info('Done.')
        table(contamdf$contaminant)
        if(plot){
          flog.info('Plotting...')
          data.pa <- transform_sample_counts(data, function(abund) 1*(abund>0))

          # data.pa.neg <- prune_samples(sample_data(data.pa)$status == 'control', data.pa)
          tmp <- paste('data.pa.neg <- prune_samples(sample_data(data.pa)$', column, ' == \'',ctrl_identifier, '\', data.pa)', sep="")
          eval(parse(text=tmp))
          # data.pa.pos <- prune_samples(sample_data(data.pa)$status == 'sample', data.pa)
          tmp <- paste('data.pa.pos <- prune_samples(sample_data(data.pa)$', column, ' == \'',spl_identifier, '\', data.pa)', sep="")
          eval(parse(text=tmp))
          # Make data.frame of prevalence in positive and negative samples
          df.pa <- data.frame(pa.pos=taxa_sums(data.pa.pos), pa.neg=taxa_sums(data.pa.neg), contaminant=contamdf$contaminant)
          p <- ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() + xlab("Prevalence (Controls)") + ylab("Prevalence (Samples)")
          ggsave(paste(output,'/prevalence_conta_exemple.png',sep=''), plot=p)
          flog.info('Done.')
        }
      } else if(method == 'combined' | method == 'either' | method == 'both'){
        flog.info('Method %s',method)

        # sample_data(data)$is.neg <- sample_data(data)$status == 'control'
        tmp <- paste('sample_data(data)$is.neg <- sample_data(data)$',column,' == \'',ctrl_identifier,'\'',sep='')
        eval(parse(text=tmp))


        if(is.null(batch)){
          fun <- paste('contamdf <- isContaminant(data, method=\"',method,'\", neg="is.neg", conc=as.numeric(sample_data(data)$',concentration,'), threshold=',threshold,')',sep='')
          eval(parse(text=fun))
        } else{
          fun <- paste('contamdf <- isContaminant(data, method=\"',method,'\", neg="is.neg", conc=as.numeric(sample_data(data)$',concentration,'),threshold=',threshold,',batch=sample_data(data)$',batch,')',sep='')
          eval(parse(text=fun))
        }

        table(contamdf$contaminant)
        flog.info('Done.')
        if(plot){
          flog.info('Plotting...')
          data.pa <- transform_sample_counts(data, function(abund) 1*(abund>0))
          # data.pa.neg <- prune_samples(sample_data(data.pa)$status == 'control', data.pa)
          tmp <- paste('data.pa.neg <- prune_samples(sample_data(data.pa)$', column, ' == \'',ctrl_identifier, '\', data.pa)', sep="")
          eval(parse(text=tmp))
          # data.pa.pos <- prune_samples(sample_data(data.pa)$status == 'sample', data.pa)
          tmp <- paste('data.pa.pos <- prune_samples(sample_data(data.pa)$', column, ' == \'',spl_identifier, '\', data.pa)', sep="")
          eval(parse(text=tmp))
          # Make data.frame of prevalence in positive and negative samples
          df.pa <- data.frame(pa.pos=taxa_sums(data.pa.pos), pa.neg=taxa_sums(data.pa.neg), contaminant=contamdf$contaminant)
          p <- ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() + xlab("Prevalence (Controls)") + ylab("Prevalence (Samples)")
          ggsave(paste(output,'/prevalence_conta_exemple.png',sep=''), plot=p)
          p <- plot_frequency(data, taxa_names(data)[sample(which(contamdf$contaminant),3)], conc="quant_reading") + xlab("DNA Concentration (PicoGreen fluorescent intensity)")
          ggsave(paste(output,'/freq_conta_exemple.png',sep=''), plot=p)
          flog.info('Done.')
        }
      } else{
        print_help(opt_parser)
        stop('Wrong method provided.', call.=FALSE)
      }

      flog.debug('Table before contaminant filtering :')
      flog.debug(show(data))
      taxToKeep0 <- row.names(contamdf[contamdf$contaminant==FALSE,])
      taxToDump0 <- row.names(contamdf[contamdf$contaminant==TRUE,])

    } else {
      flog.info('Decontam step skipped or too few control samples (less than 3)...')
      taxToDump0 = NULL
    }
  }else{taxToDump0 = NULL} #-k skip all decontam


  #CLASSIC FILTERING FREQ / PREV / NREADS
  flog.info('Frequence filtering...')
  sumTot <- sum(otu_table(data))
  freqGlobale <- apply(otu_table(data), 1, function(x){sum(x)/sumTot})

  taxToKeep1=names(freqGlobale[freqGlobale>freq])
  taxToDump1=names(freqGlobale[freqGlobale<freq])
  flog.info(paste('Min frequency is ', signif(min(freqGlobale),4),'. ',length(taxToDump1),' ASVs are under ',freq,sep=''))

  flog.info('Prevalence filtering...')
  prevdf <- apply(X = otu_table(data), MARGIN = ifelse(taxa_are_rows(data), yes = 1, no = 2), FUN = function(x){sum(x > 0)})
  prevdf <- data.frame(Prevalence = prevdf, TotalAbundance = taxa_sums(data))
  taxToKeep2 <- rownames(prevdf)[(prevdf$Prevalence >= prev)]
  taxToDump2 <- rownames(prevdf)[(prevdf$Prevalence < prev)]

  data_rel <- transform_sample_counts(data, function(x) x / sum(x) )
  flog.info('Done.')

  #UNASSIGNED TAXA
  if(unassigned==TRUE){
    flog.info('Filtering unassigned domain...')
    if(domain==TRUE){
      k <- "k__Bacteria"
      p <- "p__Bacteria_Phylum"
    }else{
      k <- "k__Fungi"
      p <- "p__Fungi_Phylum"
    }
    flog.debug('Domain : %s ; Phylum : %s',k,p)

    taxdf <- as.data.frame(data@tax_table@.Data)
    taxToKeep3 = row.names(taxdf[taxdf$Domain==k,])
    taxToDump3 = row.names(taxdf[taxdf$Domain!=k,])
    flog.info('Done.')

    flog.info('Filtering unassigned phylum...')
    taxToKeep4 = row.names(taxdf[taxdf$Phylum!=p,])
    taxToDump4 = row.names(taxdf[taxdf$Phylum==p,])
    flog.info('Done.')

    taxToKeep4 = unique(c(taxToKeep3,taxToKeep4))
    taxToDump4 = unique(c(taxToDump3,taxToDump4))

  }else{
    taxToDump4 = NULL
  }


  flog.info(paste("BEFORE FILTERING: ",nsamples(data), "samples and", ntaxa(data),"ASVs in otu table") )
  uniqTaxToDump <- unique(c(taxToDump0,taxToDump1,taxToDump2,taxToDump4))
  allASV <- taxa_names(data)
  uniqTaxToKeep <- setdiff(allASV,uniqTaxToDump)

  dataKeep <- prune_taxa(uniqTaxToKeep,data)  #filtering

  TF = list(decontam=taxToDump0,freq=taxToDump1,prev=taxToDump2,unassigned=taxToDump4)
  TF = Filter(length, TF) #ommit empty field of list TF

  flog.info('Plotting Venn diagrams...')
  venn.plot <- venn.diagram(TF, filename = NULL, col = "black",
                            fill = rainbow(length(TF)), alpha = 0.50,
                            cex = 1.5, cat.col = 1, lty = "blank",
                            cat.cex = 1.8, cat.fontface = "bold",
                            margin = 0.1, main=paste("filtered ASVs"), main.cex=2.5,
                            fontfamily ="Arial",main.fontfamily="Arial",cat.fontfamily="Arial") #cat.dist = 0.09,
  png(paste(output,'/venndiag_filtering.png',sep=''), width=20, height=20, units="cm", res=200)
  grid.draw(venn.plot)
  invisible(dev.off())

  flog.info('Generate Exclu_out table...')
  # Tests in which method each ASV is filtered.
  TABf = otu_table(prune_taxa(uniqTaxToDump,data))
  for (j in 1:length(TF)){
    TABtest = TF[[j]]
    TABtest_filt=rep(NA, length(uniqTaxToDump))
    for (i in 1:length(uniqTaxToDump)) {
      featureI = uniqTaxToDump[i]
      res=grep( paste('^',featureI,'$', sep="" ) , TABtest)
      if(length(res)>0){TABtest_filt[i]=length(res)
      }
    }
    TABf=cbind.data.frame( TABtest_filt, TABf )
    names(TABf)[1] = names(TF)[j]
  }
  TABff <- cbind(as.matrix(TABf), as.matrix(data_no_filtering@tax_table[row.names(TABf),]))
  write.table(TABff, file = paste(output,'/Exclu_out.csv',sep=''), sep = "\t", col.names=NA)


  data <- dataKeep
  #NUMBER OF READS in samples
  flog.info(paste('Filtering samples with less than ',number,' reads...',sep=''))
  data <- prune_samples(sample_sums(data) > number, data)
  flog.info('Done.')
  flog.info(paste("AFTER FILTERING: ",nsamples(data), "samples and", ntaxa(data),"ASVs in otu table") )


  save.image("debug.rdata")
  ##GENUS to remove manually
  if(!is.null(manual_cont)){
    cont_list <- unlist(strsplit(manual_cont,","))
    flog.info(paste('Removing ',cont_list, sep=''))
    ttable = data@tax_table@.Data
    taxToKeep5=NULL
    for(i in cont_list){
      if(manual_cont_rank == "ASV"){
        taxToKeep5 = c(taxToKeep5,row.names(ttable)[row.names(ttable) !=i ]  )
      }else{
        taxToKeep5 = c(taxToKeep5,row.names(ttable)[ttable[,manual_cont_rank]!=i])
      }
    }
    data <- prune_taxa(taxToKeep5, data)
  }

  ##Remove Control samples for next analysis
  fun <- paste("data <- subset_samples(data, ",column," %in% '",spl_identifier,"')",sep="")
  eval(parse(text=fun))


  flog.info('Writing raw tables.')
  write.table(cbind(otu_table(data),"Consensus Lineage" = apply(tax_table(data), 1, paste, collapse = ";"), "sequences"=as.data.frame(refseq(data)) ), paste(output,"/raw_otu-table.csv",sep=''), sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

  flog.info('Writing relative tables.')
  write.table(cbind(otu_table(data_rel),"Consensus Lineage" = apply(tax_table(data_rel), 1, paste, collapse = ";"), "sequences"=as.data.frame(refseq(data_rel))),paste(output,"/relative_otu-table.csv",sep=''), sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

  flog.info('Saving R objects.')
  save(data, data_rel, file=paste(output,'/robjects.Rdata',sep=''))

  flog.info('Generating Krona.')
  plot_krona(data, paste(output,'/krona_filtering',sep=""),'sample.id')
  flog.info('Done.')
  flog.info('Finish.')

  if(returnval){return(data)}

}
