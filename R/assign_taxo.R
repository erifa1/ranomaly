#' Assign Taxo
#'
#' Processing DADA2 algorithm on raw sequences, return raw otu table in phyloseq object.
#'
#' @param amplicon Choose amplipcon "16S" or "ITS"
#'
#' @return Return raw otu table in phyloseq object.
#' @import futile.logger
#' @import dada2
#' @import phyloseq
#' @import DECIPHER
#' @import ShortRead
#' @import Biostrings
#' @export



assign_taxo <- function(dada_res = dada_res,  out = "./idtaxa/", id_db = "/home/db/unite/UNITE_idtaxa201706.Rdata", verbose = 1, confidence = 50){


  if(verbose == 3){
    invisible(flog.threshold(DEBUG))
  } else {
    invisible(flog.threshold(INFO))
  }


  if(!dir.exists(out)){
    flog.debug('Creating output directory...')
    dir.create(out)
    flog.debug('Done.')
  }

  ## IDATAXA assignment
  flog.info(paste('Taxonomy assignation with IDTAXA',sep=''))
  dna <- DNAStringSet(getSequences(dada_res$seqtab.nochim))


  db_list <- unlist(strsplit(id_db,","))
  taxid_list=vector("list", length(db_list)+1)
  for (i in 1:length(db_list)){
    db_file <- db_list[i]
    taxid_list[[i]] <- idTaxa_assign(db_file, dna, colnames(dada_res$seqtab.export), confidence)
  }

  if(verbose == 3){
    save.image("./annot_debug.rdata")
  }

  if(length(db_list)>1){
    flog.info('Merging both annotation...')
    for (seq_name in colnames(dada_res$seqtab.export)){
      #for each seq
      if(length(taxid_list[[1]][[seq_name]]$taxon) == 0 & length(taxid_list[[2]][[seq_name]]$taxon) == 0){
        taxid_list[[3]][[seq_name]]$taxon = rep("unassigned", 8)
        taxid_list[[3]][[seq_name]]$confidence = rep(0, 8)
      } else {
        if(length(taxid_list[[1]][[seq_name]]$taxon) == length(taxid_list[[2]][[seq_name]]$taxon)){
          #if same assignation depth, test on confidence
          last_rank <- length(taxid_list[[1]][[seq_name]]$taxon)
          if(taxid_list[[1]][[seq_name]]$confidence[last_rank] > taxid_list[[2]][[seq_name]]$confidence[last_rank]){
            taxid_list[[3]][[seq_name]] <- taxid_list[[1]][[seq_name]] #DB1
          } else{
            taxid_list[[3]][[seq_name]] <- taxid_list[[2]][[seq_name]] #DB2
          }
        } else {
          #if different assignation depth
          if(length(taxid_list[[1]][[seq_name]]$taxon) > length(taxid_list[[2]][[seq_name]]$taxon)){
            taxid_list[[3]][[seq_name]] <- taxid_list[[1]][[seq_name]] #DB1
          } else{
            taxid_list[[3]][[seq_name]] <- taxid_list[[2]][[seq_name]] #DB2
          }
        }

      }

    }
    flog.info('Done.')

    flog.info('Output table...')
    #Process output table with different assignations.
    all_annot_tab = data.frame(row.names=names(taxid_list[[1]]))
    for (i in 1:3){
      if(i<3){print(db_list[i])}else{print("Final")}
      ALLassign <- sapply(taxid_list[[i]],function (id) {
        paste(id$taxon," (",round(id$confidence, digits=1),"%)",sep="",collapse="; ")
      })
      all_annot_tab[,i] <- ALLassign
      if(verbose == 3){print(head(all_annot_tab))}
    }

    seq_tab <- data.frame(sequences=as.character(dna))
    row.names(seq_tab) <- colnames(dada_res$seqtab.export)
    final_tax_table <- merge(all_annot_tab, seq_tab, by = "row.names")


    colnames(final_tax_table)=c("ASV",db_list[1], db_list[2], "FINAL", "sequences")
    write.table(final_tax_table, paste(out,"/allDB_tax_table.csv",sep=""), quote = FALSE, sep = "\t", row.names = FALSE)

    # Keep merged assignment (taxid_list[[3]])
    taxid <- sapply(taxid_list[[3]], function(x) {
      taxa <- rep(NA,7)
      assign <- x$taxon
      # print(assign)
      if(length(assign[-1])==0){
        taxa <- rep("unassigned",7)
      } else {
        taxa[1:length(assign[-1])] <- assign[-1]
      }

      return(taxa)
    })

    # taxid <- as.data.frame(taxid)

  }else{
    # One DB assignment
    taxid <- sapply(taxid_list[[1]], function(x) {
      taxa <- rep(NA,7)
      assign <- x$taxon
      # print(assign)
      if(length(assign[-1])==0){
        taxa <- rep("unassigned",7)
      } else {
        taxa[1:length(assign[-1])] <- assign[-1]
      }

      return(taxa)
    })
  }

  if(max(sapply(taxid, length)) > 7){
    LL=sapply(taxid, length)
    Tnames=names(LL[LL>7])
    for(nn in Tnames){
      taxid[[nn]] = taxid[[nn]][1:7]
    }
  }
  taxid <- as.data.frame(taxid)
  flog.info('Done.')



  # Filling taxonomy with last assigned rank.
  row.names(taxid) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

  fill_tax_table = function(x){
    RANKS = c("domain","phylum","class","order","family","genus","species")
    PREFIX = c("k__","p__","c__","o__","f__","g__","s__")
    TAX = na.omit(as.character(x))
    lastRank = length(TAX)
    if(length(TAX)==0){
      TAX=c("unassigned")
      lastRank = 1
    }
    for( i in 1:7){
      if(i <= lastRank){
        if(grepl(PREFIX[i], TAX[i], ignore.case = FALSE, perl = FALSE, fixed = TRUE)){
          TAX[i] <- sub(PREFIX[i],'',TAX[i])
        }
        x[i] = paste(PREFIX[i],TAX[i], sep="")
      } else {
        if(grepl(PREFIX[lastRank], TAX[lastRank], ignore.case = FALSE, perl = FALSE, fixed = TRUE)){
          TAX[lastRank] <- sub(PREFIX[lastRank],'',TAX[lastRank])
        }
        x[i] = paste(PREFIX[i],TAX[lastRank],'_', RANKS[i], sep="")
      }
    }
    return(x)
  }

  flog.info("Filling missing taxonomy ranks...")
  taxid <- t(apply(taxid, 2, fill_tax_table))
  flog.info('Done.')

  tax.table <- tax_table(as.matrix(taxid))

  flog.info("Check taxonomy consistency...")
  # Check for multiple ancestors at each rank, choose first occurence for each problematic taxon
  for(rank in 6:2){
    if(verbose==3){flog.info(paste(colnames(tax.table)[rank],".",sep=""))}

    stockLres = 100; nloop=1
    while(max(stockLres)>1){
      if(verbose==3){flog.info(paste("Loop",nloop,".",sep=""))}
      stockLres = NULL
      allTax = apply(tax.table[,1:rank] ,1, function(x){paste(x, collapse = ";")})
      uniqTax <- unique(allTax)

      #Test if unique Genus has unique taxonomy.
      for (i in unique(tax.table[,rank])){
        res=grep( paste(';',i,'$', sep="" ) , uniqTax)
        stockLres = c(stockLres,length(res))
        #If multiple tax for one taxon...
        if(length(res)>1){
          cat("\n\n");print(i);print(uniqTax[res]);
          tax2 <- apply(tax.table[tax.table[,rank]==i,1:rank], 1, function(x){paste(x, collapse = ";")})
          uniqTax2 <- table(tax2)
          ftax <- names(uniqTax2[order(uniqTax2,decreasing=TRUE)])[1]
          ftax <- unlist(strsplit(ftax,";"))
          print(ftax)
          #Change taxonomy with final ftax. the most common in tax.table
          for(j in row.names(tax.table[tax.table[,rank]==i,])){
            tax.table[j,] = c(ftax, tax.table[j,(rank+1):ncol(tax.table)])
          }
        }
      }
      nloop = nloop + 1
    }
  }
  flog.info("Done.")

  #Output table 2
  Tab2 = apply( as.data.frame(taxid, stringsAsFactors=FALSE) , 1 , paste , collapse = ";" )
  Tab2 <- cbind(Tab2, seq=colnames(dada_res$seqtab.nochim))
  write.table(Tab2, paste(out,"/final_tax_table.csv",sep=""), quote = FALSE, sep = "\t",
              col.names = FALSE, row.names = TRUE)

  flog.info('Saving R objects.')
  save(tax.table,  file=paste(out,'/robjects.Rdata',sep=''))
  flog.info('Finish.')

  return(tax.table)



}
