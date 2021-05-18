#' ASV taxonomic assignment with DECIPHER::IDTAXA (assign_taxo)
#'
#'
#'
#' @param dada_res Results of dada2_fun()
#' @param output Output directory
#' @param id_db Vector with list of absolute path to IDTAXA formatted reference database(s) (up to 2 databases).
#' @param verbose Verbose level. (1: quiet, 3: verbal)
#' @param confidence Bootstrap threshold 0...100
#' @param returnval Boolean to return values in console or not.
#' @param ncpu Number of cpus to use.
#'
#'
#' @return Return a taxonomy table with multiple ancestor checking and incongruence checking when 2 databases are used.
#'
#' @import futile.logger
#' @import dada2
#' @import phyloseq
#' @import DECIPHER
#' @import ShortRead
#' @import Biostrings
#' @export



assign_taxo_fun <- function(dada_res = dada_res,  output = "./idtaxa/", id_db = "/PathToDB/UNITE_idtaxa.Rdata", confidence = 50, verbose = 1, returnval = TRUE, ncpu=NULL){


  if(verbose == 3){
    invisible(flog.threshold(DEBUG))
  } else {
    invisible(flog.threshold(INFO))
  }


  if(!dir.exists(output)){
    flog.debug('Creating output directory...')
    dir.create(output)
    flog.debug('Done.')
  }

  ## IDATAXA assignment
  flog.info(paste('Taxonomy assignation with IDTAXA',sep=''))
  dna <- DNAStringSet(getSequences(dada_res$seqtab.nochim))


  db_list <- unlist(strsplit(id_db,","))
  taxid_list=vector("list", length(db_list)+1)
  for (i in 1:length(db_list)){
    db_file <- db_list[i]
    taxid_list[[i]] <- idTaxa_assign(db_file, dna, colnames(dada_res$seqtab.export), confidence, ncpu)
  }

  if(verbose == 3){
    save(taxid_list, file = "./annot_taxid_list_debug.rdata")
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
    write.table(final_tax_table, paste(output,"/allDB_tax_table.csv",sep=""), quote = FALSE, sep = "\t", row.names = FALSE)

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
  FeatNames = colnames(taxid); if(is.null(FeatNames)){FeatNames = names(taxid)}
  taxid <- t(as.data.frame(taxid))
  row.names(taxid) = FeatNames
  flog.info('Done.')

  if(verbose == 3){
    save(taxid, file = "./annot_taxid_debug.rdata")
  }


  # Filling taxonomy with last assigned rank.
  names(taxid) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  flog.info("Filling missing taxonomy ranks...")
  PREFIX = c("k__","p__","c__","o__","f__","g__","s__")

  if( length(grep("p__",taxid[,2])) == 0 ){
    taxid = fill_tax_fun(taxid, prefix = FALSE)
    taxid = as.data.frame( t(apply(taxid, 1, function(x){ paste(PREFIX, x, sep="")})), stringAsFactors = FALSE)
  }else{
    taxid = fill_tax_fun(taxid, prefix = TRUE)
  }


  flog.info('Done.')

  names(taxid) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  tax.table <- tax_table(as.matrix(taxid))


  flog.info("Check taxonomy consistency...")
  tax.table = check_tax_fun(tax.table, output = NULL)
  flog.info("Done.")

  #Output table 2
  Tab2 = apply( as.data.frame(tax.table, stringsAsFactors=FALSE) , 1 , paste , collapse = ";" )
  Tab2 <- cbind(Tab2, seq=colnames(dada_res$seqtab.nochim))
  write.table(Tab2, paste(output,"/final_tax_table.csv",sep=""), quote = FALSE, sep = "\t",
              col.names = FALSE, row.names = TRUE)

  flog.info('Saving R objects.')
  save(tax.table,  file=paste(output,'/robjects.Rdata',sep=''))
  flog.info('Finish.')

  if(returnval){return(tax.table)}



}
