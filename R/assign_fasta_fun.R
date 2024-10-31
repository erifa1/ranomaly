#' IDTAXA assign
#'
#' Same function as assign_taxo_fun but use FASTA file as input.
#'
#' @param fasta Path to fasta file, or DNAStringSet object.
#' @param output Output directory
#' @param id_db Vector with list of absolute path to IDTAXA formatted reference database(s).
#' @param verbose Verbose level. (1: quiet, 3: verbal)
#' @param confidence Bootstrap threshold 0...100
#' @param returnval Boolean to return values in console or not.
#' @param ncpu Number of cpu to use. 
#' @param prefix Vector of prefixes to use.
#' @param ranks_names Taxonomy ranks names.
#'
#'
#' @return Return a taxonomy table with multiple ancestor checking and incongruence checking when more than one databases are used.
#'
#' @import futile.logger
#' @import DECIPHER

#' @export

idtaxa_assign_fasta_fun <- function(fasta, id_db, output = "./assign_fasta/", confidence = 50, verbose = 1, 
  returnval = TRUE, ncpu=NULL, prefix = c("k__","p__","c__","o__","f__","g__","s__"), 
  ranks_names = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")){
  if(verbose == 3){
    flog.threshold(DEBUG)
  }
  if(!dir.exists(output)){
    flog.debug('Creating output directory...')
    dir.create(output, recursive = TRUE)
    flog.debug('Done.')
  }

  flog.info(paste('Taxonomy assignation with IDTAXA',sep=''))
  if(class(fasta) == "DNAStringSet"){
    dna = fasta
  } else {
    dna <- readDNAStringSet(fasta)
  }
  dna <- RemoveGaps(dna)

  db_list <- unlist(strsplit(id_db,","))
  if(any(!file.exists(c(db_list)))){stop("One or more reference DB unreachable. Check path(s).")}

  taxid_list=vector("list", length(db_list)+1)
  for (i in 1:length(db_list)){
    db_file <- db_list[i]
    taxid_list[[i]] <- idTaxa_assign(db_file, dna, names(dna), confidence, ncpu)
  }

  if(verbose == 3){
    save(taxid_list, file = "./annot_taxid_list_debug.rdata")
  }

  if(length(db_list)>1){
    # Multiple references
    flog.info('Merging all annotations...')

    Fannot = list()
    for (seq_name in names(dna)){
      flog.debug(seq_name)

      # create list
      L1=list(); L2=list()
      for (i in 1:length(db_list)){
        L1[[i]] =  taxid_list[[i]][[seq_name]]$taxon
        L2[[i]] =  taxid_list[[i]][[seq_name]]$confidence
      }

      depth1 = sapply(L1, length)
      if ( all(depth1 == 0) ){
        flog.debug("unassigned")
        Fannot[[seq_name]]$taxon = rep("unassigned", 8)
        Fannot[[seq_name]]$confidence = rep(0, 8)
      }
      else {
        if(length(unique(depth1)) == 1){
          flog.debug("best conf")
          #if same assignation depth, test on confidence
          last_rank <- length(taxid_list[[1]][[seq_name]]$taxon)
          conf1 = sapply(L2, function(x){x[last_rank]})
          Fannot[[seq_name]] <- taxid_list[[which(conf1 == max(conf1))[1]]][[seq_name]] # ..[1] in case of same depth / same conf.
        }
        else {
          if( any( table(depth1) > 1 ) ){
            # if some equal depth
            flog.debug("some equal depth, best conf")
            last_rank <- max(sapply(L1, length))
            maxDepth_list = L2[which(depth1 == max(depth1))]
            conf1 = sapply(L2, function(x){x[last_rank]}); conf1[is.na(conf1)] = 0
            Fannot[[seq_name]] <- taxid_list[[which(conf1 == max(conf1))]][[seq_name]]
          }
          else{
            #if all different assignation depth
            flog.debug("best depth")
            Fannot[[seq_name]] <- taxid_list[[which(depth1 == max(depth1))]][[seq_name]]
          }
        }
      }
      flog.debug(paste('length Fannot:',length(Fannot)))
    }
    flog.info('Done.')

    flog.info('Output table...')
    #Process output table with different assignations.
    all_annot_tab = data.frame(row.names=names(taxid_list[[1]]))
    for (i in 1:length(db_list)){

      ALLassign <- sapply(taxid_list[[i]],function (id) {
        paste(id$taxon," (",round(id$confidence, digits=1),"%)",sep="",collapse="; ")
      })
      all_annot_tab[,i] <- ALLassign

    }

    FINALassign <- sapply(Fannot,function (id) {
      paste(id$taxon," (",round(id$confidence, digits=1),"%)",sep="",collapse="; ")
    })

    all_annot_tab$FINAL = FINALassign

    seq_tab <- data.frame(sequences=as.character(dna))
    row.names(seq_tab) <- names(dna)
    final_tax_table <- merge(all_annot_tab, seq_tab, by = "row.names")

    colnames(final_tax_table)=c("ASV",db_list, "FINAL", "sequences")
    write.table(final_tax_table, paste(output,"/allDB_tax_table.csv",sep=""), quote = FALSE, sep = "\t", row.names = FALSE)

    # Keep merged assignment (Fannot)
    taxid <- sapply(Fannot, function(x) {
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


  }else{
    # One DB assignment
    taxid <- sapply(taxid_list[[1]], function(x) {
      taxa <- rep(NA,length(ranks_names))
      assign <- x$taxon
      # print(assign)
      if(length(assign[-1])==0){
        taxa <- rep("unassigned",length(ranks_names))
      } else {
        taxa[1:length(assign[-1])] <- assign[-1]
      }

      return(taxa)
    })
  }

  if(max(sapply(taxid, length)) > 7){
    # taxid is list when non homogeneous length in taxonomy
    LL=sapply(taxid, length)
    Tnames=names(LL[LL>7])
    for(nn in Tnames){
      taxid[[nn]] = taxid[[nn]][1:7]
    }
    taxnames <- names(taxid)
    taxid <- t(as.data.frame(taxid))
    row.names(taxid) <- taxnames
  }else{
    # taxid is a matrix
    taxnames <- colnames(taxid)
    taxid <- t(as.data.frame(taxid))
    row.names(taxid) <- taxnames
  }

  flog.info('Done.')

  if(verbose == 3){
    save(taxid, file = "./annot_taxid_debug.rdata")
  }

  # Filling taxonomy with last assigned rank.
  colnames(taxid) <- ranks_names
  flog.info("Filling missing taxonomy ranks...")
  # PREFIX = prefix
  # handling possible prefix
  # noprefix_taxid = taxid[grep("k__",taxid[,1], invert = TRUE),]
  # prefix_taxid = taxid[grep("k__",taxid[,1]),]

  # if(nrow(noprefix_taxid) != 0){
  #   noprefix_taxid = fill_tax_fun(noprefix_taxid, prefix = FALSE)
  #   noprefix_taxid = as.data.frame( t(apply(noprefix_taxid, 1, function(x){ paste(PREFIX, x, sep="")})), stringAsFactors = FALSE)
  #   colnames(noprefix_taxid) = colnames(taxid)
  # }

  # if(nrow(prefix_taxid) != 0){
  #   prefix_taxid = fill_tax_fun(prefix_taxid, prefix = TRUE)
  # }
  
  final_taxid = fill_tax_fun(taxid, prefix = prefix, ranks_names = ranks_names)

  tax.table = final_taxid[names(dna),]

  flog.info('Done.')

  flog.info("Check taxonomy consistency...")
  tax.tablecheck = check_tax_fun(tax.table, output = NULL, rank = length(ranks_names) - 1, ranks_names = ranks_names, verbose = 3)
  flog.info("Done.")

  #Output table 2
  tax.tablecheck[,1] <- stringr::str_replace(tax.tablecheck[,1], "d__", "k__") # replace d__ prefix by k__ for kingdom

  Tab2 = apply( as.data.frame(tax.tablecheck, stringsAsFactors=FALSE) , 1 , paste , collapse = ";" )
  Tab2 <- cbind(Tab2, seq=as.character(dna))
  write.table(Tab2, paste(output,"/final_tax_table.csv",sep=""), quote = FALSE, sep = "\t",
              col.names = FALSE, row.names = TRUE)

  flog.info('Saving R objects.')
  save(tax.tablecheck,  file=paste(output,'/robjects.Rdata',sep=''))
  flog.info('Finish.')

  if(returnval){return(tax.tablecheck)}


}
