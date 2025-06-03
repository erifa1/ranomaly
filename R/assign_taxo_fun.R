#' ASV taxonomic assignment with DECIPHER::IDTAXA (assign_taxo)
#'
#'
#'
#' @param dada_res Results of dada2_fun()
#' @param output Output directory
#' @param id_db Vector with list of absolute path to IDTAXA formatted reference database(s).
#' @param verbose Verbose level. (1: quiet, 3: verbal)
#' @param confidence Bootstrap threshold 0...100
#' @param returnval Boolean to return values in console or not.
#' @param ncpu Number of cpus to use.
#' @param prefix Vector of prefixes to use.
#' @param ranks_names Taxonomy ranks names.
#'
#'
#' @return Return a taxonomy table with multiple ancestor checking and incongruence checking when more than one databases are used.
#'
#' @import futile.logger
#' @import dada2
#' @import phyloseq
#' @import DECIPHER
#' @import ShortRead
#' @export



assign_taxo_fun <- function(dada_res = dada_res,  output = "./idtaxa/", id_db = "/PathToDB/UNITE_idtaxa.Rdata", 
  confidence = 50, verbose = 1, returnval = TRUE, ncpu=NULL, prefix = c("k__","p__","c__","o__","f__","g__","s__"), 
  ranks_names = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")){


  if(verbose == 3){
    invisible(flog.threshold(DEBUG))
  } else {
    invisible(flog.threshold(INFO))
  }



  ## IDATAXA assignment
  dna <- Biostrings::DNAStringSet(getSequences(dada_res$seqtab.nochim))
  names(dna) <- sapply(colnames(dada_res$seqtab.nochim), digest::digest, algo="md5")

  if(!all( names(dna) == rownames(dada_res$otu.table) )){stop("Seq names and ASV table row names are different")}

  tt2 = idtaxa_assign_fasta_fun(fasta = dna, id_db = id_db,
        output = output, confidence = confidence, verbose = verbose, returnval = returnval, 
        prefix = prefix, ranks_names = ranks_names)


}


#' add_blast_fun
#'
#' Add taxonomy information from taxonkit lineage table to phyloseq object.
#'
#' @param data Phyloseq object
#' @param lineage_table Lineage table from taxonkit software.
#' @param domain Domain to use (Bacteria or Fungi)
#' @param output Output directory
#' @param prefix Prefix for taxonomy to use.
#'
#' @return Phyloseq object with taxonomy table updated.
#'
#' @export


add_blast_fun <- function(data = data, lineage_table = NULL, domain = "Bacteria", output = "./add_blast/", prefix = c("k__","p__","c__","o__","f__","g__","s__")){

  ttable1 = as.matrix(tax_table(data))
  txkit2 <- txkit <- rio::import(lineage_table)

  if(domain == "Bacteria"){
    idx_ok = c(which(grepl("Bacteria", txkit[,4])), which(grepl("Archaea", txkit[,4])))
  }else if(domain == "Fungi"){
    idx_ok = c(which(grepl("Fungi", txkit[,4])))
    txkit[idx_ok,5] = gsub("Eukaryota;", "Fungi;",txkit[idx_ok,5])
  }else(
    stop("Domain not recognized")
  )
  txkit2 = txkit[idx_ok, ]

  txkit3 = txkit2[!duplicated(txkit2[,1]),]
  row.names(txkit3) = txkit3[,1]

  # browser()
  # Resépare la taxonomie en colonnes.
  list0 <- strsplit(txkit3[,5], ";")
  max_len <- max(sapply(list0, length))
  list1 <- lapply(list0, function(x) {
  length(x) <- max_len  # complète automatiquement avec NA
  return(x)
  })
  df0 <- ttable_txkit <- t(as.data.frame(list1))
  row.names(ttable_txkit) = txkit3[,1]

  # ttable_txkit = t(as.data.frame(strsplit(txkit3[,5], ";"), stringsAsFactors = FALSE))
  # row.names(ttable_txkit) = txkit3[,1]
  ttable_txkit[which(ttable_txkit == "")] = NA

  ttable_txkit2 = na.omit(ttable_txkit)

  #Remplace les taxonomies si plus complètes avec blast.
  for(i in row.names(ttable_txkit2)){
    if(length(grep("_[Ss]pecies", ttable1[i,7])) != 0 | length(grep("_unassigned", ttable1[i,7])) != 0){
      print(i)
      print(ttable_txkit2[i,])

      PREFIX = prefix
      Ftax = paste(PREFIX, ttable_txkit2[i,], sep="")
      print(rbind(ttable1[i,],Ftax))

      ttable1[i,] = sub(" ", "_", Ftax)
    }

  }

  filltable = fill_tax_fun(ttable1, prefix = TRUE)
  check1 = check_tax_fun(filltable, output = NULL)


  tax_table(data0) = as.matrix(check1)

  data <- data0

  if(!is.null(output)){
    dir.create(output, showWarnings = FALSE, recursive = TRUE)
    data_unassigned <- prune_taxa(as.vector(tax_table(data)[,1] == "unassigned"), data)
    Biostrings::writeXStringSet(refseq(data_unassigned), filepath = glue::glue("{output}/unassigned_asv.fasta"), format = "fasta", width = 10000 )

    if(all(taxa_names(data) == row.names(check1))){
      check1$seq = as.character(refseq(data))
      write.csv(check1, file = glue::glue("{output}/check_final_tax.csv"))
    }else(
      stop("ASV names are different")
    )
  }


  return(data)
}
