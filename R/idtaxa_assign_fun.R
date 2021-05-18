#' IDTAXA assign
#'
#' Additionnal function for assign_taxo and assign_fasta function.
#'
#' @param db_file Path to IDTAXA formatted database file.
#' @param dna A ‘DNAStringSet’ of unaligned sequences.
#' @param asv_names sequences IDs in same order.
#' @param confidence Bootstrap threshold 0...100
#' @param ncpu Number of cpu to use
#'
#' @return return taxonomic assignment of given sequences.
#' @import futile.logger
#' @import DECIPHER

#' @export


idTaxa_assign = function(db_file, dna, asv_names, confidence, ncpu = NULLS){
  flog.info(paste('Using database ',db_file,sep=''))
  toto <- load(db_file)
  ids <- IdTaxa(dna, trainingSet, strand="both", processors=ncpu, verbose=TRUE)
  names(ids) <- asv_names
  flog.info("Confidence filtering...")
  IDCONF = as.numeric(confidence)
  ids2 = ids
  for (i in 1:length(ids2)){
    ID = ids2[[i]]
    ids2[[i]]$taxon = ID$taxon[ID$confidence>IDCONF]
    ids2[[i]]$confidence = ID$confidence[ID$confidence>IDCONF]
  }
  # deleting taxon field that starts with unclassified.
  for (seq_name in asv_names){
    for(i in 2:8){
      if(!is.na(ids2[[seq_name]]$taxon[i])){
        if(startsWith(ids2[[seq_name]]$taxon[i], "unclassified_")){
          # print(i)
          # print(length(ids2[[seq_name]]$taxon))
          # print(ids2[[seq_name]]$taxon[-c(i:length(ids2[[seq_name]]$taxon))])
          ids2[[seq_name]]$taxon <- ids2[[seq_name]]$taxon[-c(i:length(ids2[[seq_name]]$taxon))]
          ids2[[seq_name]]$confidence <- ids2[[seq_name]]$confidence[-c(i:length(ids2[[seq_name]]$confidence))]
        }
      }
    }
  }

  return(ids2)
}
