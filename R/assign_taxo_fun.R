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
#'
#'
#' @return Return a taxonomy table with multiple ancestor checking and incongruence checking when more than one databases are used.
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



  ## IDATAXA assignment
  dna <- DNAStringSet(getSequences(dada_res$seqtab.nochim))
  names(dna) <- sapply(colnames(dada_res$seqtab.nochim), digest::digest, algo="md5")

  if(!all( names(dna) == rownames(dada_res$otu.table) )){stop("Seq names and ASV table row names are different")}

  tt2 = idtaxa_assign_fasta_fun(fasta = dna, id_db = id_db,
        output = output, confidence = confidence, verbose = verbose, returnval = returnval)


}
