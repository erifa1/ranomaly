#' Generate Phyloseq Object
#'
#'
#'
#' @param dada_res Results of dada2_fun()
#' @param taxtable Results of assign_taxo_fun()
#' @param tree Results of generate_tree_fun()
#' @param metadata Path of metadata file (tab separated with header)
#' @param output Output directory
#' @param verbose Verbose level. (1: quiet, 3: verbal)
#' @param returnval Boolean to return values in console or not.
#'
#' @return Return a complete phyloseq object.
#' @import dada2
#' @import phyloseq
#' @import DECIPHER
#' @import ShortRead
#' @import Biostrings
#' @import futile.logger
#' @import digest
#'
#' @export


# Generate Phyloseq

generate_phyloseq_fun <- function(dada_res = dada_res, taxtable = tax.table, tree = tree, metadata = "",
                                  output = "./phyloseq/", verbose = 1, returnval = TRUE){

  flog.info("Loading sample metadata..")
  sampledata <- read.table(metadata, sep="\t",header=TRUE)
  rownames(sampledata) <- sampledata$sample.id
  print(rownames(sampledata))
  sample.metadata <- sample_data(sampledata)
  flog.info('Done.')

  if(length(setdiff(colnames(dada_res$otu.table), sampledata$sample.id)) > 0){
    flog.info(setdiff(colnames(dada_res$otu.table), sampledata$sample.id))
    stop("ERROR: number of samples in metadata differ from otu table.")
  }

  flog.info("Sequences..")
  sequences <- getSequences(dada_res$seqtab.nochim)
  names(sequences) <- sapply(sequences,digest,algo='md5')
  data <- phyloseq(dada_res$otu.table, tax.table, sample.metadata, phy_tree(tree), DNAStringSet(sequences))
  data_rel <- transform_sample_counts(data, function(x) x / sum(x) )

  if(!dir.exists(output)){
    flog.info('Creating output directory...')
    dir.create(output)
    flog.info('Done.')
  }
  flog.info('Saving R objects.')
  save(data, data_rel, file=paste(output,'/robjects.Rdata',sep=''))
  flog.info('Finish.')

  if(returnval){return(data)}

}
