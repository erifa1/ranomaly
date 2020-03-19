#' Generate Tree
#'
#' Generate Tree
#'
#' @param dada_res output from dada2_fun
#'
#' @return Return raw otu table in phyloseq object.
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

generate_phyloseq_fun <- function(otutable = dada_res, taxtable = tax.table, tree = tree, metadata = "",
                                  output = "./phyloseq/", verbose = 1, returnval = TRUE){

  flog.info("Loading sample metadata..")
  sampledata <- read.table(metadata, sep="\t",header=TRUE)
  rownames(sampledata) <- sampledata$sample.id
  print(rownames(sampledata))
  sample.metadata <- sample_data(sampledata)
  flog.info('Done.')

  if(length(setdiff(colnames(otutable$otu.table), sampledata$sample.id)) > 0){
    flog.info(setdiff(colnames(otutable$otu.table), sampledata$sample.id))
    stop("ERROR: number of samples in metadata differ from otu table.")
  }


  sequences <- getSequences(otutable$seqtab.nochim)
  names(sequences) <- sapply(sequences,digest,algo='md5')
  data <- phyloseq(otutable$otu.table, tax.table, sample.metadata, phy_tree(tree), DNAStringSet(sequences))
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
