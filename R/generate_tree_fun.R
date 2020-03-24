#' Generate Phylogenetic Tree
#'
#'
#'
#' @param dada_res Results of dada2_fun()
#' @param output Output directory
#' @param returnval Boolean to return values in console or not.
#'
#' @return Return a formatted tree object ready to use in phyloseq.
#'
#' @import dada2
#' @import phyloseq
#' @import DECIPHER
#' @import ShortRead
#' @import Biostrings
#' @import futile.logger
#' @import digest
#' @import phangorn
#' @export


# DADA2 function

generate_tree_fun <- function(dada_res = dada_res, output = "./tree", returnval=TRUE){


  flog.info('Generating tree...')
  sequences <- getSequences(dada_res$seqtab.nochim)
  names(sequences) <- sapply(sequences,digest,algo='md5')
  flog.info('Aligning sequences...')
  alignment <- AlignSeqs(DNAStringSet(sequences),anchor=NA,processors=NULL)
  flog.info('Creating distance matrices...')
  phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
  dm <- dist.ml(phang.align)
  flog.info('Neigbour joining...')
  treeNJ <- NJ(dm) # Note, tip order != sequence order
  fit = pml(treeNJ, data=phang.align)

  flog.info('GTR...')
  fitGTR <- update(fit, k=4, inv=0.2)
  # fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "stochastic", control = pml.control(epsilon = 1e-08, maxit = 10,trace = 0))

  # detach("package:phangorn", unload=TRUE)
  flog.info('Done.')

  tree <- fitGTR$tree

  if(!dir.exists(output)){
    flog.info('Creating output directory...')
    dir.create(output)
    flog.info('Done.')
  }
  flog.info('Saving R objects.')
  save(tree, file=paste(output,'/robjects.Rdata',sep=''))
  flog.info('Finish.')

  if(returnval){return(tree)}


}
