#' Create Phyloseq object from external datas.
#'
#'
#' @param otutable Tabulated ASVtable file path.
#' @param taxtable Tabulated taxonomy table file path.
#' @param seq Tabulated sequence file path.
#' @param metadata Tabulated metadata file path.
#' @param generateTree Boolean to generate the phylogenetic tree.
#' @param output Output directory
#' @param returnval Boolean to return values in console or not.
#'
#' @return Return a phyloseq object
#'
#' @export



csv2phyloseq_fun <- function(otutable = NULL, taxtable = NULL, seq = NULL, metadata = NULL, generateTree = FALSE, output = "./csv2phyloseq/", returnval = TRUE){

  flog.info("Input...")
  flog.info("...OTUtable...")
  otable <- read.table(otutable, sep="\t", h=TRUE)
  row.names(otable)=otable[,1]
  flog.info("...TAXtable...")
  ttable <- as.matrix(read.table(taxtable, sep="\t", h=TRUE, stringsAsFactors=FALSE))
  row.names(ttable)=ttable[,1]
  flog.info("...Metadatas...")
  mtable <- read.table(metadata, sep="\t", h=TRUE)
  row.names(mtable)=mtable[,1]
  flog.info("...Sequences...")
  stable <- read.table(seq, sep="\t", h=TRUE)
  sequences1 = DNAStringSet(stable[,-1])
  names(sequences1)=stable[,1]


  if(generateTree){
    flog.info('Generating tree...')
    sequences <- sequences1
    # names(sequences) <- sapply(sequences,digest,algo='md5')
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

    data <- phyloseq(otu_table(otable[,-1], taxa_are_rows=TRUE), sample_data(mtable), tax_table(ttable[,-1]), sequences1, tree) # phy_tree(tree),
  }else{
    data <- phyloseq(otu_table(otable[,-1], taxa_are_rows=TRUE), sample_data(mtable), tax_table(ttable[,-1]), sequences1)
  }

  if(!dir.exists(output)){
    flog.info('Creating output directory...')
    dir.create(output)
    flog.info('Done.')
  }
  flog.info('Saving R objects.')
  save(data, file=paste(output,'/phyloseq_object.Rdata',sep=''))

  if(returnval){return(data)}

  flog.info('Finish.')

}
