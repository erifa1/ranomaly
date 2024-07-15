#' Generate Phyloseq Object
#'
#'
#'
#' @param dada_res Results of dada2_fun()
#' @param tax.table Results of assign_taxo_fun()
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
#' @import tools
#' @import vroom
#' @import here
#' @import readxl
#' @import readr
#'
#' @export


# Generate Phyloseq

generate_phyloseq_fun <- function(dada_res = dada_res, tax.table = tax.table, tree = NULL, metadata = "",
                                  output = "./phyloseq/", verbose = 1, returnval = TRUE){


  if(verbose == 3){
    invisible(flog.threshold(DEBUG))
  } else {
    invisible(flog.threshold(INFO))
  }

  flog.info("Loading sample metadata..")
  if(is.data.frame(metadata)){
    sampledata <- metadata
  }else if(tools::file_ext(metadata) %in% c('xls', 'xlsx')){
    sampledata <- readxl::read_excel(path=metadata, sheet=1, col_names=T)
  } else if (tools::file_ext(metadata) %in% c('csv', 'tsv')){
    sampledata <- vroom::vroom(file=metadata, locale = readr::locale(decimal_mark = ",", encoding = "UTF-8"), show_col_types = FALSE)
  } 
  # sampledata <- read.table(metadata, sep="\t",header=TRUE)
  sampledata <- as.data.frame(sampledata)
  rownames(sampledata) <- sampledata$sample.id

  sample.metadata <- sample_data(sampledata)
  flog.info('Done.')

  if(length(setdiff(colnames(dada_res$otu.table), sampledata$sample.id)) > 0){
    message("ERROR: sample names in otu table differs from sample metadata.")
    print(setdiff(colnames(dada_res$otu.table), sampledata$sample.id))
    stop("Please check metadata table.",call. = FALSE)
  }

  if(length(setdiff(sampledata$sample.id, colnames(dada_res$otu.table))) > 0){
    message("WARNING: One or more samples in metadata are not present in otu table:")
    print(setdiff(sampledata$sample.id, colnames(dada_res$otu.table)))

    sample.metadata <- sample.metadata[colnames(dada_res$otu.table),]
  }

  flog.info("Sequences..")
  sequences <- getSequences(dada_res$seqtab.nochim)
  flog.debug('Generate MD5 ids...')
  names(sequences) <- sapply(sequences,digest,algo='md5')

  if(!is.null(tax.table)){

    if(is.null(tree)){
      flog.info('Building phyloseq object without tree...')

      data <- phyloseq(dada_res$otu.table[,row.names(sample.metadata)], tax_table(as.matrix(tax.table)), sample.metadata, DNAStringSet(sequences))
    } else{
      flog.info('Building phyloseq object with tree...')
      data <- phyloseq(dada_res$otu.table[,row.names(sample.metadata)], tax_table(as.matrix(tax.table)), sample.metadata, phy_tree(tree), DNAStringSet(sequences))
    }

  }else{

     if(is.null(tree)){
      flog.info('Building phyloseq object without tree...')
      data <- phyloseq(dada_res$otu.table[,row.names(sample.metadata)], sample.metadata, DNAStringSet(sequences))
    } else{
      flog.info('Building phyloseq object with tree...')
      data <- phyloseq(dada_res$otu.table[,row.names(sample.metadata)], sample.metadata, phy_tree(tree), DNAStringSet(sequences))
    }


  }



  flog.info('Done.')
  data_rel <- transform_sample_counts(data, function(x) x / sum(x) )

  if(!dir.exists(output)){
    flog.info('Creating output directory...')
    dir.create(output, recursive = TRUE, showWarnings = FALSE)
    flog.info('Done.')
  }
  flog.info('Saving R objects.')
  save(data, file=paste(output,'/robjects.Rdata',sep=''))
  save(data_rel, file=paste(output,'/robjects_relative.Rdata',sep=''))
  flog.info('Finish.')

  if(returnval){return(data)}

}
