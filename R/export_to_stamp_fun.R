#' Export to stamp
#'
#' Export 2 text file to use with STAMP
#'
#' @param data a phyloseq object (output from decontam or generate_phyloseq)
#' @param output Output directory
#' @param correc If TRUE, correct metadata to replace most common special characters, save the new file in meta_stampOK.tsv.
#'
#' @return Export 2 text files ready to use with STAMP.
#'
#' @references https://beikolab.cs.dal.ca/software/STAMP
#'
#' @import phyloseq
#' @import psadd
#' @import futile.logger
#'
#' @export


# Decontam Function

export_to_stamp_fun <- function(data = data, output = "./stamp/", correc = FALSE){


  ranknames <- colnames(tax_table(data))
  flog.info("Creating Stamp table")
  stamp.table <- merge(tax_table(data), otu_table(data), by="row.names", all.x=TRUE)
  # stamp.table["Row.names"] <- NULL
  colnames(stamp.table)[colnames(stamp.table)=="Row.names"] <- "cluster"
  stamp.table <- stamp.table[,c(ranknames,"cluster",names(stamp.table)[9:ncol(stamp.table)])]
  flog.info('Done.')

  if(!dir.exists(output)){
    flog.info('Creating output directory...')
    dir.create(output)
    flog.info('Done.')
  }

  meta_stamp <- get_variable(data)
  if(correc==TRUE){
    flog.info('Correcting special characters ...')
    meta_stamp <- data.frame(lapply(meta_stamp, remove_non_ascii), stringsAsFactors = FALSE)
    names(meta_stamp) <- remove_non_ascii(names(meta_stamp))
  }

  flog.info('Saving ...')
  write.table(stamp.table,paste(output,'/table_stamp.tsv',sep=''),sep="\t",row.names=FALSE, quote=FALSE)
  write.table(meta_stamp,paste(output,'/meta_stamp.tsv',sep=''),sep="\t",row.names=FALSE, quote=FALSE)
  flog.info('Done.')

  flog.info('Finish.')

}
