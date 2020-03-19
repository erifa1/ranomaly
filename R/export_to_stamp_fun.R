#' Export to stamp
#'
#' Export 2 text file to use with STAMP
#'
#' @param dada_res output from dada2_fun
#'
#' @return Return raw otu table in phyloseq object.
#' @import phyloseq
#' @import psadd
#' @import futile.logger
#'
#' @export


# Decontam Function

export_to_stamp_fun <- function(data = data, output = "./stamp/", correc = TRUE){


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
  flog.info('Saving ...')
  write.table(stamp.table,paste(output,'/table_stamp.tsv',sep=''),sep="\t",row.names=FALSE, quote=FALSE)
  write.table(sample_data(data),paste(output,'/meta_stamp.tsv',sep=''),sep="\t",row.names=FALSE, quote=FALSE)
  flog.info('Done.')

  #Replace most common french special characters
  if(correc==TRUE){
    flog.info('Correcting special characters ...')
    command1 <- paste('sed "s/é/e/g;s/è/e/g;s/É/E/g;s/à/a/g;s/ê/e/g;s/â/a/g;s/û/u/g;s/ï/i/g;s/°//g;s/ô/o/g" ',output,'/meta_stamp.tsv > ',output,'meta_stampOK.tsv',sep="")
    invisible(system(command1, intern=TRUE))
  }

  flog.info('Finish.')

}
