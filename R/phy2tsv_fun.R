#' Phyloseq to TSV
#'
#' Export otu table with taxonomy and sequences, + metadata in separated file.
#'
#' @param data output from decontam or generate_phyloseq
#' @param output Output directory
#' @param rank Taxonomic rank used to agglomerate otu table.
#'
#' @return Export tabulated otu table and metadata
#'
#' @import phyloseq
#'
#' @export


# Decontam Function

phy2tsv_fun <- function(data = data, output = "./tsv_table/", rank = "ASV"){

  if(!dir.exists(output)){
    dir.create(output, recursive=TRUE)
  }


  flog.info('Generating table ...')
  stable <- as.data.frame(sample_data(data))
  write.table(stable, paste(output,"/metadata_table.csv",sep=""), sep="\t", quote=FALSE, col.names=NA)
  if(any(rank == rank_names(data))){
    flog.info(paste(rank,' ...'))
    data_genus <- tax_glom(data, rank)
    ttable <- data_genus@tax_table@.Data
    otable <- as.data.frame(otu_table(data_genus))
    # refseq1 <- as.data.frame(refseq(data_genus)); names(refseq1)="seq"
  }else{
    if(rank=="ASV"){
      flog.info(paste('ASV ...'))
      ttable <- data@tax_table@.Data
      otable <- as.data.frame(otu_table(data))
      refseq1 <- as.data.frame(refseq(data)); names(refseq1)="seq"

    }else{flog.info('Choose rank name among:'); print(rank_names(data))}
  }
  if(!any(rownames(ttable) == rownames(otable))){flog.info("Different order in otu table and tax table");quit()}
  TT = cbind(otable,ttable,refseq1)
  # TT = cbind(otable,ttable)
  flog.info(paste('Saving ...'))
  write.table(TT, paste(output,"/otu_table_",rank,".csv",sep=""), sep="\t", quote=FALSE, col.names=NA)
  flog.info(paste('Done ...'))

}
