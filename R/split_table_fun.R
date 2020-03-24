#' Split phyloseq object
#'
#'
#' @param data output from decontam or generate_phyloseq
#' @param output Output directory
#' @param column1 Column name of factor to split phyloseq object with.
#'
#' @return Export new folders with splited phyloseq Robjects
#' @import phyloseq
#'
#' @export


# Decontam Function

split_table_fun <- function(data = data, output = "./", column1 = ""){

  print(data)
  vector <- levels(as.factor( data.frame(sample_data(data)[,column1])[,1]) )
  print(vector)

  data_tmp <- data
  for (var in vector){
    data <- data_tmp
    fun <- paste('data <- subset_samples(data, ',column1,'=="',var,'")',sep='')
    eval(parse(text=fun))
    if(!dir.exists(paste(output,'/split_',column1,'_',var,sep=''))){
      dir.create(paste(output,'/split_',column1,'_',var,sep=''), recursive = TRUE)
    }
    data <- prune_taxa(taxa_sums(data) >= 1, data)
    save(data, file=paste(output,'/split_',column1,'_',var,'/robjects.Rdata',sep=''))
    taxa.string <- apply(tax_table(data), 1, paste, collapse = ";")
    write.table(cbind(otu_table(data),"Consensus Lineage" = taxa.string),paste(output,'split_',column1,'_',var,"/raw_otu-table.csv",sep=''), sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

  }



}
