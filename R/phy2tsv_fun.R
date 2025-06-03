#' Phyloseq to TSV
#'
#' Export otu table with taxonomy and sequences, + metadata in separated file.
#'
#' @param data a phyloseq object (output from decontam or generate_phyloseq)
#' @param output Output directory
#' @param rank Taxonomy rank to merge features that have same taxonomy at a certain taxonomic rank (among rank_names(data), or 'ASV' for no glom)
#' @param relative Output table with relative abundances. 
#'
#' @return Export tabulated otu table and metadata
#'
#' @import phyloseq
#'
#' @export


# Decontam Function

phy2tsv_fun <- function(data = data, output = "./tsv_table/", rank = "ASV", relative = FALSE){

  if(!dir.exists(output)){
    dir.create(output, recursive=TRUE)
  }


  flog.info('Generating table ...')
  stable <- as.data.frame(as.matrix(sample_data(data)))
  write.table(stable, paste(output,"/metadata_table.csv",sep=""), sep="\t", quote=FALSE, col.names=NA)
  if(any(rank == rank_names(data))){
    flog.info(paste(rank,' ...'))
    data_genus <- tax_glom(data, rank)

    if(relative){
      normf = function(x){ x/sum(x) }
      data_genus <- transform_sample_counts(data_genus, normf)
    }

    ttable <- data_genus@tax_table@.Data
    otable <- as.data.frame(otu_table(data_genus))
    if(!is.null(data@refseq)){
      refseq1 <- as.data.frame(refseq(data_genus)); names(refseq1)="seq"
    }
  }else{
    if(rank=="ASV"){
      flog.info(paste('ASV ...'))
      if(relative){
        normf = function(x){ x/sum(x) }
        data <- transform_sample_counts(data, normf)
      }      
      ttable <- data@tax_table@.Data
      otable <- as.data.frame(otu_table(data))
      if(!is.null(data@refseq)){
        refseq1 <- as.data.frame(refseq(data)); names(refseq1)="seq"
      }
    }else{flog.info('Choose rank name among:'); print(rank_names(data))}
  }
  if(!any(rownames(ttable) == rownames(otable))){flog.info("Different order in otu table and tax table");quit()}
  if(!is.null(data@refseq)){
    TT = cbind(otable,ttable,refseq1)
  }else{
    TT = cbind(otable,ttable)
  }

  flog.info(paste('Saving ...'))
  write.table(TT, paste(output,"/otu_table_",rank,".csv",sep=""), sep="\t", quote=FALSE, col.names=NA)
  flog.info(paste('Done ...'))

}



#' taxa_counts
#'
#' Function to count the number of taxa at each rank on a phyloseq object
#'
#' @param data a phyloseq object
#' @param except exclude taxa (eg. unclassified|unknown)
#'
#' @return List with counts and details
#'
#' @import phyloseq
#'
#' @export

# taxa_nums function
taxa_counts <- function(data = NULL, except = NULL){

  if(is.null(data)|class(data)!="phyloseq"){
    stop("Phyloseq object required...")
  }

  ttable = tax_table(data)

  list1 = list()
  for(i in 1:7){
    if(is.null(except)){
        list1[[i]] = as.character(unique(ttable[,i]))
        }else{
        list1[[i]] = as.character(unique(ttable[grep(except, ttable[,i], value = FALSE, invert = TRUE),i]))
        }
  }
  list1[[8]] = as.character(taxa_names(data))
  names(list1) = c(rank_names(data), "ASV")

  list2 = list()
  list2$counts = sapply(list1, length)
  list2$details = list1

  return(list2)

}