#' Update metadata
#'
#'
#' @param data a phyloseq object (output from decontam or generate_phyloseq)
#' @param metadata Tabulated metadata file path.
#' @param output Output directory
#' @param returnval Boolean to return values in console or not.
#'
#' @return Return a phyloseq object with updated metadata.
#'
#' @export


update_metadata_fun <- function(data = data, output = "./updated_physeq/", metadata = NULL, returnval = TRUE){

  if (is.null(metadata)){
    print_help(opt_parser)
    stop("You must provide metadata file path.", call.=FALSE)
  } else{
    flog.info("Loading sample metadata..")
    sampledata <- read.table(metadata, sep="\t",header=TRUE)
    rownames(sampledata) <- sampledata$sample.id
    sample.metadata <- sample_data(sampledata)
    flog.info('Done.')
  }

  if(!dir.exists(output)){
    flog.info('Creating output directory...')
    dir.create(output)
    flog.info('Done.')
  }
  print(sample.metadata)

  sample_data(data) <- sample.metadata

  save(data, file=paste(output,'/robjects.Rdata',sep=''))

  if(returnval){data}

}
