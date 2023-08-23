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
#' @import readxl


update_metadata_fun <- function(data = data, output = "./updated_physeq/", metadata = NULL, returnval = TRUE){

  if (is.null(metadata)){
    print_help(opt_parser)
    stop("You must provide metadata file path.", call.=FALSE)
  } else{
    flog.info("Loading sample metadata..")
    if(tools::file_ext(metadata) %in% c('xls', 'xlsx')){
      sampledata <- readxl::read_excel(path=metadata, sheet=1, col_names=T)
      sampledata <- as.data.frame(sampledata)
    } else if (tools::file_ext(metadata) %in% c('csv', 'tsv')){
      sampledata <- vroom::vroom(file=metadata, delim="\t", locale = readr::locale(decimal_mark = ",", encoding = "UTF-8"))
    }
    rownames(sampledata) <- sampledata[,'sample.id']

    sample.metadata <- sample_data(sampledata)
    flog.info('Done.')
  }

  if(!dir.exists(output)){
    flog.info('Creating output directory...')
    dir.create(output)
    flog.info('Done.')
  }

  sample_data(data) <- sample.metadata

  save(data, file=paste(output,'/robjects.Rdata',sep=''))

  if(returnval){data}

}
