#' subset_fastx
#'
#' Allows subset fastq or fasta files at a given threshold. This function can convert fastq to fasta.
#'
#' @param path Path to the fastq files directory
#' @param format input format, fasta or fastq format are allowed.
#' @param outformat output format, fasta or fastq format are allowed.
#' @param output Path to the output directory
#' @param nbseq Number of sequences to output per fastq files; if NULL no subset.
#' @param ncores Number of CPU used to process.
#' @param compress set on TRUE to get .fastx.gz compressed files
#' @param verbose Print information while processing
#' @param random Pick randomly 'nbseq' sequences.
#' @param seed Fix seed with numeric to reproduce random picking (eg. 123)
#'
#' @return Output fastx files in the output directory
#'
#' @import doParallel
#' @import foreach
#' @export

subset_fastx <- function(path = NULL, format = "fastq", outformat = "fastq", output = "./subset_fastq/", nbseq = 10000, ncores = 3, compress=FALSE, verbose=FALSE, random = FALSE, seed = NULL){
  
  if(is.null(path)){stop("Require path to fastq files directory...")
  }else if(!dir.exists(path)){stop("Path does not exist...")}
  if(format != "fasta" & format != "fastq"){stop("fastq and fasta format are allowed...")}

  registerDoParallel(ncores)
  if(!dir.exists(output)){dir.create(output, recursive=TRUE)}
  L1 = list.files(path, pattern = glue::glue("*{format}*"), full.names = TRUE)
  L2 = tools::file_path_sans_ext(list.files(path, pattern = glue::glue("*{format}*"), 
    full.names = FALSE), compression = TRUE)

  X=NULL
  foreach (i=1:length(L1)) %dopar% {
    outname <- glue::glue("{output}/{L2[i]}.{outformat}")
    if (compress) {
      outname <- glue::glue("{outname}.gz")
    }
    X <- Biostrings::readDNAStringSet(L1[i], format=format, with.qualities=TRUE)
    if(verbose){print(L2[i]);print(length(X))}
    if(is.null(nbseq)){
      Biostrings::writeXStringSet(X, outname, format = outformat, compress=compress)
      return("Done")
    }

    if(length(X)>nbseq){
      if(random){
        if(!is.null(seed)){set.seed(seed)}
        rdm = sample(1:length(X),nbseq)
        Xout = X[rdm]
      }else{
        Xout = X[1:nbseq]
      }
    }else{Xout = X}
    Biostrings::writeXStringSet(Xout, outname, format = outformat, compress=compress)
  }
  return("Done")

}
