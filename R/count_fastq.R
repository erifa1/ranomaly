#' count_seq
#'
#' Count the number of reads in fastq files, .gz files are allowed.
#' 
#' @param path Path to the fastq files directory
#' @param pattern Regular expression to match the fastq files
#'
#' @return  A vector of sequence counts
#'
#' @import Biostrings
#' @import foreach
#' @import doParallel
#' @export

count_seq <- function(path, pattern = ".*R1.*fastq.*"){
    L1 = list.files(path, pattern = pattern, full.names = TRUE)

    if(length(L1) == 0){
        stop(glue::glue("No fastq files found in the directory {path}, pattern {pattern}"))
    }

    # names of samples 
    snames <- dir(path) %>% grep(pattern ,., value = TRUE) %>% as.data.frame %>% 
    tidyr::separate(., col = 1, into = "sample", sep = "_", remove = TRUE) %>%
    as.vector %>% unlist()
    res <- NULL
    foreach (i=1:length(L1)) %dopar% {
        res[i] <- count_fastq_sequences(L1[i])
    } %>% unlist() %>% as.numeric() -> res1

    names(res1) <- snames
    return(res1)
}

#' count_fastq_sequences
#'
#' Count the number of sequences in a fastq file.

count_fastq_sequences <- function(file_path) {
  # Ouvrir le fichier en mode lecture
  con <- file(file_path, "r")
  count <- 0
  
  # Lire le fichier ligne par ligne
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    # Vérifier si la ligne est une ligne de séparateur '+'
    if (grepl("^@", line)) {
      count <- count + 1
    }
    # Sauter les trois lignes suivantes (séquence, séparateur, qualité)
    readLines(con, n = 3, warn = FALSE)
  }
  
  # Fermer le fichier
  close(con)
  
  return(count)
}
