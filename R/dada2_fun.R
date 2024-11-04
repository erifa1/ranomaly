#' DADA2 process (dada2_fun())
#'
#' Processing DADA2 algorithm on raw sequences (IonTorrent or Illumina / single or paired end), return raw otu table with representative sequence of ASV.
#'
#' @param path Read files folder path
#' @param outpath output .Rdata file name
#' @param cutadapt Use of cutadapt to trim primers based on their sequences, f_ and r_primer are needed (ambiguous nucleotides allowed)
#' @param maxEE Maximum expected error in reads (in filterAndTrim function) (default: 5)
#' @param dadapool option for dada function (FALSE, TRUE or "pseudo"), default is "pseudo". See ? dada.
#' @param f_trunclen Forward read truncate length (only for paired end 16S)
#' @param r_trunclen Reverse read truncate length (only for paired end 16S)
#' @param f_primer Forward primer sequence (mandatory if cutadapt = TRUE),  5' -> 3'.
#' @param r_primer Reverse primer sequence (mandatory if cutadapt = TRUE), 5' -> 3'.
#' @param plot Plot all test or not
#' @param compress Reads files are compressed (.gz)
#' @param extension Paired reads extension files for R1. (default: _R1.fastq), use "_R1.fastq.gz" for compressed files.
#' @param verbose Verbose level. (1: quiet, 3: verbal)
#' @param torrent_single Boolean to choose between Illumina Paired End SOP or Torrent Single End SOP. default: FALSE
#' @param trim_l Trim left size, to trim forward primer without cutadapt.
#' @param trim_r Trim right size, to trim reverse primer without cutadapt.
#' @param returnval Boolean to return values in console or not.
#' @param paired Boolean for Illumina Paired End Reads.
#' @param orient_torrent Forward primer sequence to orient all reads to same strand (only unambiguous nucleotides).
#' @param n_cpu Number of CPU to use in multithread.
#'
#' @return Return raw otu table in phyloseq object and export it in an Rdata file.
#'
#' @examples 
#' \dontrun{
#' See https://forgemia.inra.fr/umrf/ranomaly/-/wikis/home#dada2-usage-according-to-raw-data-type
#' }
#' @import dada2
#' @import psadd
#' @import ShortRead
#' @import ggplot2
#' @import futile.logger
#' @import phyloseq
#' @importFrom stringr str_remove str_replace
#' @export

# DADA2 function

dada2_fun <- function(path = "", outpath = "./dada2_out/", cutadapt = FALSE, maxEE = 5, f_trunclen = 240, r_trunclen = 240, dadapool = "pseudo",
                      f_primer = "GCATCGATGAAGAACGCAGC", r_primer = "TCCTCCGCTTWTTGWTWTGC", plot = FALSE, compress = FALSE, extension = "_R1.fastq",
                      verbose = 1, torrent_single = FALSE,returnval = TRUE, paired = TRUE, trim_l=15, trim_r=0, orient_torrent = NULL, n_cpu=6){
  if(torrent_single == TRUE & is.null(orient_torrent)){stop("Need forward primer to orient TORRENT reads...")}

  if(verbose == 3){
    invisible(flog.threshold(DEBUG))
  } else {
    invisible(flog.threshold(INFO))
  }

  flog.info('Reads path :')
  flog.info(path)
  if( length(list.files(path, full.names = TRUE)) == 0){
    stop(glue::glue("No files in the input directory, or this directory does not exist. ({path})"))
  }

  if(summary(file(dir(path,full.names=TRUE)[1]))$class == "gzfile" & compress == FALSE){
    stop(glue::glue("Files are compressed, please set compress = TRUE."))
  }

  flog.info("Creating directory.")
  if(!dir.exists(outpath)){
    dir.create(outpath, recursive = TRUE)
  }
  flog.info('Done.')

  flog.info("Extension...")
  extension2 = stringr::str_replace(extension, "1", "2")


  if(paired == TRUE){
    flog.info('###ILLUMINA PAIRED END SOP')

  if(compress==TRUE & length(grep(".gz", extension)) == 0){
    stop(glue::glue("Argument compress = TRUE, please add '.gz' to extension argument ({extension}.gz)."))
  }

    flog.info('Loading files...')
    fnFs <- sort(list.files(path, pattern = extension, full.names = TRUE))
    fnRs <- sort(list.files(path, pattern = extension2, full.names = TRUE))

    rawCounts <- count_seq(path, pattern = ".*R1.fastq.*")

    flog.debug("File list...")
    flog.debug(length(fnFs))
    flog.debug(length(fnRs))
    flog.info('Done.')

    if(plot){
      flog.info('Quality plot ...')
      pf <- plotQualityProfile(fnFs[1:2])
      pr <- plotQualityProfile(fnRs[1:2])
      ggsave(paste(outpath,'/qual_plot_f.png',sep=''), plot=pf)
      ggsave(paste(outpath,'/qual_plot_r.png',sep=''), plot=pr)
      flog.info('Done.')
    }




    if(cutadapt){

      flog.info(glue::glue('DADA2 Trim primers with cutadapt based on primers sequences ...'))

      # Primers sequences
      FWD <- f_primer
      REV <- r_primer

      # Generate all primers orientations
      allOrients <- function(primer) {
        # Create all orientations of the input sequence
        # require(Biostrings)
        dna <- Biostrings::DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
        orients <- c(Forward = dna, Complement = IRanges::reverse(reverseComplement(dna)), Reverse = IRanges::reverse(dna),  # bug with complement() function
        RevComp = reverseComplement(dna))
        return(sapply(orients, toString))  # Convert back to character vector
      }
      flog.info('Primer orientation')
      FWD.orients <- allOrients(FWD)
      REV.orients <- allOrients(REV)
      #FWD.orients

      flog.info('Remove Ns from reads')
      # Remove Ns from reads
      fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
      fnRs.filtN <- file.path(path, "filtN", basename(fnRs))

      # print(fnFs.filtN)
      # print(fnRs.filtN)

      flog.info('filterAndTrim...')
      if(! dir.exists(paste(path,'/filtN',sep=''))){
        filterAndTrim(fwd = fnFs, filt = fnFs.filtN, rev = fnRs, filt.rev = fnRs.filtN, maxN = 0, multithread=n_cpu, verbose=TRUE, rm.phix = TRUE, compress=compress)
      }else{
        flog.info('Filtered files exist, skipping...')
      }

      flog.info('Done.')

      # Search primers in reads.
      primerHits <- function(primer, fn) {
        # Counts number of reads in which the primer is found
        nhits <- Biostrings::vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
        return(sum(nhits > 0))
      }

      flog.info('Primer hits: ')

      primerhits <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]])
      )
      print(primerhits)

      flog.info('Done.')
      # Remove primers
      flog.info('Runnning cutadapt...')

      cutadapt <- system("which cutadapt", intern=TRUE)   # CHANGE ME to the cutadapt path on your machine
      flog.debug(cutadapt)
      path.cut <- file.path(path, "cutadapt")
      fnFs.cut <- file.path(path.cut, basename(fnFs))
      flog.debug(head(fnFs))
      fnRs.cut <- file.path(path.cut, basename(fnRs))
      if(!dir.exists(path.cut)){
        dir.create(path.cut)
        FWD.RC <- dada2:::rc(FWD)
        REV.RC <- dada2:::rc(REV)
        # Trim FWD and the reverse-complement of REV off of R1 (forward reads)
        R1.flags <- paste("-g", FWD, "-a", REV.RC)
        # Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
        R2.flags <- paste("-G", REV, "-A", FWD.RC)
        # Run Cutadapt
        pb <- txtProgressBar(min = 0, max = length(fnFs), style = 3)
        for(i in seq_along(fnFs)) {
          setTxtProgressBar(pb, i)
          if(verbose == 3){
            system2(cutadapt, args = c(R1.flags, R2.flags, "--discard-untrimmed", "-n", 2, # -n 2 required to remove FWD and REV from reads
            "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
            fnFs.filtN[i], fnRs.filtN[i]), stdout="", stderr="") # input files
          } else{
            system2(cutadapt, args = c(R1.flags, R2.flags, "--discard-untrimmed", "-n", 2, # -n 2 required to remove FWD and REV from reads
            "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
            fnFs.filtN[i], fnRs.filtN[i]), stdout=NULL, stderr=NULL) # input files
          }
        }
        close(pb)
      } else{
        flog.info('Cutadapt files exists. Skipping...')
      }

      flog.info('Done.')
      # Sanity check
      flog.info('Primer left: ')
      primerhits2 <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))
      print(primerhits2)

      # Forward and reverse fastq filenames have the format:
      cutFs <- sort(list.files(path.cut, pattern = extension, full.names = TRUE))
      cutRs <- sort(list.files(path.cut, pattern = extension2, full.names = TRUE))

      # Extract sample names, assuming filenames have format:

      sample.names <- unname(sapply(cutFs, get.sample.name))
      #head(sample.names)

      if(compress){
        filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
        filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
      }else{
        filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
        filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))
      }

      names(filtFs) <- sample.names
      names(filtRs) <- sample.names

      flog.debug(length(cutFs))
      flog.debug(length(filtFs))
      flog.debug(length(cutRs))
      flog.debug(length(filtRs))

      out0 <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(maxEE,maxEE), #truncLen=c(f_trunclen,r_trunclen),
      truncQ = 2, minLen = 50, compress = compress, multithread=n_cpu)  # on windows, set multithread = FALSE
      #head(out)
      row.names(out0) = sample.names
      out <- as.data.frame(out0) %>%tibble::rownames_to_column(var = "sample.id")


      filtFs_out <- grep("F_filt",list.files(glue::glue("{path}/filtered/"), full.names = TRUE), value = TRUE)
      filtRs_out <- grep("R_filt",list.files(glue::glue("{path}/filtered/"), full.names = TRUE), value = TRUE)

      # samples_names <- grep(".fastq",list.files("./data_arch_links/", full.names = FALSE), value = TRUE)

    } else {

      flog.info(glue::glue('DADA2 Trim primers on {trim_l} bases ...'))

      flog.debug(length(fnFs))
      flog.debug(length(fnRs))

      sample.names <- unname(sapply(fnFs, get.sample.name))
      flog.debug(head(sample.names))
      # sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

      if(compress){
        filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
        filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
      }else{
        filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
        filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))
      }


      names(filtFs) <- sample.names
      names(filtRs) <- sample.names


      flog.info('Filtering reads...')

      out0 <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(f_trunclen,r_trunclen),
      maxN=0, maxEE=c(maxEE,maxEE), truncQ=2, rm.phix=TRUE, trimLeft=trim_l,
      compress=compress, multithread=n_cpu)
      row.names(out0) = sample.names
      out <- as.data.frame(out0) %>%tibble::rownames_to_column(var = "sample.id")


      filtFs_out <- grep("F_filt",list.files(glue::glue("{path}/filtered/"), full.names = TRUE), value = TRUE)
      filtRs_out <- grep("R_filt",list.files(glue::glue("{path}/filtered/"), full.names = TRUE), value = TRUE)

      # samples_names <- grep(".fastq",list.files(path, full.names = FALSE), value = TRUE)
      # print(samples_names)

      flog.info('Done.')
    }

    #COMMON
    flog.info('Learning error model...')
    errF <- learnErrors(filtFs_out, multithread=n_cpu)
    errR <- learnErrors(filtRs_out, multithread=n_cpu)
    flog.info('Done.')


    if(plot){
      flog.info('Plotting 2 ...')
      pf2 <- plotErrors(errF, nominalQ=TRUE)
      pr2 <- plotErrors(errR, nominalQ=TRUE)
      # print("plotOK")
      ggsave(paste(outpath,'/err_plot_f.png',sep=''), plot=pf2)
      ggsave(paste(outpath,'/err_plot_f.png',sep=''), plot=pr2)
      flog.info('Done.')
    }
    # save.image("debug.rdata")

    mergers <- vector("list", length(sample.names))
    names(mergers) <- sample.names;
    stockFs=NULL; stockRs=NULL
    getN <- function(x) sum(getUniques(x))

    for(sam in sample.names) {
      flog.info(paste('Processing sample ',sam))
      flog.info('Dereplicating fastq...')
      derepFs <- derepFastq(filtFs[[sam]], verbose=TRUE)
      derepRs <- derepFastq(filtRs[[sam]], verbose=TRUE)
      flog.info('Done.')
      flog.info('dada2...')
      dadaFs <- dada(derepFs, err=errF, multithread=n_cpu, pool=dadapool, selfConsist=FALSE)
      stockFs <- c(stockFs, getN(dadaFs))
      dadaRs <- dada(derepRs, err=errR, multithread=n_cpu, pool=dadapool, selfConsist=FALSE)
      stockRs <- c(stockRs,getN(dadaRs))
      flog.info('Done.')
      flog.info('Merging pairs...')
      merger <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=FALSE)
      flog.info('Done.')
      mergers[[sam]] <- merger
    }

    #
    #
    # flog.info('Dereplicating fastq...')
    # derepFs <- derepFastq(filtFs_out, verbose=TRUE)
    # derepRs <- derepFastq(filtRs_out, verbose=TRUE)
    #
    # flog.info('Done.')
    # flog.info('dada2...')
    # dadaFs <- dada(derepFs, err=errF, multithread=n_cpu, pool=dadapool, selfConsist=FALSE)
    # stockFs <- sapply(dadaFs, getN)
    # dadaRs <- dada(derepRs, err=errR, multithread=n_cpu, pool=dadapool, selfConsist=FALSE)
    # stockRs <- sapply(dadaRs, getN)
    # flog.info('Done.')
    # flog.info('Merging pairs...')
    # mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=FALSE)
    # flog.info('Done.')


    seqtab <- makeSequenceTable(mergers)
    flog.info('Removing chimeras...')
    seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=n_cpu, verbose=TRUE)
    flog.debug(sum(seqtab.nochim))
    flog.debug(sum(seqtab))

    flog.debug(sum(seqtab.nochim)/sum(seqtab))
    flog.info('Done.')


    track0 <- cbind.data.frame(stockFs, stockRs, sapply(mergers, getN), rowSums(seqtab.nochim))
    if(compress==TRUE){
      rownames(track0) <- stringr::str_remove(rownames(track0), "_F_filt.fastq.gz")
    }else{
      rownames(track0) <- stringr::str_remove(rownames(track0), "_F_filt.fastq")
    }

    # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

    track <- track0 %>% tibble::rownames_to_column(var="sample.id")


    final_track <- out %>% dplyr::left_join(y = track, by = "sample.id") %>% dplyr::mutate(input = rawCounts, .after = 1)
    colnames(final_track) <- c("sample.id", "rawcounts", "primer filtered", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
    head(final_track)

    write.table(final_track, paste(outpath,"/read_tracking.csv",sep=''), sep="\t", row.names=FALSE, quote=FALSE)


    seqtab.export <- seqtab.nochim
    colnames(seqtab.export) <- sapply(colnames(seqtab.export), digest::digest, algo="md5")

    otu.table <- phyloseq::otu_table(t(seqtab.export), taxa_are_rows = TRUE)
    if(compress==TRUE){
      colnames(otu.table) <- stringr::str_remove(colnames(otu.table), "_F_filt.fastq.gz")
    }else{
      colnames(otu.table) <- stringr::str_remove(colnames(otu.table), "_F_filt.fastq")
    }

    flog.info('Writing raw tables.')
    write.table(cbind(t(seqtab.export), "Sequence" = colnames(seqtab.nochim)), paste(outpath,"/raw_otu-table.csv",sep=''), sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

    flog.info('Writing fasta.')
    uniquesToFasta(seqtab.nochim, fout=paste(outpath,'/rep-seqs.fna',sep=''), ids=colnames(seqtab.export))

    dada_res=list()
    dada_res$seqtab.nochim = seqtab.nochim
    dada_res$seqtab.export = seqtab.export
    dada_res$otu.table = otu.table
    flog.info('Saving R objects.')
    save(dada_res, file=paste(outpath,'/robjects.Rdata',sep=''))
    flog.info('Finish.')

    if(returnval) {return(dada_res)}


  }else{
    flog.info('### SINGLE END SOP')
    rawCounts <- count_seq(path, pattern = ".*fastq.*")

    if(compress==TRUE){
      fnFs <- sort(list.files(path, pattern = ".fastq.gz", full.names = TRUE))
      fastq.names <- sort(list.files(path, pattern = ".fastq.gz", full.names = FALSE))
    }else{
      fnFs <- sort(list.files(path, pattern = ".fastq", full.names = TRUE))
      fastq.names <- sort(list.files(path, pattern = ".fastq", full.names = FALSE))
    }
    #sample name
    sample.names <- unname(sapply(fastq.names, get.sample.name))
    if(compress==TRUE){
      filtFs <- file.path(path, "filtered", paste0(sample.names, "_filt.fastq.gz"))
    }else{
      filtFs <- file.path(path, "filtered", paste0(sample.names, "_filt.fastq"))
    }

    # Cutadapt

      if(cutadapt){
          flog.info(glue::glue('DADA2 Trim primers with cutadapt based on primers sequences ...'))

      # Primers sequences
      FWD <- f_primer
      REV <- r_primer

      # Generate all primers orientations
      allOrients <- function(primer) {
        # Create all orientations of the input sequence
        # require(Biostrings)
        dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
        orients <- c(Forward = dna, Complement = IRanges::reverse(reverseComplement(dna)), Reverse = IRanges::reverse(dna),  # bug with complement() function
        RevComp = reverseComplement(dna))
        return(sapply(orients, toString))  # Convert back to character vector
      }
      flog.info('Primer orientation')
      FWD.orients <- allOrients(FWD)
      REV.orients <- allOrients(REV)
      #FWD.orients

      flog.info('Remove Ns from reads')
      # Remove Ns from reads
      fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory


      # print(fnFs.filtN)

      flog.info('filterAndTrim & orient reads...')
      if(! dir.exists(paste(path,'/filtN',sep=''))){
        filterAndTrim(fwd = fnFs, filt = fnFs.filtN, maxN = 0, multithread=n_cpu, verbose=TRUE,
          rm.phix = TRUE, compress=compress, orient.fwd = orient_torrent)
      }else{
        flog.info('Filtered files exist, skipping...')
      }

      flog.info('Done.')

      # Search primers in reads.
      primerHits <- function(primer, fn) {
        # Counts number of reads in which the primer is found
        nhits <- Biostrings::vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
        return(sum(nhits > 0))
      }

      flog.info('Primer hits: ')
      primerhits <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]])
      )
      print(primerhits)

      flog.info('Done.')
      # Remove primers
      flog.info('Runnning cutadapt...')

      cutadapt <- system("which cutadapt", intern=TRUE)   # CHANGE ME to the cutadapt path on your machine
      flog.debug(cutadapt)
      path.cut <- file.path(path, "cutadapt")
      fnFs.cut <- file.path(path.cut, basename(fnFs))
      flog.debug(head(fnFs))

      if(!dir.exists(path.cut)){
        dir.create(path.cut)
        FWD.RC <- dada2:::rc(FWD)
        REV.RC <- dada2:::rc(REV)
        # Trim FWD and the reverse-complement of REV off of R1 (forward reads)
        R1.flags <- paste("-g", FWD, "-a", REV.RC, "-m 100")

        # Run Cutadapt
        pb <- txtProgressBar(min = 0, max = length(fnFs), style = 3)
        for(i in seq_along(fnFs)) {
          setTxtProgressBar(pb, i)
          if(verbose == 3){
            system2(cutadapt, args = c(R1.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
            "-o", fnFs.cut[i], # output files
            fnFs.filtN[i], "--rc", "--discard-untrimmed"), stdout="", stderr="") # input files
          } else{
            system2(cutadapt, args = c(R1.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
            "-o", fnFs.cut[i], # output files
            fnFs.filtN[i], "--rc", "--discard-untrimmed"), stdout=NULL, stderr=NULL) # input files
          }
        }
        close(pb)
      } else{
        flog.info('Cutadapt files exists. Skipping...')
      }

      flog.info('Done.')
      # Sanity check
      flog.info('Primer left: ')
      primerhits2 <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]))
      print(primerhits2)

      if(compress){
        cutFs <- sort(list.files(path.cut, pattern = ".fastq.gz", full.names = TRUE))
      }else{
        cutFs <- sort(list.files(path.cut, pattern = ".fastq", full.names = TRUE))
      }

      sample.names <- unname(sapply(cutFs, get.sample.name))

      if(compress){
        filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
      }else{
        filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
      }

      names(filtFs) <- sample.names
      # print(filtFs)

      flog.debug(length(cutFs))
      flog.debug(length(filtFs))
      flog.info('filterAndTrim...')
      out0 <- filterAndTrim(fwd = cutFs, filt = filtFs, maxN = 0, multithread=n_cpu, rm.phix = TRUE,
        , maxEE = maxEE , minLen = 100, compress=compress)


      row.names(out0) = sample.names
      out <- as.data.frame(out0) %>%tibble::rownames_to_column(var = "sample.id")
      flog.info('Done.')

      filtFs_out <- grep("F_filt",list.files(glue::glue("{path}/filtered/"), full.names = TRUE), value = TRUE)

    }else{
      flog.info(glue::glue('DADA2 Trim primers on {trim_l}(FWD) and {trim_r}(REV) bases ...'))

      if(torrent_single == TRUE){
        out0 <- filterAndTrim(fwd = fnFs, filt = filtFs, maxN = 0, multithread=n_cpu, verbose=TRUE, rm.phix = TRUE,
          , maxEE = maxEE , minLen = 100, compress=compress, trimLeft=trim_l, trimRight=trim_r, orient.fwd = orient_torrent)
      }else{
        out0 <- filterAndTrim(fwd = fnFs, filt = filtFs, maxN = 0, multithread=n_cpu, verbose=TRUE, rm.phix = TRUE,
          , maxEE = maxEE , minLen = 100, compress=compress, trimLeft=trim_l, trimRight=trim_r )
      }
      row.names(out0) = sample.names
      out <- as.data.frame(out0) %>%tibble::rownames_to_column(var = "sample.id")
      flog.info('Done.')
    }

    flog.info('Learning error model...')
    errF <- learnErrors(filtFs, multithread=n_cpu)
    flog.info('Done.')

    if(plot){
      flog.info('Plotting 2 ...')
      pf2 <- plotErrors(errF, nominalQ=TRUE)
      ggsave(paste(outpath,'/err_plot_f.png',sep=''), plot=pf2)
      flog.info('Done.')
    }

    stockFs=NULL
    getN <- function(x) sum(getUniques(x))


    flog.info('Dereplicating fastq...')
    if(compress==TRUE){
      filtFs <- sort(list.files(file.path(path, "filtered/"), pattern = ".fastq.gz$", full.names = TRUE))
      name1 <- sort(list.files(file.path(path, "filtered/"), pattern = ".fastq.gz$", full.names = FALSE))
      sample.names <- names(filtFs) <- sapply( stringr::str_split(name1, "_"), "[[", 1)
    } else{
      filtFs <- sort(list.files(file.path(path, "filtered/"), pattern = ".fastq$", full.names = TRUE))
      name1 <- sort(list.files(file.path(path, "filtered/"), pattern = ".fastq$", full.names = FALSE))
      sample.names <- names(filtFs) <- sapply( stringr::str_split(name1, "_"), "[[", 1)
    }

    # names(filtFs) <- sample.names
    derepFs <- derepFastq(filtFs, verbose=TRUE)
    flog.info('Done.')

    flog.info('dada2...')

    if(torrent_single == TRUE){
      dadaFs <- dada(derepFs, err=errF, multithread=n_cpu, pool=dadapool, selfConsist=FALSE, HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)
    }
    else{
      dadaFs <- dada(derepFs, err=errF, multithread=n_cpu, pool=dadapool, selfConsist=FALSE)
    }


    if(length(filtFs)<2){
      stockFs <- getUniques(dadaFs)
    }else{
      stockFs <- sapply(dadaFs, getN)
    }

    flog.info('Done.')

    seqtab <- makeSequenceTable(dadaFs)

    flog.info('Removing chimeras...')
    seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=n_cpu, verbose=TRUE)

    if(length(filtFs)<2){
      track <- c(out, sum(stockFs), rowSums(seqtab.nochim)[1])
      names(track) <- c("input", "filtered", "denoisedF", "nonchim")
    }else{
      rownames(out) <- out$sample.id
      nn <- row.names(out)
      # print(names(stockFs))
      track <- cbind.data.frame(out, stockFs[nn], rowSums(seqtab.nochim)[nn]) %>% dplyr::mutate(input = rawCounts, .after = 1)
      colnames(track) <- c("sample.id", "rawcounts", "primer filtered", "filtered", "denoisedF", "nonchim")
      head(track)
    }
    flog.info('Writing table ...')
    write.table(track, paste(outpath,"/read_tracking.csv",sep=''), sep="\t", row.names=FALSE, quote=FALSE)

    seqtab.export <- seqtab.nochim
    colnames(seqtab.export) <- sapply(colnames(seqtab.export), digest::digest, algo="md5")

    otu.table <- phyloseq::otu_table(t(seqtab.export), taxa_are_rows = TRUE)
    # colnames(otu.table) = sample.names

    flog.info('Writing raw tables.')
    write.table(cbind(t(seqtab.export), "Sequence" = colnames(seqtab.nochim)), paste(outpath,"/raw_otu-table.csv",sep=''), sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

    flog.info('Writing fasta.')
    uniquesToFasta(seqtab.nochim, fout=paste(outpath,'/rep-seqs.fna',sep=''), ids=colnames(seqtab.export))

    dada_res=list()
    dada_res$seqtab.nochim = seqtab.nochim
    dada_res$seqtab.export = seqtab.export
    dada_res$otu.table = otu.table
    flog.info('Saving R objects.')
    save(dada_res, file=paste(outpath,'/robjects.Rdata',sep=''))


    flog.info('Finish.')
    if(returnval) {return(dada_res)}


  }


}


#' Get sample name
#'
#' Extract sample names, assuming filenames have format names_miscinfo.fastq
#'
#' @param fname a parameter
#' @keywords internal

get.sample.name <- function(fname){
  tt <- strsplit(basename(fname), "_")[[1]][1]
  return(tt)
}

#' Frequency filter on dada_fun output
#'
#' Apply filter to eliminate rare ASVs, allowing to reduce the size of the dataset and to remove potential contaminants.
#' @param dada_res a dada_fun output
#' @param freq a frequency threshold to filter rare ASVs
#' @param top a number of top ASVs to keep
#' @return Filtered dada_fun output
#' 
#' @examples
#' \dontrun{
#' dada_filter(dada_res = dada_res, freq = 0.00005)
#' }
#' @export

dada_filter <- function(dada_res = dada_res, freq = 0.00005, top = NULL){

  seqtab.nochim <- dada_res$seqtab.nochim %>% t() %>% as.data.frame() %>%
      dplyr::mutate( sumASV =  apply(., 1, sum) ) %>%
      dplyr::mutate( freqASV = sumASV / sum(sumASV) ) %>%
      dplyr::filter( freqASV > 0.00005 ) %>%
      dplyr::select( -sumASV, -freqASV ) %>% t()

  if(!is.null(top)){
      seqtab.nochim <- seqtab.nochim[,1:top]
  }

  seqtab.export <- seqtab.nochim
  colnames(seqtab.export) <- sapply(colnames(seqtab.export), digest::digest, algo="md5")

  otu.table <- phyloseq::otu_table(t(seqtab.export), taxa_are_rows = TRUE)

  dada_res2 <- list(seqtab.nochim = seqtab.nochim, seqtab.export = seqtab.export, otu.table = otu.table)

  return(dada_res2)
}

