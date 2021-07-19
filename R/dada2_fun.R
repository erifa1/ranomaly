#' DADA2 process (dada2_fun())
#'
#' Processing DADA2 algorithm on raw sequences, return raw otu table with representative sequence of ASV.
#'
#' @param amplicon Choose amplipcon "16S" or "ITS"
#' @param path Read files folder path
#' @param outpath output .Rdata file name
#' @param dadapool option for dada function (FALSE, TRUE or "pseudo"), default is "pseudo". See ? dada.
#' @param f_trunclen Forward read tuncate length (only for paired end 16S)
#' @param r_trunclen Reverse read tuncate length (only for paired end 16S)
#' @param f_primer Forward primer sequence (only for ITS)
#' @param r_primer Reverse primer sequence (only for ITS)
#' @param plot Plot all test or not
#' @param compress Reads files are compressed (.gz)
#' @param verbose Verbose level. (1: quiet, 3: verbal)
#' @param torrent_single Boolean to choose between Illumina Paired End SOP or Torrent Single End SOP. default: FALSE
#' @param trim_l Trim left size.
#' @param trim_r Trim right size.
#' @param returnval Boolean to return values in console or not.
#' @param paired Boolean for Illumina Paired End Reads.
#' @param orient_torrent Forward primer sequence to orient all reads to same strand.
#'
#' @return Return raw otu table in phyloseq object and export it in an Rdata file.
#'
#' @import dada2
#' @import psadd
#' @import ShortRead
#' @import Biostrings
#' @import ggplot2
#' @import futile.logger
#' @import digest
#' @import phyloseq
#' @export

# DADA2 function

dada2_fun <- function(amplicon = "16S", path = "", outpath = "./dada2_out/", f_trunclen = 240, r_trunclen = 240, dadapool = "pseudo",
                      f_primer = "GCATCGATGAAGAACGCAGC", r_primer = "TCCTCCGCTTWTTGWTWTGC", plot = FALSE, compress = FALSE, verbose = 1,
                      torrent_single = FALSE,returnval = TRUE, paired = TRUE, trim_l=15, trim_r=0, orient_torrent = NULL){
  if(torrent_single == TRUE & is.null(orient_torrent)){stop("Need forward primer to orient TORRENT reads...")}

  if(verbose == 3){
    invisible(flog.threshold(DEBUG))
  } else {
    invisible(flog.threshold(INFO))
  }

  flog.info('Read path :')
  flog.info(path)

  flog.info("Creating directory.")
  if(!dir.exists(outpath)){
    dir.create(outpath, recursive = TRUE)
  }
  flog.info('Done.')

  if(paired == TRUE){
    flog.info('###ILLUMINA PAIRED END SOP')

    if(compress==TRUE){
      flog.info('Loading compress files...')
      fnFs <- sort(list.files(path, pattern = "_R1.fastq.gz", full.names = TRUE))
      fnRs <- sort(list.files(path, pattern = "_R2.fastq.gz", full.names = TRUE))
    }else{
      flog.info('Loading flat files...')
      fnFs <- sort(list.files(path, pattern = "_R1.fastq$", full.names = TRUE))
      fnRs <- sort(list.files(path, pattern = "_R2.fastq$", full.names = TRUE))
    }

    flog.debug("File list...")
    flog.debug(length(fnFs))
    flog.debug(length(fnRs))
    flog.info('Done.')

    if(amplicon=="ITS"){
      flog.info('DADA2 ITS')

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
      fnRs.filtN <- file.path(path, "filtN", basename(fnRs))

      # print(fnFs.filtN)
      # print(fnRs.filtN)

      flog.info('filterAndTrim...')
      if(! dir.exists(paste(path,'/filtN',sep=''))){
        filterAndTrim(fwd = fnFs, filt = fnFs.filtN, rev = fnRs, filt.rev = fnRs.filtN, maxN = 0, multithread = TRUE, verbose=TRUE, rm.phix = TRUE, compress=compress)
      }else{
        flog.info('Filtered files exist, skipping...')
      }

      flog.info('Done.')

      # Search primers in reads.
      primerHits <- function(primer, fn) {
        # Counts number of reads in which the primer is found
        nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
        return(sum(nhits > 0))
      }

      flog.info('Primer hits: ')

      rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]])
      )

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
        for(i in seq_along(fnFs)) {
          if(verbose == 3){
            system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
            "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
            fnFs.filtN[i], fnRs.filtN[i]), stdout="", stderr="") # input files
          } else{
            system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
            "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
            fnFs.filtN[i], fnRs.filtN[i]), stdout=NULL, stderr=NULL) # input files
          }
        }
      } else{
        flog.info('Cutadapt files exists. Skipping...')
      }

      flog.info('Done.')
      # Sanity check
      flog.info('Primer left: ')
      rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

      # Forward and reverse fastq filenames have the format:
      cutFs <- sort(list.files(path.cut, pattern = "_R1.fastq", full.names = TRUE))
      cutRs <- sort(list.files(path.cut, pattern = "_R2.fastq", full.names = TRUE))

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
      print(head(cutFs))
      flog.debug(length(filtFs))
      print(head(filtFs))
      flog.debug(length(cutRs))
      print(head(cutRs))
      flog.debug(length(filtRs))
      print(head(filtRs))

      out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2),
      truncQ = 2, minLen = 50, compress = FALSE, multithread = TRUE)  # on windows, set multithread = FALSE
      #head(out)

    } else if(amplicon=="16S"){

      flog.info('DADA2 16S')
      # save.image("debug.rdata")
      if(plot){
        flog.info('Plotting 1 ...')
        pf <- plotQualityProfile(fnFs[1:2])
        pr <- plotQualityProfile(fnRs[1:2])
        ggsave(paste(outpath,'/qual_plot_f.png',sep=''), plot=pf)
        ggsave(paste(outpath,'/qual_plot_r.png',sep=''), plot=pr)
        flog.info('Done.')
      }

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


      out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(f_trunclen,r_trunclen),
      maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, trimLeft=trim_l,
      compress=compress, multithread=TRUE)

      flog.info('Done.')

    }

    #COMMON
    flog.info('Learning error model...')
    errF <- learnErrors(filtFs, multithread=TRUE)
    errR <- learnErrors(filtRs, multithread=TRUE)
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


    flog.info('Dereplicating fastq...')
    derepFs <- derepFastq(filtFs, verbose=TRUE)
    derepRs <- derepFastq(filtRs, verbose=TRUE)

    flog.info('Done.')
    flog.info('dada2...')
    dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool=dadapool, selfConsist=FALSE)
    stockFs <- sapply(dadaFs, getN)
    dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool=dadapool, selfConsist=FALSE)
    stockRs <- sapply(dadaRs, getN)
    flog.info('Done.')
    flog.info('Merging pairs...')
    mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=FALSE)
    flog.info('Done.')


    seqtab <- makeSequenceTable(mergers)

    flog.info('Removing chimeras...')
    seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
    flog.debug(sum(seqtab.nochim))
    flog.debug(sum(seqtab))

    flog.debug(sum(seqtab.nochim)/sum(seqtab))
    flog.info('Done.')


    track <- cbind.data.frame(out, stockFs, stockRs, sapply(mergers, getN), rowSums(seqtab.nochim))
    # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
    colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
    rownames(track) <- sample.names
    head(track)

    write.table(track, paste(outpath,"/read_tracking.csv",sep=''), sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)


    seqtab.export <- seqtab.nochim
    colnames(seqtab.export) <- sapply(colnames(seqtab.export), digest::digest, algo="md5")

    otu.table <- phyloseq::otu_table(t(seqtab.export), taxa_are_rows = TRUE)

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

    if(torrent_single == TRUE){
      out <- filterAndTrim(fwd = fnFs, filt = filtFs, maxN = 0, multithread = TRUE, verbose=TRUE, rm.phix = TRUE,
        , maxEE = 5 , minLen = 100, compress=compress, trimLeft=trim_l, trimRight=trim_r, orient.fwd = orient_torrent)
    }else{
      out <- filterAndTrim(fwd = fnFs, filt = filtFs, maxN = 0, multithread = TRUE, verbose=TRUE, rm.phix = TRUE,
        , maxEE = 5 , minLen = 100, compress=compress, trimLeft=trim_l, trimRight=trim_r )
    }
    row.names(out) = sample.names

    flog.info('Learning error model...')
    errF <- learnErrors(filtFs, multithread=TRUE)
    flog.info('Done.')

    if(plot){
      flog.info('Plotting 2 ...')
      pf2 <- plotErrors(errF, nominalQ=TRUE)
      ggsave(paste(outpath,'/err_plot_f.png',sep=''), plot=pf2)
      flog.info('Done.')
    }

    mergers <- vector("list", length(sample.names))
    names(mergers) <- sample.names;
    stockFs=NULL; stockRs=NULL
    getN <- function(x) sum(getUniques(x))


    flog.info('Dereplicating fastq...')
    if(compress==TRUE){
      filtFs <- sort(list.files(file.path(path, "filtered/"), pattern = ".fastq.gz$", full.names = TRUE))
    } else{
      filtFs <- sort(list.files(file.path(path, "filtered/"), pattern = ".fastq$", full.names = TRUE))
    }

    names(filtFs) <- sample.names
    derepFs <- derepFastq(filtFs, verbose=TRUE)
    flog.info('Done.')

    flog.info('dada2...')

    if(torrent_single == TRUE){
      dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool=dadapool, selfConsist=FALSE, HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)
    }
    else{
      dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool=dadapool, selfConsist=FALSE)
    }


    if(length(filtFs)<2){
      stockFs <- getUniques(dadaFs)
    }else{
      stockFs <- sapply(dadaFs, getN)
    }

    flog.info('Done.')

    seqtab <- makeSequenceTable(dadaFs)

    flog.info('Removing chimeras...')
    seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

    if(length(filtFs)<2){
      track <- c(out, sum(stockFs), rowSums(seqtab.nochim)[1])
      names(track) <- c("input", "filtered", "denoisedF", "nonchim")
    }else{
      nn = row.names(out)
      print(rownames(stockFs))
      track <- cbind.data.frame(out, stockFs[nn], rowSums(seqtab.nochim)[nn])
      colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
      head(track)
    }
    # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
    write.table(track, paste(outpath,"/read_tracking.csv",sep=''), sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

    seqtab.export <- seqtab.nochim
    colnames(seqtab.export) <- sapply(colnames(seqtab.export), digest::digest, algo="md5")

    otu.table <- phyloseq::otu_table(t(seqtab.export), taxa_are_rows = TRUE)
    colnames(otu.table) = sample.names

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
