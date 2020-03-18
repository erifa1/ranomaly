#' DADA2 process
#'
#' Processing DADA2 algorithm on raw sequences, return raw otu table in phyloseq object.
#'
#' @param amplicon Choose amplipcon "16S" or "ITS"
#'
#' @return Return raw otu table in phyloseq object.
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

dada2_fun <- function(amplicon = "16S", path = "", outpath = "./dada2_out/", f_trunclen = 240, r_trunclen = 240,
                      f_primer = "GCATCGATGAAGAACGCAGC", r_primer = "TCCTCCGCTTWTTGWTWTGC", plot = FALSE, compress = FALSE, verbose = 1){

  if(verbose == 3){
    invisible(flog.threshold(DEBUG))
  } else {
    invisible(flog.threshold(INFO))
  }




  flog.info('Done.')
  wd <- getwd()



  flog.info("Creating directory.")
  if(!dir.exists(outpath)){
    dir.create(outpath)
  }
  flog.info('Done.')
  ## Reads path
  # path <- opt$path
  #list.files(path)

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
  flog.debug(fnFs)
  flog.debug(fnRs)
  flog.info('Done.')

  if(amplicon=="ITS"){
    flog.info('DADA2 ITS')

    # Primers sequences
    FWD <- f_primer
    REV <- r_primer

    # Generate all primers orientations
    allOrients <- function(primer) {
      # Create all orientations of the input sequence
      require(Biostrings)
      dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
      orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna),
                   RevComp = reverseComplement(dna))
      return(sapply(orients, toString))  # Convert back to character vector
    }
    FWD.orients <- allOrients(FWD)
    REV.orients <- allOrients(REV)
    #FWD.orients

    # Remove Ns from reads
    fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
    fnRs.filtN <- file.path(path, "filtN", basename(fnRs))

    print(fnFs.filtN)

    flog.info('filterAndTrim...')
    if(! dir.exists(paste(path,'/filtN',sep=''))){
      filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE, verbose=TRUE, rm.phix = TRUE, compress=compress)
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
    get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
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

    sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

    if(compress){
      filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
      filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
    }else{
      filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
      filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))
    }


    names(filtFs) <- sample.names
    names(filtRs) <- sample.names

    # print(head(fnFs))
    # print(head(fnRs))
    # print(head(filtFs))
    # print(head(filtRs))
    #
    # print(head(f_trunclen))
    # print(head(r_trunclen))
    #
    # print(compress)

# -------------------------------------------------------------------------


    flog.info('Filtering reads...')


    out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(f_trunclen,r_trunclen),
                         maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, trimLeft=15,
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
    print("plotOK")
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
    dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool=FALSE, selfConsist=FALSE)
    stockFs <- c(stockFs, getN(dadaFs))
    dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool=FALSE, selfConsist=FALSE)
    stockRs <- c(stockRs,getN(dadaRs))
    flog.info('Done.')
    flog.info('Merging pairs...')
    merger <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=FALSE)
    flog.info('Done.')
    mergers[[sam]] <- merger
  }

  seqtab <- makeSequenceTable(mergers)

  flog.info('Removing chimeras...')
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
  flog.debug(sum(seqtab.nochim))
  flog.debug(sum(seqtab))

  flog.debug(sum(seqtab.nochim)/sum(seqtab))
  # flog.info(sum(seqtab.nochim)/sum(seqtab))
  flog.info('Done.')


  track <- cbind.data.frame(out, stockFs, stockRs, sapply(mergers, getN), rowSums(seqtab.nochim))
  # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  rownames(track) <- sample.names
  head(track)

  write.table(track, paste(outpath,"/read_tracking.csv",sep=''), sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

  # suppressMessages(library(digest))
  # suppressMessages(library(phyloseq))



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
  return(dada_res)


  flog.info('Saving R objects.')
  save(dada_res, file=paste(outpath,'/robjects.Rdata',sep=''))
  flog.info('Finish.')

}

