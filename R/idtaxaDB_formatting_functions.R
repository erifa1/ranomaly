#' Fill taxonomy table (idtaxa training functions)
#'
#'
#'
#' @param taxtable Taxonomy table in tabulated format.
#' @param prefix TRUE if there are prefix like c("k__","p__","c__","o__","f__","g__","s__")
#'
#' @return Return the same taxonomy without empty field, last known ranks are informed.
#'
#' @import futile.logger
#' @import DECIPHER
#' @export

fill_tax_fun <- function(taxtable = taxtable, prefix = TRUE){

  fill_tax_table = function(x){
    RANKS = c("domain","phylum","class","order","family","genus","species")
    if(prefix){
      PREFIX = c("k__","p__","c__","o__","f__","g__","s__")
    }else{
      PREFIX = c("","","","","","","")
    }

    TAX = x

    if(length(na.omit(TAX))==0){
      fTAX = paste(PREFIX,"unassigned_", RANKS, sep="")
      return(fTAX)
    }

    for( i in 1:7){
      if(is.na(TAX[i])){
        if(grepl(PREFIX[i-1], TAX[i-1], ignore.case = FALSE, perl = FALSE, fixed = TRUE)){
          TAX[i-1] <- sub(PREFIX[i-1],'',TAX[i-1])
        }
        x[i] = paste(PREFIX[i],TAX[i-1],'_', RANKS[i], sep="")
        TAX[i] = paste(PREFIX[i],TAX[i-1], sep="")
      }
      # print(x)
    }
    return(x)
  }


  filltable <- tt2 <- as.data.frame( t(apply(taxtable, 1, fill_tax_table)) ,stringsAsFactors = FALSE )
  names(filltable) = c("Domain","Phylum","Class","Order","Family","Genus","Species")
  rownames(filltable) =  rownames(taxtable)

  return(filltable)

}


#' Check taxonomy table (idtaxa training functions)
#'
#' @param taxtable Output from fill_tax_table function.
#' @param output output path
#' @param verbose verbose mode
#' @param returnval Boolean to return values in console or not.

#' @return Return the same taxonomy without incongruence like multiple ancestors, by replacing problematic taxonomy by the most common.
#'
#' @export



check_tax_fun <- function(taxtable = taxtable, output = NULL, verbose=3, returnval = TRUE){

  print("Check taxonomy consistency...")
  # Check for multiple ancestors at each rank, choose first occurence for each problematic taxon
  for(rank in 6:2){
    if(verbose==3){flog.info(paste(colnames(taxtable)[rank],".",sep=""))}

    stockLres = 100; nloop=1
    while(max(stockLres)>1){
      if(verbose==3){flog.info(paste("Loop",nloop,".",sep=""))}
      stockLres = NULL
      allTax = apply(taxtable[,1:rank] ,1, function(x){paste(x, collapse = ";")})
      uniqTax <- unique(allTax)

      #Test if unique Genus has unique taxonomy.
      for (i in unique(taxtable[,rank])){
        # if(verbose==3){flog.info(paste(i,".",sep=""))}
        res=grep( paste(';',i,'$', sep="" ) , uniqTax)
        stockLres = c(stockLres,length(res))
        #If multiple tax for one taxon...
        if(length(res)>1){
          cat("\n");print(i);
                      print(uniqTax[res]);

          tax2 <- apply(taxtable[taxtable[,rank]==i,1:rank], 1, function(x){paste(x, collapse = ";")})
          uniqTax2 <- table(tax2)
          ftax <- names(uniqTax2[order(uniqTax2,decreasing=TRUE)])[1]
          ftax <- unlist(strsplit(ftax,";"))
          cat( glue::glue( "CORRECTED :
                      {paste(ftax, collapse = ';')}" )
          ); cat("\n")
          #Change taxonomy with final ftax. the most common in taxtable
          for(j in row.names(taxtable[taxtable[,rank]==i,])){
            taxtable[j,] = c(ftax, taxtable[j,(rank+1):ncol(taxtable)])
          }
        }
      }
      nloop = nloop + 1
    }
  }
  if(!is.null(output)){
    print("OUTPUT")
    write.table(taxtable, output, sep = "\t", col.names=NA)
  }
  flog.info("Done.")
  if(returnval){
    return(taxtable)
  }

}


#' Taxid processing (idtaxa training functions)
#'
#' @param taxtable Output from check_tax_fun function.
#' @param output output path

#' @return Return the taxid object.
#'
#' @export



taxid_fun <- function(taxtable = taxtable, output = "./taxid.txt"){

  taxa <- setNames(c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
                   c("k__", "p__", "c__", "o__", "f__", "g__", "s__"))
  # ranks <- setNames( split(check1, seq(nrow(check1))) , rownames(check1) )
  ranks <- split(as.matrix(taxtable), seq(nrow(taxtable)))

  ## PARSE RANKS taxid
  count <- 1L
  groups <- "Root"
  index <- -1L
  level <- 0L
  rank <- "rootrank"
  pBar <- txtProgressBar(style=3)
  for (i in seq_along(ranks)) {
    for (j in seq_along(ranks[[i]])) {
      #print(c(i,j))
      rank_level <- taxa[substring(ranks[[i]][j], 1, 3)]
      #print(rank_level)
      group <- ranks[[i]][j] 			#substring(, 4)#
      #print(group)
      w <- which(groups==group & rank==rank_level)  #Verifie si le groupe est déjà présent
      if (length(w) > 0) {  #si oui check le parent et next
        parent <- match(ranks[[i]][j - 1], #substring(, 4)
                        groups)
        if (j==1 || any((parent - 1L)==index[w]))
          next # already included
      }

      count <- count + 1L
      groups <- c(groups, group)
      if (j==1) { #condition pour les index level rank
        index <- c(index, 0)
      } else {
        parent <- match(ranks[[i]][j - 1], #substring(, 4)
                        groups)
        index <- c(index,
                   parent - 1L)
      }
      level <- c(level, j)
      rank <- c(rank, taxa[j])
    }

    setTxtProgressBar(pBar, i/length(ranks))
  }
  groups <- gsub("^[ ]+", "", groups)
  groups <- gsub("[ ]+$", "", groups)

  if(!is.null(output)){
    writeLines(paste(0:(length(index) - 1L), groups, index, level, rank, sep="*"),
               con=file(output))
  }

  taxid = cbind.data.frame(0:(length(index) - 1L), groups, index, level, rank)
  names(taxid) <- c("Index", "Name", "Parent", "Level", "Rank")
  print(head(taxid))

  return(taxid)
}


#' IDTAXA db training (idtaxa training functions)
#'
#' @param taxtable Output from check_tax_fun function.
#' @param taxid Output from taxid_fun function
#' @param seqs FASTA files of reference sequence
#' @param outputDIR output directory path
#' @param outputDBname db name
#' @param returnval Boolean to return values in console or not.

#' @return Save rdata files with the formatted idtaxa reference database.
#'
#' @export


idtaxa_traindb <- function(taxtable = taxtable, taxid = taxid, seqs = "", prunedb = NULL, outputDIR = "./",
                           outputDBname = "newDB.rdata", returnval = FALSE){

  dna <- readDNAStringSet(seqs)

  if( any(rownames(check1) != names(dna)) ){
    stop("sequence IDS and taxonomy IDS do not exactly match... Check IDS and order.")
  }

  cat("\tTaxonomy\n")
  taxonomy <- paste( "Root", apply(taxtable, 1, paste, collapse = "; "), sep="; ")


  if(!is.null(prunedb)){
    cat("\n\tPruningDB\n")
    #Prune DB
    groups <- taxonomy
    groupCounts <- table(groups)
    u_groups <- names(groupCounts) # unique groups
    #length(u_groups) # number of groups

    #Pruning DB
    maxGroupSize <- prunedb # max sequences per label (>= 1)
    remove <- logical(length(dna))
    for (i in which(groupCounts > maxGroupSize)) {
      index <- which(groups==u_groups[i])
      keep <- sample(length(index), maxGroupSize)
      remove[index[-keep]] <- TRUE
    }
    # number of sequences eliminated
    cat("\n",paste("Number of sequences eliminated:",sum(remove)),"\n")

    prune_dna = dna[!remove]
    prune_tax = taxonomy[!remove]

    cat("\n",paste("Remaining groups:",length(prune_tax)),"\n")
    cat("\n",paste("Remaining sequences:",length(prune_dna)),"\n")

    cat("\nTrain classifier\n")
    trainingSet <- LearnTaxa(prune_dna, prune_tax)
    trainingSet

  }else{
    cat("\tNo Pruning\n")
    cat("\nTrain classifier\n")
    # train the classifier
    trainingSet <- LearnTaxa(dna, taxonomy, rank=taxid)
    trainingSet
  }

  if(!dir.exists(outputDIR)){
    dir.create(outputDIR, recursive = TRUE)
  }
  cat(glue::glue("\tSaving {paste(outputDIR, outputDBname,sep='/')} ...\n"))
      save(trainingSet, file = paste(outputDIR, outputDBname,sep="/"))

}





