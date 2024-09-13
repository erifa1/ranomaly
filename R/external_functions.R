#' qdaptools::mtabulate Tabulate Frequency Counts for Multiple Vectors 
#' 
#' Similar to \code{\link[base]{tabulate}} that works on multiple vectors.
#' 
#' @param vects A \code{\link[base]{vector}}, \code{\link[base]{list}}, or 
#' \code{\link[base]{data.frame}} of named/unnamed vectors.
#' @keywords tabulate frequency internal
#' @seealso \code{\link[base]{tabulate}}, \code{\link[qdapTools]{counts2list}}
#' @return Returns a \code{\link[base]{data.frame}} with columns equal to 
#' number of unique elements and the number of rows equal to the the original 
#' length of the \code{\link[base]{vector}}, \code{\link[base]{list}}, or 
#' \code{\link[base]{data.frame}} (length equals ncols in 
#' \code{\link[base]{data.frame}}).  If list of vectors is named 
#' these will be the rownames of the dataframe.
#' @author Joran Elias and Tyler Rinker <tyler.rinker@@gmail.com>.
#' @references \url{https://stackoverflow.com/a/9961324/1000343}
#' @examples 
#' mtabulate(list(w=letters[1:10], x=letters[1:5], z=letters))
#' mtabulate(list(mtcars$cyl[1:10]))
#' 
#' ## Dummy coding
#' mtabulate(mtcars$cyl[1:10])
#' mtabulate(CO2[, "Plant"])
#' 
#' dat <- data.frame(matrix(sample(c("A", "B"), 30, TRUE), ncol=3))
#' mtabulate(dat)
#' t(mtabulate(dat))
#' counts2list(mtabulate(dat))
mtabulate <- function(vects) { 
    lev <- sort(unique(unlist(vects)))
    dat <- do.call(rbind, lapply(vects, function(x, lev){ 
        tabulate(factor(x, levels = lev, ordered = TRUE),
        nbins = length(lev))}, lev = lev))
    colnames(dat) <- sort(lev) 
    data.frame(dat, check.names = FALSE)
}


## param drop.na logical.  If \code{TRUE} \code{NA} columns (elements that 
## contained just an \code{NA}) will be dropped.## @author akrun of StackOverflow and Tyler Rinker <tyler.rinker@@gmail.com>.
## references url{http://stackoverflow.com/a/32753233/1000343}
## mtabulate <- function(vects, drop.na = TRUE) {
## 
##     x <- y <- . <- NULL
##     vects <- as.list(vects)
##     if (is.null(names(vects))) names(vects) <- seq_along(vects)
##     dat <- data.table::data.table(
##         x = names(vects),
##         y = vects, 
##         stringsAsFactors = FALSE
##     )
##     dat$y <- relist(unlist(dat$y), skeleton=dat$y)
##     data.table::setDT(dat)
##     dat <- dat[, .(y = unlist(y)), by = x]
##     out <- suppressMessages(data.table::dcast(dat, x ~ y, fun=length, drop=FALSE, fill=0))
##     out2 <- as.data.frame(out[, -1, with=FALSE])
##     rownames(out2) <- out[[1]]
##     if (isTRUE(drop.na)) out2[, "NA"] <- NULL
##     out2[names(vects), ]
## }
## }

#' gtools Order or Sort strings with embedded numbers so that the numbers are in the
#' correct order
#'
#' These functions sort or order character strings containing embedded numbers
#' so that the numbers are numerically sorted rather than sorted by character
#' value.  I.e. "Aspirin 50mg" will come before "Aspirin 100mg".  In addition,
#' case of character strings is ignored so that "a", will come before "B" and
#' "C".
#'
#' I often have character vectors (e.g. factor labels), such as compound and
#' dose, that contain both text and numeric data.  This function is useful for
#' sorting these character vectors into a logical order.
#'
#' It does so by splitting each character vector into a sequence of character
#' and numeric sections, and then sorting along these sections, with numbers
#' being sorted by numeric value (e.g. "50" comes before "100"), followed by
#' characters strings sorted by character value (e.g. "A" comes before "B")
#' \emph{ignoring case} (e.g. 'A' has the same sort order as 'a').
#'
#' By default, sort order is ascending, empty strings are sorted to the front,
#' and \code{NA} values to the end.  Setting \code{descending=TRUE} changes the
#' sort order to descending and reverses the meanings of \code{na.last} and
#' \code{blank.last}.
#'
#' Parsing looks for decimal numbers unless \code{numeric.type="roman"}, in
#' which parsing looks for roman numerals, with character case specified by
#' \code{roman.case}.
#'
#' @aliases mixedsort mixedorder
#' @param x Vector to be sorted.
#' @param decreasing logical.  Should the sort be increasing or decreasing?
#' Note that \code{descending=TRUE} reverses the meanings of \code{na.last} and
#' \code{blanks.last}.
#' @param na.last for controlling the treatment of \code{NA} values.  If
#' \code{TRUE}, missing values in the data are put last; if \code{FALSE}, they
#' are put first; if \code{NA}, they are removed.
#' @param blank.last for controlling the treatment of blank values.  If
#' \code{TRUE}, blank values in the data are put last; if \code{FALSE}, they
#' are put first; if \code{NA}, they are removed.
#' @param numeric.type either "decimal" (default) or "roman".  Are numeric
#' values represented as decimal numbers (\code{numeric.type="decimal"}) or as
#' Roman numerals (\code{numeric.type="roman"})?
#' @param roman.case one of "upper", "lower", or "both".  Are roman numerals
#' represented using only capital letters ('IX') or lower-case letters ('ix')
#' or both?
#' @param scientific logical. Should exponential notation be allowed for numeric values.
#' @return \code{mixedorder} returns a vector giving the sort order of the
#' input elements. \code{mixedsort} returns the sorted vector.
#' @author Gregory R. Warnes \email{greg@@warnes.net}
#' @seealso \code{\link[base]{sort}}, \code{\link[base]{order}}
#' @keywords univar manip
#' @examples
#'
#' ## compound & dose labels
#' Treatment <- c(
#'   "Control", "Aspirin 10mg/day", "Aspirin 50mg/day",
#'   "Aspirin 100mg/day", "Acetomycin 100mg/day",
#'   "Acetomycin 1000mg/day"
#' )
#'
#' ## ordinary sort puts the dosages in the wrong order
#' sort(Treatment)
#'
#' ## but mixedsort does the 'right' thing
#' mixedsort(Treatment)
#'
#' ## Here is a more complex example
#' x <- rev(c(
#'   "AA 0.50 ml", "AA 1.5 ml", "AA 500 ml", "AA 1500 ml",
#'   "EXP 1", "AA 1e3 ml", "A A A", "1 2 3 A", "NA", NA, "1e2",
#'   "", "-", "1A", "1 A", "100", "100A", "Inf"
#' ))
#'
#' mixedorder(x)
#'
#' mixedsort(x) # Notice that plain numbers, including 'Inf' show up
#' # before strings, NAs at the end, and blanks at the
#' # beginning .
#'
#'
#' mixedsort(x, na.last = TRUE) # default
#' mixedsort(x, na.last = FALSE) # push NAs to the front
#'
#'
#' mixedsort(x, blank.last = FALSE) # default
#' mixedsort(x, blank.last = TRUE) # push blanks to the end
#'
#' mixedsort(x, decreasing = FALSE) # default
#' mixedsort(x, decreasing = TRUE) # reverse sort order
#'
#' ## Roman numerals
#' chapters <- c(
#'   "V. Non Sequiturs", "II. More Nonsense",
#'   "I. Nonsense", "IV. Nonesensical Citations",
#'   "III. Utter Nonsense"
#' )
#' mixedsort(chapters, numeric.type = "roman")
#'
#' ## Lower-case Roman numerals
#' vals <- c(
#'   "xix", "xii", "mcv", "iii", "iv", "dcclxxii", "cdxcii",
#'   "dcxcviii", "dcvi", "cci"
#' )
#' (ordered <- mixedsort(vals, numeric.type = "roman", roman.case = "lower"))
#' roman2int(ordered)
#'
#' ## Control scientific notation for number matching:
#' vals <- c("3E1", "2E3", "4e0")
#'
#' mixedsort(vals) # With scientfic notation
#' mixedsort(vals, scientific = FALSE) # Without scientfic notation
#' @keywords internal
mixedsort <- function(x,
                      decreasing = FALSE,
                      na.last = TRUE,
                      blank.last = FALSE,
                      numeric.type = c("decimal", "roman"),
                      roman.case = c("upper", "lower", "both"),
                      scientific = TRUE) {
  ord <- mixedorder(x,
    decreasing = decreasing,
    na.last = na.last,
    blank.last = blank.last,
    numeric.type = numeric.type,
    roman.case = roman.case,
    scientific = scientific
  )
  x[ord]
}

#' @rdname mixedsort
#' @keywords internal
mixedorder <- function(x,
                       decreasing = FALSE,
                       na.last = TRUE,
                       blank.last = FALSE,
                       numeric.type = c("decimal", "roman"),
                       roman.case = c("upper", "lower", "both"),
                       scientific = TRUE) {
  # - Split each each character string into an vector of strings and
  #   numbers
  # - Separately rank numbers and strings
  # - Combine orders so that strings follow numbers

  numeric.type <- match.arg(numeric.type)
  roman.case <- match.arg(roman.case)

  if (length(x) < 1) {
    return(NULL)
  } else if (length(x) == 1) {
    return(1)
  }

  if (!is.character(x)) {
    return(order(x, decreasing = decreasing, na.last = na.last))
  }

  delim <- "\\$\\@\\$"

  if (numeric.type == "decimal") {
    if (scientific) {
      regex <- "((?:(?i)(?:[-+]?)(?:(?=[.]?[0123456789])(?:[0123456789]*)(?:(?:[.])(?:[0123456789]{0,}))?)(?:(?:[eE])(?:(?:[-+]?)(?:[0123456789]+))|)))"
    } # uses PERL syntax
    else {
      regex <- "((?:(?i)(?:[-+]?)(?:(?=[.]?[0123456789])(?:[0123456789]*)(?:(?:[.])(?:[0123456789]{0,}))?)))"
    } # uses PERL syntax

    numeric <- function(x) as.numeric(x)
  }
  else if (numeric.type == "roman") {
    regex <- switch(roman.case,
      "both"  = "([IVXCLDMivxcldm]+)",
      "upper" = "([IVXCLDM]+)",
      "lower" = "([ivxcldm]+)"
    )
    numeric <- function(x) roman2int(x)
  }
  else {
    stop("Unknown value for numeric.type: ", numeric.type)
  }

  nonnumeric <- function(x) {
    ifelse(is.na(numeric(x)), toupper(x), NA)
  }

  x <- as.character(x)

  which.nas <- which(is.na(x))
  which.blanks <- which(x == "")

  ####
  # - Convert each character string into an vector containing single
  #   character and  numeric values.
  ####

  # find and mark numbers in the form of +1.23e+45.67
  delimited <- gsub(regex,
    paste(delim, "\\1", delim, sep = ""),
    x,
    perl = TRUE
  )

  # separate out numbers
  step1 <- strsplit(delimited, delim)

  # remove empty elements
  step1 <- lapply(step1, function(x) x[x > ""])

  # create numeric version of data
  suppressWarnings(step1.numeric <- lapply(step1, numeric))

  # create non-numeric version of data
  suppressWarnings(step1.character <- lapply(step1, nonnumeric))

  # now transpose so that 1st vector contains 1st element from each
  # original string
  maxelem <- max(sapply(step1, length))

  step1.numeric.t <- lapply(
    1:maxelem,
    function(i) {
      sapply(
        step1.numeric,
        function(x) x[i]
      )
    }
  )

  step1.character.t <- lapply(
    1:maxelem,
    function(i) {
      sapply(
        step1.character,
        function(x) x[i]
      )
    }
  )

  # now order them
  rank.numeric <- sapply(step1.numeric.t, rank)
  rank.character <- sapply(
    step1.character.t,
    function(x) as.numeric(factor(x))
  )

  # and merge
  rank.numeric[!is.na(rank.character)] <- 0 # mask off string values

  rank.character <- t(
    t(rank.character) +
      apply(matrix(rank.numeric), 2, max, na.rm = TRUE)
  )

  rank.overall <- ifelse(is.na(rank.character), rank.numeric, rank.character)

  order.frame <- as.data.frame(rank.overall)
  if (length(which.nas) > 0) {
    if (is.na(na.last)) {
      order.frame[which.nas, ] <- NA
    } else if (na.last) {
      order.frame[which.nas, ] <- Inf
    } else {
      order.frame[which.nas, ] <- -Inf
    }
  }

  if (length(which.blanks) > 0) {
    if (is.na(blank.last)) {
      order.frame[which.blanks, ] <- NA
    } else if (blank.last) {
      order.frame[which.blanks, ] <- 1e99
    } else {
      order.frame[which.blanks, ] <- -1e99
    }
  }

  order.frame <- as.list(order.frame)
  order.frame$decreasing <- decreasing
  order.frame$na.last <- NA

  retval <- do.call("order", order.frame)

  return(retval)
}




#' ranacapa: Make a rarefaction curve using ggplot2
#' @param physeq_object A phyloseq class object, from which abundance data are extracted
#' @param step Step Size for sample size in rarefaction curves
#' @param label Default `NULL`. Character string. The name of the variable to map to text labels on the plot. Similar to color option but for plotting text.
#' @param color Default `NULL`. Character string. The name of the variable to map to the colors in the plot. This can be a sample variables among the set returned by sample_variables(physeq_object) or taxonomic rank, among the set returned by rank_names(physeq_object)
#' @param plot default `TRUE`. Logical. Should the graph be plotted
#' @param parallel default `FALSE`. Logical. Should rarefaction be parallelized
#' @param se default `TRUE`. Logical. Should standard errors be calculated.
#' @examples
#' \dontrun{
#' ggrare(physeq_object, step = 20, se = TRUE)
#' }
#' @keywords internal

ggrare <- function (physeq_object, step = 10, label = NULL, color = NULL, 
    plot = TRUE, parallel = FALSE, se = TRUE) 
{
    x <- methods::as(phyloseq::otu_table(physeq_object), "matrix")
    if (phyloseq::taxa_are_rows(physeq_object)) {
        x <- t(x)
    }
    tot <- rowSums(x)
    S <- rowSums(x > 0)
    nr <- nrow(x)
    rarefun <- function(i) {
        cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
        n <- seq(1, tot[i], by = step)
        if (n[length(n)] != tot[i]) {
            n <- c(n, tot[i])
        }
        y <- vegan::rarefy(x[i, , drop = FALSE], n, se = se)
        if (nrow(y) != 1) {
            rownames(y) <- c(".S", ".se")
            return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
        }
        else {
            return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
        }
    }
    if (parallel) {
        out <- parallel::mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
    }
    else {
        out <- lapply(seq_len(nr), rarefun)
    }
    df <- do.call(rbind, out)
    if (!is.null(phyloseq::sample_data(physeq_object, FALSE))) {
        sdf <- methods::as(phyloseq::sample_data(physeq_object), 
            "data.frame")
        sdf$Sample <- rownames(sdf)
        data <- merge(df, sdf, by = "Sample")
        labels <- data.frame(x = tot, y = S, Sample = rownames(x))
        labels <- merge(labels, sdf, by = "Sample")
    }
    if (length(color) > 1) {
        data$color <- color
        names(data)[names(data) == "color"] <- deparse(substitute(color))
        color <- deparse(substitute(color))
    }
    if (length(label) > 1) {
        labels$label <- label
        names(labels)[names(labels) == "label"] <- deparse(substitute(label))
        label <- deparse(substitute(label))
    }
    p <- ggplot2::ggplot(data = data, ggplot2::aes_string(x = "Size", 
        y = ".S", group = "Sample", color = color))
    p <- p + ggplot2::labs(x = "Sequence Sample Size", y = "Species Richness")
    if (!is.null(label)) {
        p <- p + ggplot2::geom_text(data = labels, ggplot2::aes_string(x = "x", 
            y = "y", label = label, color = color), size = 4, 
            hjust = 0)
    }
    p <- p + ggplot2::geom_line()
    if (se) {
        p <- p + ggplot2::geom_ribbon(ggplot2::aes_string(ymin = ".S - .se", 
            ymax = ".S + .se", color = NULL, fill = color), alpha = 0.2)
    }
    if (plot) {
        plot(p)
    }
    invisible(p)
}



#'@title Pairwise multilevel comparison using adonis (from https://github.com/pmartinezarbizu/pairwiseAdonis)
#'
#'@description This is a wrapper function for multilevel pairwise comparison
#' using adonis() from package 'vegan'. The function returns adjusted p-values using p.adjust().
#'
#'@param x Data frame (the community table), or "dist" object (user-supplied distance matrix).
#'
#'@param factors Vector (a column or vector with the levels to be compared pairwise).
#'
#'@param sim.function Function used to calculate the similarity matrix,
#' one of 'daisy' or 'vegdist' default is 'vegdist'. Ignored if x is a distance matrix.
#'
#'@param sim.method Similarity method from daisy or vegdist, default is 'bray'. Ignored if x is a distance matrix.
#'
#'@param p.adjust.m The p.value correction method, one of the methods supported by p.adjust(),
#' default is 'bonferroni'.
#'
#'@param reduce String. Restrict comparison to pairs including these factors. If more than one factor, separate by pipes like  reduce = 'setosa|versicolor'
#' 
#' @param perm Number of permutations, default is 999.
#'
#'@return Table with the pairwise factors, Df, SumsOfSqs, F-values, R^2, p.value and adjusted p.value.
#'
#'@author Pedro Martinez Arbizu & Sylvain Monteux
#'
#'@examples
#' data(iris)
#' pairwise.adonis(iris[,1:4],iris$Species)
#'
#'
#'#similarity euclidean from vegdist and holm correction
#' pairwise.adonis(x=iris[,1:4],factors=iris$Species,sim.function='vegdist',
#' sim.method='euclidian',p.adjust.m='holm')
#' 
#'#identical example using a distance matrix as an input
#' dist_matrix=vegan::vegdist(iris[,1:4],method="euclidean")
#' pairwise.adonis(dist_matrix,factors=iris$Species,
#' p.adjust.m='holm')
#'
#'#similarity manhattan from daisy and bonferroni correction
#' pairwise.adonis(x=iris[,1:4],factors=iris$Species,sim.function='daisy',
#' sim.method='manhattan',p.adjust.m='bonferroni')
#' 
#'#Restrict comparison to only some factors
#'pairwise.adonis(iris[,1:4],iris$Species, reduce='setosa')
#'
#'#for more than one factor separate by pipes
#'pairwise.adonis(iris[,1:4],iris$Species, reduce='setosa|versicolor')
#' 
#'
#'@export pairwise.adonis
#'@importFrom stats p.adjust
#'@importFrom utils combn
#'@importFrom vegan adonis2 vegdist
#'@importFrom cluster daisy



pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni',reduce=NULL,perm=999)
{
  
  co <- combn(unique(as.character(factors)),2)
  pairs <- c()
  Df <- c()
  SumsOfSqs <- c()
  F.Model <- c()
  R2 <- c()
  p.value <- c()
  
  
  for(elem in 1:ncol(co)){
    if(inherits(x, 'dist')){
      x1=as.matrix(x)[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),
                      factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))]
      }
  
    else  (
      if (sim.function == 'daisy'){
            x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
        } 
      else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    )
    
    x2 = data.frame(Fac = factors[factors %in% c(co[1,elem],co[2,elem])])
    
    ad <- vegan::adonis2(x1 ~ Fac, data = x2,
                 permutations = perm);
    pairs <- c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    Df <- c(Df,ad$Df[1])
	SumsOfSqs <- c(SumsOfSqs,ad$SumOfSqs[1])
	F.Model <- c(F.Model,ad$F[1]);
    R2 <- c(R2,ad$R2[1]);
    p.value <- c(p.value,ad$`Pr(>F)`[1])
  }
  p.adjusted <- p.adjust(p.value,method=p.adjust.m)
  
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  pairw.res <- data.frame(pairs,Df,SumsOfSqs,F.Model,R2,p.value,p.adjusted,sig)
  
  if(!is.null(reduce)){
    pairw.res <- subset (pairw.res, grepl(reduce,pairs))
    pairw.res$p.adjusted <- p.adjust(pairw.res$p.value,method=p.adjust.m)
    
    sig = c(rep('',length(pairw.res$p.adjusted)))
 	sig[pairw.res$p.adjusted <= 0.1] <-'.'
	sig[pairw.res$p.adjusted <= 0.05] <-'*'
	sig[pairw.res$p.adjusted <= 0.01] <-'**'
	sig[pairw.res$p.adjusted <= 0.001] <-'***'
    pairw.res <- data.frame(pairw.res[,1:7],sig)
  }
  class(pairw.res) <- c("pwadonis", "data.frame")
  return(pairw.res)
} 


### Method summary
summary.pwadonis = function(object, ...) {
  cat("Result of pairwise.adonis:\n")
  cat("\n")
  print(object, ...)
  cat("\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
}
## end of method summary