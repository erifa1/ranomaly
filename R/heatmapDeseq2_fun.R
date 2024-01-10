#' heatmap for deseq2_fun output
#'
#' @param desq output of deseq2_fun
#' @param phys the phyloseq object used for deseq2 analyse
#' @param var Factor to test.
#' @param workingRank Taxonomy rank used in deseq2_fun function. 
#' @param output Output directory
#' @return Return a list object with heatmap plots.
#' @import futile.logger
#' @import stringr
#' @import phyloseq
#' @import pheatmap
#' @export
heatmapDeseq2_fun <- function(desq, phys, var, workingRank, output = "./deseq/heatmap/") {
  dir.create(output, recursive = TRUE)
  pheat = list()
  for (i in names(desq)) {
    tableDESeqSig <- desq[[i]][["table"]] |> subset(padj < 0.05) |> na.omit()

    # remove unused taxa
    psTaxaRelSig <- prune_taxa(row.names(tableDESeqSig),
                               transform_sample_counts(phys, function(x) x/sum(x)))

    # remove unused sample
    keepSamp <<- stringr::str_split(i, "_vs_", simplify = TRUE)

    psTaxaRelSig2 <- paste("psTaxaRelSigSubset <<- subset_samples(psTaxaRelSig, ",
                           var, " %in% keepSamp)", sep = "")
    eval(parse(text = psTaxaRelSig2))
    psTaxaRelSigSubset <- subset_taxa(psTaxaRelSigSubset,
                                      taxa_sums(psTaxaRelSigSubset) != 0)
    matrix <- otu_table(psTaxaRelSigSubset) |>
      data.frame() |>
      as.matrix()
    if (!taxa_are_rows(psTaxaRelSigSubset)) {
      colnames(matrix) <- tax_table(psTaxaRelSigSubset)[, workingRank] |>
        as.character()
    } else {
      rownames(matrix) <- tax_table(psTaxaRelSigSubset)[, workingRank] |>
        as.character()
    }

    metadataSub <- sample_data(psTaxaRelSigSubset) |>
      data.frame()
    # print(metadataSub)

    if (length(var) == 1) {
      annotationCol <- data.frame(
        var1 = metadataSub[, var[1]] |>
          as.vector() |>
          as.factor(),
        check.names = FALSE)
    } else {
      annotationCol <- data.frame(
        var1 = metadataSub[, var[1]] |>
          as.vector() |>
          as.factor(),
        var2 = metadataSub[, var[2]] |>
          as.vector() |>
          as.factor(),
        check.names = FALSE)
    }

    rownames(annotationCol) <- rownames(metadataSub)
    colnames(annotationCol) <- var

    annotationRow <- data.frame(
      Phylum = tax_table(psTaxaRelSigSubset)[, "Phylum"] |>
        as.factor()
    )
    if (taxa_are_rows(psTaxaRelSigSubset)) {
      rownames(annotationRow) <- rownames(matrix)
    } else {
      rownames(annotationRow) <- colnames(matrix)
    }
    title <- i 

    pheat[[i]] <- pheatmap::pheatmap(matrix, scale = "row",
                           annotationCol = annotationCol,
                           annotationRow = annotationRow,
                           main = title)

    png(glue::glue("{output}/pheatmap-{i}.png"), , width = 1200, height = 800, res = 150)
    print(pheat[[i]])
    dev.off()

  }
  return(pheat)
}
