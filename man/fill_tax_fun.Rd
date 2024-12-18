% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/idtaxaDB_formatting_functions.R
\name{fill_tax_fun}
\alias{fill_tax_fun}
\title{Fill taxonomy table (idtaxa training functions)}
\usage{
fill_tax_fun(
  taxtable = taxtable,
  prefix = c("k__", "p__", "c__", "o__", "f__", "g__", "s__"),
  ranks_names = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
)
}
\arguments{
\item{taxtable}{Taxonomy table in tabulated format.}

\item{prefix}{Vector of prefixes to use like c("k__","p__","c__","o__","f__","g__","s__"), NULL if no prefix.}

\item{ranks_names}{Taxonomy ranks names}
}
\value{
Return the same taxonomy without empty field, last known ranks are informed.
}
\description{
Fill taxonomy table (idtaxa training functions)
}
