#' core_soft_fun
#'
#' @description Define core microbiome, soft microbiome and transitory microbiome.
#'
#' @param data A phyloseq object
#' @param fact Factor to test (must be in sample_variables(data))
#' @param group Choose which level of the factor, if NULL generate a list for each level.
#' @param freq frequence threshold of microbiome::core_members function
#' @param prev prevalence threshold of microbiome::core_members function
#' @param rank Taxonomy rank.
#'
#' @importFrom microbiome transform core_members
#'
#' @export


core_soft_fun <- function(data = NULL, fact = NULL, group = NULL, freq = 0.001, prev = 0.5, rank = "ASV"){
  if(is.null(data)|is.null(fact)){stop("Require a phyloseq object, factor and group...")}
  if(! any( sample_variables(data) == fact ) ){stop(glue::glue("{fact} is not in sample_variables(phyloseq_object)..."))}

  if(rank == "ASV"){
    glomps <- data
  }else{
    glomps <- tax_glom(data, taxrank = rank)
    taxa_names(glomps) = glomps@tax_table@.Data[,rank]
  }


  fun2 <- function(glomps = glomps, fact = fact, group = group){
    LL = list()
    fun <- glue::glue("sub_phy <- subset_samples(glomps, {fact} == '{group}')")
    eval(parse(text=fun))
    sub_phy <- prune_taxa(taxa_sums(sub_phy) > 0, sub_phy)
    sub_phy.rel <- microbiome::transform(sub_phy, "compositional")

    #core
    LL$core <- core.taxa.standard <- core_members(sub_phy.rel, detection = freq, prevalence = prev)
    #soft in at least 1 sample
    core.plus <- core_members(sub_phy.rel, detection = freq, prevalence = 1/nsamples(data))
    LL$soft = setdiff(core.plus, core.taxa.standard)
    # transit
    LL$transit = setdiff(taxa_names(sub_phy.rel), core.plus)

    return(LL)
  }

  if(!is.null(group)){
    return(fun2(glomps = glomps, fact = fact , group = group))

  }else{
    LL2 = list()
    for(group in levels(as.data.frame(as.matrix(sample_data(data)))[,fact])){
      nlist = gsub(" ", "_" ,group)
      fun <- glue::glue("LL2${nlist} <- fun2(glomps = glomps, fact = fact, group = '{group}')")
      eval(parse(text=fun))
    }
    return(LL2)
  }

}


#' List generator for venndiagram
#'
#' @description Output list to generate venn diagram with VENNFUN function.
#'
#' @param inlist Output of core_soft_fun() function
#' @param part which part of microbiome to output (core, soft or transit)
#'
#' @export

list_venn_fun <- function(inlist = NULL, part="core"){
  if(! any(part == c("core", "soft", "transit"))){stop("Choose 'core', 'soft' or 'transit' for 'part' argument")}

  outlist = list()
  for( i in names(inlist)){
    fun = glue::glue("outlist${i} = inlist${i}${part}")
    eval(parse(text=fun))
  }
  return(outlist)
}
