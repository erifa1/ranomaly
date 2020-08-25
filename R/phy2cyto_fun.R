#' Export files for cytoscape.
#'
#'
#' @param data a phyloseq object (output from decontam or generate_phyloseq)
#' @param output Output directory
#' @param column1 A metadata column (among sample_variables(data)).
#' @param repl Replicates column.
#' @param verbose Verbose level. (1: quiet, 3: verbal)
#'
#' @return Export files ready to use with Cytoscape
#'
#' @references https://cytoscape.org/
#'
#' @export



phy2cyto_fun <- function(data = data, output = "./cytoscape/", column1 = NULL, repl = NULL, verbose = 1){

  if(!dir.exists(output)){
    dir.create(output, recursive=TRUE)
  }


  data1 <- data
  otab = as.data.frame(otu_table(data1))
  sdat <- sdata <- as.data.frame(as.matrix(sample_data(data1)))
  fact <- levels(as.factor(sdata[,column1]))
  # repl <- "temps"


  if(is.null(repl)){
    flog.info('No replicates, performing links ...')
    #Test1 VUE globale sans tenir compte des réplicats
    # Interaction table
    sif_tab=NULL
    for(env in fact ){
      flog.debug(print(env))
      eval(parse(text=glue("tt <- subset_samples(data1, {column1} %in% '{env}')") ))   ###
      flog.debug(print(nsamples(tt)))

      #ASV present dans au moins 3 échantillons, autre filtre sur abondance?
      tt2 = prune_taxa( apply(otu_table(tt), 1, function(x){length(which(x>0))>3}), tt)
      flog.debug(print(ntaxa(tt2)))
      #Effectif des asv dans la table < 20
      otab2 <- as.data.frame(otu_table(tt2))
      otab2[otab2<20] = 0
      otu_table(tt2) = otu_table(otab2, taxa_are_rows=TRUE)
      tt2 <- prune_taxa(taxa_sums(tt2) > 0, tt2)

      tabf = cbind(taxa_names(tt2), rep(glue("type_{env}"), ntaxa(tt2)), rep(env, ntaxa(tt2)))
      sif_tab = rbind(sif_tab, tabf)
    }
    flog.info('Output ...')
    write.table(as.data.frame(sif_tab), glue("{output}/sif_tab1.tsv"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

    # Node table
    src = unique(sif_tab[,3])
    trgt = sif_tab[,1]
    freq = taxa_sums(data1) / sum(taxa_sums(data1))

    node_table=cbind.data.frame(c(src, trgt), c(rep("source", length(src)), rep("cible", length(trgt))), c(rep(NA, length(src)), freq[trgt]) )
    names(node_table)=c("names","attribute", "freq")
    node_table[node_table$attribute=="source","freq"] = 0.75*max(na.omit(node_table$freq))

    node_table$tax = "source"
    node_table[node_table$attribute == "cible","tax"] = tax_table(data1)[as.character(node_table[node_table$attribute == "cible",1]),"Phylum"]


    write.table(as.data.frame(node_table), glue("{output}/node_tab1.tsv"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

  }else{
    flog.info('Replicates, performing links ...')
    save(list = ls(all.names = TRUE), file = "debug_cyto1.rdata", envir = environment())
    sif_tab = NULL
    for(asv in taxa_names(data1)){
      # print(asv)
      notu = otab[asv,]
      nsamples = names(notu)[notu>20]  #nombre d'occurence de l'otu
      if(length(nsamples) == 1){flog.debug(print(asv))}

      sdat2 = sdat[nsamples,]
      rep1 = table(as.character(sdat[nsamples,repl]))
      reps = names(rep1)   #[which(rep1>1)]
      # if(max(rep1[which(rep1>1)]) == 1){print(asv)}

      # si aucun partage 1 seul replenvironnement
      if(length(rep1)==1){

        srcs <- as.character(unique(sdat2[,column1]))
        if(rep1==1){
          flog.debug(glue("{asv} exclusif 1 seul réplicat 1 env"))
          LINKS = c(asv, glue("type_{srcs}"), srcs)
        }else{
          # print(c(asv, "debug"))
          LINKS = cbind(rep(asv, length(srcs)), glue("type_{srcs}"), srcs)
        }
      } else{
        #sinon partage et création d'un lien pour chaque environnement source.
        LINKS = NULL
        for(rep in reps){
          sdat3 = sdat2[sdat2[,repl]==rep,]
          srcs = as.character(unique(sdat3[,column1]))
          # Lien seulement si l'asv présent dans même sample à 2 environnements ou plus
          if(length(srcs)>1){
            link = cbind(rep(asv, length(srcs)), glue("type_{srcs}"), srcs)
            # print(glue("{rep} {srcs} envs"))
            # print(link)
            LINKS = rbind(LINKS, link)
          }
        }
      }
      sif_tab = rbind(sif_tab, LINKS)

    }
    save(list = ls(all.names = TRUE), file = "debug_cyto.rdata", envir = environment())
    flog.info('Output ...')
    # LINK table output
    #Compte le nombre d'occurence de chaque lien (environnement - ASV)
    tt=apply(sif_tab, 1, paste, collapse="-")
    # Unique + counts de chaque lien
    sif_tabF = cbind(unique(sif_tab), table(tt), names(table(tt)) )
    colnames(sif_tabF) = c("asv", "type", "srcs", "nb_lien", "rownames")

    length(unique(sif_tabF[,1]))
    length(taxa_names(data1))

    sif_tabF <- as.data.frame(sif_tabF)

    # Metriques links
    #nombres de liens pour chaque OTU
    table(sif_tabF[,"asv"])
    table(sif_tabF[,"srcs"])


    write.table(as.data.frame(sif_tabF), glue("{output}/sif_tab2.tsv"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

    # Node table output
    src = unique(sif_tabF[,3])
    trgt = sif_tabF[,1]
    freq = taxa_sums(data1) / sum(taxa_sums(data1))

    # Taille des noeud en fonction de la frequence, noeud "groupe" 0.75 de la freq max.
    node_table=cbind.data.frame(c(as.character(src), as.character(trgt)),
      c(rep("source", length(src)), rep("cible", length(trgt))),
      c(rep(NA, length(src)), freq[trgt]) )
    names(node_table)=c("names","attribute", "freq")

    node_table[node_table$attribute=="source","freq"] = 0.75*max(na.omit(node_table$freq))

    node_table$tax = "source"
    node_table[node_table$attribute == "cible","tax"] = tax_table(data1)[as.character(node_table[node_table$attribute == "cible",1]),"Phylum"]

    write.table(as.data.frame(node_table), glue("{output}/node_tab2.tsv"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

  }


}
