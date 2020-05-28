# Rscript ANOMALY functions:

library(ranomaly)
#system.file("reads", "", package="ranomaly")  # Reads folder
#system.file("supdata", "sample-metadata.csv", package="ranomaly") # metadata

#setwd("")

dada_res = dada2_fun(path=system.file("reads", "", package="ranomaly"), dadapool = "pseudo", compress=TRUE, plot=FALSE)

tax.table = assign_taxo_fun(dada_res = dada_res, id_db = "~/bank/silva/SILVA_SSU_r132_March2018.RData" )

tree = generate_tree_fun(dada_res)

data = generate_phyloseq_fun(dada_res = dada_res, taxtable = tax.table, tree = tree, metadata = system.file("supdata", "sample-metadata.csv", package="ranomaly"))

data = decontam_fun(data = data, domain = "Bacteria", column = "type", ctrl_identifier = "control", spl_identifier = "sample", number = 100)

export_to_stamp_fun(data = data)

phy2tsv_fun(data = data)

split_table_fun(data, column1 = "lot")

save.image("test_image.rdata")

# load("test_image.rdata")

bars_fun(data = data, column1 = "temps", column2 = "lot", rank=c("Phylum,Genus"), rare = "lot")

diversity_alpha_fun(data = data, output = "./plot_div_alpha/", column1 = "temps", column2 = "lot",
                    column3 = "", supcovs = "", measures = c("Observed","Shannon","Simpson","InvSimpson") )

diversity_beta_fun(data = data, output = "./plot_div_beta/", glom = "ASV", column1 = "temps", column2 = "lot", covar ="")


diversity_beta_light(psobj = data, rank = "ASV", col = "temps", cov="lot", dist0 = "unifrac", ord = "MDS")


metacoder_fun(data = data, output = "./metacoder", column1 = "temps_lot", column2 = "", rank = "Genus",
                          signif = TRUE, plottrees = TRUE, min ="10", comp = "T6_lot1~T6_lot3,T9_lot1~T9_lot3")

deseq2_fun(data = data, output = "./deseq/", column1 = "temps_lot", verbose = 1, rank = "Genus", comp = "T6_lot1~T6_lot3,T9_lot1~T9_lot3")


metagenomeseq_fun(data = data, output = "./metagenomeseq/", column1 = "temps_lot", verbose = 1, rank = "Genus", comp = "T6_lot1~T6_lot3,T9_lot1~T9_lot3")

TABF = aggregate_fun(data = data, metacoder = "./metacoder/metacoder_signif_Genus.csv", deseq = "./deseq/", mgseq = "./metagenomeseq/", output = "./aggregate_diff/",
                          column1 = "temps_lot", column2 = NULL, verbose = 1, rank = "Genus", comp = "T6_lot1~T6_lot3,T9_lot1~T9_lot3")
head(TABF)



#other TOOLS
ASVenn_fun(data = data, output = "./ASVenn/", rank = "ASV",
                            column1 = "lot", subset = "", lvls = "", krona = "",
                            shared = TRUE)

csv2phyloseq_fun(otutable = "otutable.csv", taxtable = "taxo.csv",
                             seq = "seq.csv", metadata = "metadata.csv", output = "./csv2phyloseq/")

update_metadata_fun(data = data, output = "./updated_physeq/", metadata = "./sample-metadata_up1.csv", )

phy2cyto_fun(data = data, output = "./cytoscape/", column1 = "lot", repl = NULL, verbose = 1)


ttax1 = idtaxa_assign_fun("seqs.fna", id_db = "SILVA_SSU_r132_March2018.RData")



# IDTAXA database formatting

taxtable <- read.table(system.file("supdata", "gtdb_1k.tax", package="ranomaly"), sep="\t", stringsAsFactors=FALSE, h = TRUE)
taxtable[taxtable == ""] = NA
row.names(taxtable) = taxtable[,1]	#rownames must be unique
ttable_ok = taxtable[,-1]

#Fill tax table
filltable = fill_tax_fun(ttable_ok)
#Check tax table
check1 = check_tax_fun(filltable, output = NULL)
#Generate taxid
taxid <- taxid_fun(taxtable = check1, output = NULL)

#Train IDTAXA db
idtaxa_traindb(taxtable = check1, taxid = taxid, seqs = system.file("supdata", "gtdb_1k.fna", package="ranomaly"), prunedb = NULL, outputDIR = "./", outputDBname = "newDB.rdata", returnval = FALSE)



