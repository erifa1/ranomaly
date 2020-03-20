# Rscript ANOMALY pipeline:
setwd("~/Repository/LRF/00_erifa_bak/ranomaly_test")

dada_res = dada2_fun(path="~/Repository/LRF/00_erifa_bak/ranomaly_test/reads", compress=TRUE, plot=TRUE)

tax.table = assign_taxo(dada_res = dada_res, id_db = "/home/erifa/bank/silva/SILVA_SSU_r132_March2018.RData")

tree = generate_tree_fun(dada_res)


data = generate_phyloseq_fun(otutable = dada_res, taxtable = tax.table, tree = tree, metadata = "./sample-metadata.csv")


data = decontam_fun(data = data, domain = TRUE, column = "type", ctrl_identifier = "control", spl_identifier = "sample", number = 100)

export_to_stamp_fun(data = data)

phy2tsv_fun(data = data)

split_table_fun(data, column1 = "lot")

save.image("test_image.rdata")


bars_fun(data = data, column1 = "temps", column2 = "lot", rank=c("Phylum,Genus"), rare = "lot")


diversity_alpha_fun(data = data, output = "./plot_div_alpha/", column1 = "temps", column2 = "lot",
                    column3 = "", supcovs = "", measures = c("Observed","Shannon","Simpson","InvSimpson") )


diversity_beta_fun(data = data, output = "./plot_div_beta/", glom = "ASV", column1 = "temps", column2 = "lot", covar ="")
