# Rscript ANOMALY pipeline:
setwd("~/Repository/LRF/00_erifa_bak/ranomaly_test")

dada_res = dada2_fun(path="~/Repository/LRF/00_erifa_bak/ranomaly_test/reads", compress=TRUE, plot=TRUE)

tax.table = assign_taxo(dada_res = dada_res, id_db = "/home/erifa/bank/silva/SILVA_SSU_r132_March2018.RData")
