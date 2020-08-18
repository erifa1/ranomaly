
devtools::load_all("~/home-local-ssd/repository/ranomaly/")

dada_res = dada2_fun(path="~/home-local-ssd/projets/anomaly/reads_celine_anomaly", dadapool = "pseudo", compress=TRUE, plot=FALSE)

tax.table = assign_taxo_fun(dada_res = dada_res, id_db = "~/home-local-ssd/bank/SILVA_SSU_r132_March2018.RData" )

tree = generate_tree_fun(dada_res)

data = generate_phyloseq_fun(dada_res = dada_res, taxtable = tax.table, tree = tree, metadata = "~/home-local-ssd/projets/anomaly/sample_metadata.csv")

data

data_decontam = decontam_fun(data = data, domain = "Bacteria", number = 100, prev = 2, freq = 5e-05, skip =FALSE,
  column = "type", ctrl_identifier = "control", spl_identifier = "sample")
#column = "type", ctrl_identifier = "control", spl_identifier = "sample"

phy2tsv_fun(data = data_decontam)

# avec témoin sans decontam
bars_fun2(data = data, Ord1 = "souche_temps", Fact1 = "souche_temps", rank="Genus", top = 20)

bars_fun2(data = data_decontam, Ord1 = "souche_temps", Fact1 = "souche_temps", rank="Genus", top = 20)



alpha1 = diversity_alpha_fun(data = data_decontam, output = "./plot_div_alpha/", column1 = "souche", column2 = "temps",
                    column3 = "", supcovs = "", measures = "Observed" )
plotly::ggplotly(alpha1$plot) %>%
  layout(boxmode = "group")


beta1 = diversity_beta_light(psobj = data_decontam, rank = "ASV", col = "souche_temps", cov=NULL, dist0 = "bray", ord = "MDS")
plotly::ggplotly(beta1$plot)


### ICI
out1 = metacoder_fun(data = data, output = "./metacoder", column1 = "souche_temps", rank = "Family",
                          signif = TRUE, plottrees = TRUE, min ="10", comp = "sauvage_t50~mutant_t50,sauvage_t0~mutant_t0")

out2 = deseq2_fun(data = data, output = "./deseq/", column1 = "souche_temps", verbose = 1, rank = "Family", comp = "sauvage_t50~mutant_t50,sauvage_t0~mutant_t0")


out3 = metagenomeseq_fun(data = data, output = "./metagenomeseq/", column1 = "souche_temps", verbose = 1, rank = "Family", comp = "sauvage_t50~mutant_t50,sauvage_t0~mutant_t0")

outF = aggregate_fun(data = data, metacoder = "./metacoder/metacoder_signif_Family.csv", deseq = "./deseq/", mgseq = "./metagenomeseq/", output = "./aggregate_diff/",
                          column1 = "souche_temps", column2 = NULL, verbose = 1, rank = "Genus", comp = "sauvage_t50~mutant_t50,sauvage_t0~mutant_t0")
head(outF)
