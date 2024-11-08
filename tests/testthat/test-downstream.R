test_that("downstream", {
  temp_file <- tempfile()
    
    load(glue::glue("{test_path('testdata')}/toy_physeq.Rdata"))
    expect_true(class(data) == "phyloseq")

    data_decontam = decontam_fun(data = data, domain = "Bacteria", column = "type", ctrl_identifier = "control", spl_identifier = "sample", number = 100, krona= FALSE, output = glue::glue("{temp_file}/analyses/decontam_out/"))
    expect_true(class(data_decontam) == "phyloseq")


    #rarefaction

    rare1 <- rarefaction(data)
    expect_true(class(rare1)[1] == "plotly")

    # plot 
    pabs = bars_fun(data = data, top = 10, Ord1 = "strain_time", rank="Genus", relative = FALSE, outfile = glue::glue("{temp_file}/analyses/05_all_res/plot_compo_rel.html"))
    expect_true(class(pabs)[1] == "plotly")


    # alpha 
    alpha_strain_time = diversity_alpha_fun( data, column1 = "strain_time", output = glue::glue("{temp_file}/analyses/05_all_res/plot_div_alpha_strain_time/"))
    expect_true(class(alpha_strain_time) == "list")


    # beta
    beta_strain_time = diversity_beta_light(data, col = "strain_time", dist0 = "bray", ord0 = "MDS", output=glue::glue("{temp_file}/analyses/05_all_res/plot_div_beta_strain_time/", tests = TRUE))
    expect_true(class(beta_strain_time) == "list")

    # diff

    metacoder1 = metacoder_fun(data = data, output = glue::glue("{temp_file}/analyses/05_all_res/diff/metacoder"), column1 = "strain_time", rank = "Family", signif = TRUE, plottrees = TRUE, min ="10", comp = "wildtype_t50~mutant_t50,wildtype_t0~mutant_t0")

    expect_true(class(metacoder1) == "list")

    deseq1 = deseq2_fun(data = data, output = glue::glue("{temp_file}/analyses/05_all_res/diff/deseq"), column1 = "strain_time", verbose = 1, rank = "Family", comp = "wildtype_t50~mutant_t50,wildtype_t0~mutant_t0")

    expect_true(class(deseq1) == "list")

    metaG1 = metagenomeseq_fun(data = data, output = glue::glue("{temp_file}/analyses/05_all_res/diff/metagenomeseq"), column1 = "strain_time", verbose = 1, rank = "Family", comp = "wildtype_t50~mutant_t50,wildtype_t0~mutant_t0")

    expect_true(class(metaG1) == "list")

    resF = aggregate_fun(data = data, metacoder = glue::glue("{temp_file}/analyses/05_all_res/diff/metacoder/metacoder_signif_Family.csv"), deseq = glue::glue("{temp_file}/analyses/05_all_res/diff/deseq"), mgseq = glue::glue("{temp_file}/analyses/05_all_res/diff/metagenomeseq"), output = glue::glue("{temp_file}/analyses/05_all_res/diff/aggregate_diff/"), column1 = "strain_time", column2 = NULL, verbose = 1, rank = "Family", comp = "wildtype_t50~mutant_t50,wildtype_t0~mutant_t0")

    expect_true(class(resF) == "list")

    #PLSDA 

    plsda_res <- plsda_fun(data = data, output = glue::glue("{temp_file}/analyses/05_all_res/plsda_family/"), column1 = "strain_time", rank = "Family")

    expect_true(class(plsda_res) == "list")

    # Others 

    heatmap_plot = heatmap_fun(data = data, column1 = "strain_time", top = 20, output = glue::glue("{temp_file}/analyses/05_all_res/plot_heatmap/"), rank = "Species")

    expect_true(class(heatmap_plot) == "list")

    outvenn = ASVenn_fun(data = data, output = glue::glue("{temp_file}/analyses/05_all_res/ASVenn/"), rank = "ASV", column1 = "strain_time", shared = TRUE)
    expect_true(class(outvenn) == "list")

    phy2tsv_fun(data = data, output = glue::glue("{temp_file}/analyses/05_all_res/tsv_table/"), rank = "ASV", relative = FALSE)

    expect_true(all(file.exists(glue::glue("{temp_file}/analyses/05_all_res/tsv_table/metadata_table.csv"), glue::glue("{temp_file}/analyses/05_all_res/tsv_table/otu_table_ASV.csv"))))
    

#   unlink(test_path("tmp"), recursive = TRUE, force = TRUE)
    unlink(testthat::test_path("tree"), recursive = TRUE, force = TRUE)
    unlink(testthat::test_path("phyloseq"), recursive = TRUE, force = TRUE)
    unlink(testthat::test_path("_snaps"), recursive = TRUE, force = TRUE)
    unlink(testthat::test_path("treealignment.Rdata"), recursive = TRUE, force = TRUE)
})