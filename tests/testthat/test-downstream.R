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

#   unlink(test_path("tmp"), recursive = TRUE, force = TRUE)
    unlink(testthat::test_path("tree"), recursive = TRUE, force = TRUE)
    unlink(testthat::test_path("phyloseq"), recursive = TRUE, force = TRUE)
    unlink(testthat::test_path("_snaps"), recursive = TRUE, force = TRUE)
    unlink(testthat::test_path("treealignment.Rdata"), recursive = TRUE, force = TRUE)
})