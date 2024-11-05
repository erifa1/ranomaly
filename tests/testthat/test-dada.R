test_that("dada & co", {
  temp_file <- tempfile()

    chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    print(chk)
    if (nzchar(chk) && chk == "TRUE") {
      # use 2 cores in CRAN/Travis/AppVeyor
      num_workers <- 2L
    } else {
      # use all cores in devtools::test()
      num_workers <- parallel::detectCores()
    }
    
    #dada
    input_dir = test_path("testdata", "fastq")
    output_dir = glue::glue("{temp_file}/analyses/01_results_dada2/")

    dada_res = dada2_fun(cutadapt = FALSE, trim_l = 18, trim_r = 22, path = input_dir, outpath = output_dir, plot = FALSE, compress = TRUE, verbose = 1, paired = TRUE, torrent_single = FALSE ,returnval = TRUE, extension = "_R1.fastq.gz", n_cpu= num_workers)

    expect_true(class(dada_res) == "list")
    expect_true(dir.exists(output_dir))

    # assign 
    download.file("https://nextcloud.inrae.fr/s/YHi3fmDdEJt5cqR/download?path=%2F&files=GTDB_bac120_arc122_IDTAXA.rdata", destfile= glue::glue("{temp_file}/GTDB_bac120_arc122_IDTAXA.rdata"))

    output_dir = glue::glue("{temp_file}/analyses/02_idtaxa/")
    tax.tablecheck = assign_taxo_fun(dada_res = dada_res, id_db = glue::glue("{temp_file}/GTDB_bac120_arc122_IDTAXA.rdata"), verbose = 1 , output = output_dir)
    
    expect_true(class(tax.tablecheck) == "data.frame")
    expect_true(file.exists(glue::glue("{output_dir}/robjects.Rdata") ))

    # tree
    tree = generate_tree_fun(dada_res)
    expect_true(class(tree) == "phylo")


    # phyloseq
    data = generate_phyloseq_fun(dada_res = dada_res, tax.table = tax.tablecheck, tree = tree, metadata = glue::glue("{test_path('testdata', 'csv2phyloseq')}/sample_metadata.csv"))
    expect_true(class(data) == "phyloseq")

#   unlink(test_path("tmp"), recursive = TRUE, force = TRUE)
})