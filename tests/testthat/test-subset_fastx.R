
test_that("multiplication works", {
  
  output_dir = test_path("testdata", "fastq")

  subset_fastx(path = output_dir, format = "fastq", outformat = "fastq", output = "./subset_fastq50/", nbseq = 50, ncores = 10, compress=TRUE, verbose=FALSE, random = FALSE, seed = NULL)

  expect_true(dir.exists(output_dir))

  fastq_files <- list.files("./subset_fastq50/", pattern = "\\.fastq.gz$", full.names = TRUE)
  expect_true(length(fastq_files) > 0) 
})
