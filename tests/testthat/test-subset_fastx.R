
test_that("utils", {
  temp_file <- tempfile()
  input_dir = test_path("testdata", "fastq")

  subset_fastx(path = input_dir, format = "fastq", outformat = "fastq", output = glue::glue("./{temp_file}/subset_fastq50/"), nbseq = 50, ncores = 10, compress=TRUE, verbose=FALSE, random = FALSE, seed = NULL)

  expect_true(dir.exists(input_dir))

  fastq_files <- list.files(glue::glue("./{temp_file}/subset_fastq50/"), pattern = "\\.fastq.gz$", full.names = TRUE)
  expect_true(length(fastq_files) > 0) 

  unlink(test_path("tmp"), recursive = TRUE, force = TRUE)
})
