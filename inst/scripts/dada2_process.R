# DADA2 bash launcher

library(futile.logger)
suppressMessages(library("optparse"))

option_list = list(
  make_option(c("-a", "--amplicon"), type="character", default="16S",
              help="16S or ITS, options change according to amplicon.", metavar="character"),
  make_option(c("-p", "--path"), type="character", default=NULL,
              help="Read files folder path", metavar="path"),
  make_option(c("-o", "--out"), type="character", default="./dada2_out/",
              help="output .Rdata file name [default= %default]", metavar="path"),
  make_option(c("-f", "--f_trunclen"), type="integer", default=240,
              help="Forward read tuncate length. [default= %default]", metavar="character"),
  make_option(c("-r", "--r_trunclen"), type="integer", default=240,
              help="Reverse read tuncate length. [default= %default]", metavar="character"),
  make_option(c("-5", "--f_primer"), type="character", default="GCATCGATGAAGAACGCAGC",
              help="Forward primer sequence [default= %default] (Only for ITS)", metavar="character"),
  make_option(c("-3", "--r_primer"), type="character", default="TCCTCCGCTTWTTGWTWTGC",
              help="Reverse primer sequence [default= %default] (Only for ITS)", metavar="character"),
  make_option(c("-q","--plot"), type="logical", default=FALSE,
              help='Plot all test. [default = %default]'),
  make_option(c("-c","--compress"), type="logical", default=FALSE,
              help='Files compress (.gz). [default= %default]'),
  make_option(c("-v", "--verbose"), type="integer", default=1,
              help="Verbose level. [default= %default]", metavar="character"),
  make_option(c("-x","--returnval"), type="logical", default=FALSE,
              help='Return results list in terminal, set to FALSE in bash mode, results are outputted in ./dada_out/robjects.Rdata. [default= %default]')
);

opt_parser = OptionParser(option_list=option_list, add_help_option = TRUE);
opt = parse_args(opt_parser);


if (is.null(opt$amplicon)){
  stop("Please choose which amplicon to analyse with -a option (16S or ITS).", call.=FALSE)
}

if (is.null(opt$path) & opt$amplicon=="16S"){
  print_help(opt_parser)
  stop("You must provide the reads files folder path.", call.=FALSE)
}

suppressMessages(library(ranomaly))

dada2_fun(opt$amplicon, path = opt$path, outpath = opt$out, f_trunclen = opt$f_trunclen, r_trunclen = opt$r_trunclen,
                    f_primer = opt$f_primer, r_primer = opt$r_primer, plot = opt$plot, compress = opt$compress, verbose = opt$verbose, returnval = opt$returnval)



