suppressMessages(library(futile.logger))
suppressMessages(library("optparse"))

option_list = list(
  make_option(c("-r", "--rdata"), type="character", default=NULL,
              help="Rdata file path.", metavar="path"),
  make_option(c("-o", "--out"), type="character", default="./idtaxa/",
              help="output .Rdata file name [default= %default]", metavar="path"),
  # make_option(c("-d", "--r_db"), type="character", default="/home/db/unite/sh_general_release_dynamic_01.12.2017.fasta",
  #             help="RDP database fasta file path [default= %default]", metavar="fasta"),
  # make_option(c("-a", "--assign"), type="character", default="idtaxa",
  #             help="Assignation program to use. [default= %default]", metavar="character"),
  make_option(c("-i", "--id_db"), type="character", default="/home/db/unite/UNITE_idtaxa201706.Rdata",
              help="IDTAXA database rdata file path [default= %default]. For double database, separate by comma.", metavar="rdata"),
  make_option(c("-v", "--verbose"), type="integer", default=1,
              help="Verbose level. [default= %default]", metavar="character"),
  make_option(c("-c", "--confidence"), default="50", , type="integer",
              help="Bootstrap threshold 0...100 [default %default%]"),
  make_option(c("-x","--returnval"), type="logical", default=FALSE,
              help='Return results list in terminal, set to FALSE in bash mode, results are outputted in ./dada_out/robjects.Rdata. [default= %default]')
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if(opt$verbose == 3){
  invisible(flog.threshold(DEBUG))
} else {
  invisible(flog.threshold(INFO))
}

if (is.null(opt$rdata)){
  print_help(opt_parser)
  stop("You must provide the Rdata file path.", call.=FALSE)
} else{
  flog.info("Loading data...")
  load(opt$rdata)
  flog.info('Done.')
}

print(opt$returnval);

suppressMessages(library(ranomaly))

# need dada2 output named dada_res in rdata
assign_taxo(dada_res = dada_res,  out = opt$out, id_db = opt$id_db, verbose = opt$verbose, confidence = opt$confidence, returnval = opt$returnval)



