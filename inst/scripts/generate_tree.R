library(futile.logger)
suppressMessages(library("optparse"))

option_list = list(
  make_option(c("-r", "--rdata"), type="character", default=NULL,
              help="Rdata file path.", metavar="path"),
  make_option(c("-o", "--out"), type="character", default="./tree",
              help="output .Rdata file name [default= %default]", metavar="path")
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$rdata)){
  print_help(opt_parser)
  stop("You must provide the Rdata file path.", call.=FALSE)
} else{
  flog.info("Loading data...")
  load(opt$rdata)
  flog.info('Done.')
}

suppressMessages(library(ranomaly))

# robjects need to contain dada2 output named dada_res.
generate_tree_fun(dada_res)
