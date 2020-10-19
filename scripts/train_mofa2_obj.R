args <- commandArgs(trailingOnly = TRUE)

infile <- args[1]
outfile <- args[2]

if (length(args) != 2) {
    stop("This script need an input and output file as argument.")
}

if (! file.exists(infile)) {
    stop("Input file does not exist.")
}

untrained_mofa_obj <- readRDS(input_file)

library(MOFA2)
library(reticulate)

py_config()

run_mofa(untrained_mofa_obj, outfile=output_file)
