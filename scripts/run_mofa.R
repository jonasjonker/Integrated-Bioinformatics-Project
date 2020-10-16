args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
    stop("This script need an input and output file as argument.")
}

if (! file.exists(args[1])) {
    stop("Input file does not exist.")
}

input_file <- args[1]
output_file <- args[2]

untrained_mofa_obj <- readRDS(input_file)

run_mofa(untrained_mofa_obj, outfile=output_file)

library(MOFA2)
