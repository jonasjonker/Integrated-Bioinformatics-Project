###############################################################################
# Create arrowfiles from ... files.
#
# To install dependancies run:
#     devtools::install_github("GreenleafLab/ArchR", ref="release_1.0.0", repos = BiocManager::repositories())
#     ArchR::installExtraPackages()
# 
# Or use docker:
#     jonasjonker/archr-notebook:latest
#
# Data used can be found at:
#     https://support.10xgenomics.com/single-cell-multiome-atac-gex/datasets
###############################################################################

# dataset locations
input_files  <- c("../data/pbmc_unsorted_3k_atac_fragments.tsv.gz", 
                  "../data/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz",
                  "../data/pbmc_unsorted_10k_atac_fragments.tsv.gz",
                  "../data/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz")

# preconditions
if (!dir.exists("../data")) {
    stop("../data doesn't exist. Create directory or change working directory.")
}
if (!all(file.exists(input_files))) {
    stop("At least one of the specified files doesn't exist.")
}

# load libraries
library(ArchR)

# run parameters
addArchRThreads(threads = 10)
addArchRGenome("hg38")

# input vectors 
sample_names <- c("3k_unsorted",
                  "3k_sorted",
                  "10k_unsorted",
                  "10k_sorted")
output_names <- c("../data/3k_unsorted",
                  "../data/3k_sorted",
                  "../data/10k_unsorted",
                  "../data/10k_sorted")

# create arrowfiles
ArrowFiles <- createArrowFiles(
  inputFiles      = input_files,
  sampleNames     = sample_names,
  outputNames     = output_names,
  filterTSS       = 4,            # Don't set this too high initially
  filterFrags     = 1000, 
  addTileMat      = TRUE,
  addGeneScoreMat = TRUE,
  subThreading    = FALSE         # ArchR doesn't lock hdf5 file properly 
                                  # which causes problems when recombining
                                  # tmp hdf5 blocks.
                                  # https://github.com/GreenleafLab/ArchR/issues/248
)
