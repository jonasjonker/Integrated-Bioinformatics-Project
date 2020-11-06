# install and load ArchR
devtools::install_github("GreenleafLab/ArchR", ref="release_1.0.0", repos = BiocManager::repositories()) # release branch
# devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories()) # master
library(ArchR)

# install archr dependencies not installed by default
ArchR::installExtraPackages()

# set number of threads (usually set to half the number of available cores)
addArchRThreads(threads = 1)

# add genome
addArchRGenome("hg38")

# assign paths to input fragment files
unsorted_3k_fragment <- "/home/james/Documents/leuven/second-year/IBP/ArchR/data/3k_unsorted.tsv.gz"
sorted_3k_fragment <- "/home/james/Documents/leuven/second-year/IBP/ArchR/data/3k_sorted.tsv.gz"
unsorted_10k_fragment <- "/home/james/Documents/leuven/second-year/IBP/ArchR/data/10k_unsorted.tsv.gz"
sorted_10k_fragment <- "/home/james/Documents/leuven/second-year/IBP/ArchR/data/10k_sorted.tsv.gz"

# assigning files to character vector for arrowfile creation
input_files <- c(unsorted_3k_fragment, sorted_3k_fragment, unsorted_10k_fragment, sorted_10k_fragment)
sample_names <- c("3k_unsorted", "3k_sorted", "10k_unsorted", "10k_sorted")
output_names <- c("3k_unsorted", "3k_sorted", "10k_unsorted", "10k_sorted")

# creating arrowfiles
ArrowFiles <- createArrowFiles(
  inputFiles = input_files,
  sampleNames = sample_names,
  outputNames = output_names,
  filterTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
