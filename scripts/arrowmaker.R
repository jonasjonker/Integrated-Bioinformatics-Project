'Make arrowfile(s) from fragment file(s).

Usage:
  arrowmaker.R <data_dir>

Options:
  -h --help     Show this screen.
  --version     Show version.
' -> doc
library(docopt)
args <- docopt(doc, version = 'arrowmaker 1.0')

# preconditions
if (!dir.exists(args$data_dir)) stop(args$data_dir, " doesn't exist.")

# load libraries
library(ArchR)
library(tools)

# run parameters
addArchRThreads(threads = 10)
addArchRGenome("hg38")
 
# file locations
all_files    <- list.files(args$data_dir) 
input_files  <- all_files[endsWith(all_files, ".tsv.gz")]
input_paths  <- file.path(args$data_dir, input_files)
sample_names <- file_path_sans_ext(input_files, compression = TRUE)
output_paths <- file.path(args$data_dir, "ArrowFiles", sample_names)

# create arrowfiles
ArrowFiles <- createArrowFiles(
  inputFiles      = input_paths,
  sampleNames     = sample_names,
  outputNames     = output_paths,
  filterTSS       = 4,            # minimum required Transcription start site (TSS) enrichment score 
                                  # Don't set this too high initially
  filterFrags     = 1000, 
  addTileMat      = TRUE,
  addGeneScoreMat = TRUE,
  subThreading    = FALSE         # ArchR doesn't lock hdf5 file properly 
                                  # which causes problems when recombining
                                  # tmp hdf5 blocks.
                                  # https://github.com/GreenleafLab/ArchR/issues/248
)
