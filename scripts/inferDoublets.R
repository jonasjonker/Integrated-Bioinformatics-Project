'Find Doublets.

Usage:
  inferDoublets.R <data_dir>

Options:
  -h --help     Show this screen.
  --version     Show version.
' -> doc
library(docopt)
args <- docopt(doc, version = '1.0')

# preconditions
if (!dir.exists(args$data_dir)) stop(args$data_dir, " doesn't exist.")

# load libraries
library(ArchR)

# run parameters
addArchRThreads(threads = 12)
addArchRGenome("hg38")

# file locations
all_files  <- list.files(args$data_dir)
ArrowFiles <- file.path(args$data_dir, all_files[endsWith(all_files, ".arrow")])

# find doublets
doubScores <- addDoubletScores(
  input     = ArrowFiles,
  k         = 10,      # count k cells near a "pseudo-doublet".
  knnMethod = "UMAP",  # embedding for nearest neighbor search.
  LSIMethod = 1
)
